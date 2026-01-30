from scipy.integrate import solve_ivp
import numpy as np
import sk_dsp_comm.sigsys as ss
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider
import matplotlib.animation as animation
from scipy.interpolate import griddata

# ==================== SIMULATION CODE ====================

def sigmoid(u):
    return np.tanh(u)

def N_oscillators(t, y, N, h_ex_rand, h_in_rand, 
                        coupling_matrix_EE, coupling_matrix_EI, cconst_E, cconst_I, pars, sr, time_stop, pert, pert_osc_list):
        
    tau_ex, tau_in, c2, c4 = pars
    
    time_index = int(t*sr)
    
    if time_index >= time_stop*sr:
        dydt = np.zeros(2*N)
        return dydt

    # Separate Variables
    y_ex = y[:-1:2]
    y_in = y[1::2]
    dy_ex, dy_in = np.zeros(N), np.zeros(N)
    dydt = np.zeros(2*N)

    for osc in np.arange(N):
        coup_EE = cconst_E*sigmoid(sum(coupling_matrix_EE[:, osc] * y_ex))
        coup_EI = cconst_I*sigmoid(sum(coupling_matrix_EI[:, osc] * y_ex))

        if osc in pert_osc_list:
            dy_ex[osc] = (pert[time_index] - y_ex[osc] - c2*sigmoid(y_in[osc]) + coup_EE)*tau_ex 
            dy_in[osc] = (h_in_rand[osc]   - y_in[osc] - c4*sigmoid(y_in[osc]) + coup_EI)*tau_in
        else:
            dy_ex[osc] = (h_ex_rand[osc]   - y_ex[osc] - c2*sigmoid(y_in[osc]) + coup_EE)*tau_ex 
            dy_in[osc] = (h_in_rand[osc]   - y_in[osc] - c4*sigmoid(y_in[osc]) + coup_EI)*tau_in
            
    # Combine Variables
    dydt[:-1:2] = dy_ex
    dydt[1::2] = dy_in
    
    return dydt

def get_hexagon_centers_sorted(layers):
    """
    Generate hexagon centers sorted by layer (distance from origin)
    Returns list of (x, y, q, r, s, distance) tuples
    """
    center_spacing = 1.0
    centers = []
    
    for q in range(-layers+1, layers):
        for r in range(-layers+1, layers):
            s = -q - r
            if max(abs(q), abs(r), abs(s)) <= layers-1:
                x = center_spacing * (3/2 * q)
                y = center_spacing * (np.sqrt(3)/2 * q + np.sqrt(3) * r)
                distance = max(abs(q), abs(r), abs(s))
                centers.append((x, y, q, r, s, distance))
    
    # Sort by distance (layer), then by angle
    centers_sorted = sorted(centers, key=lambda c: (c[5], np.arctan2(c[1], c[0])))
    return centers_sorted
  

def create_connectivity_matrix(layers):
    """
    Create connectivity/adjacency matrix for hexagonal lattice
    Returns NxN matrix where entry (i,j) = 1 if hexagons i and j are neighbors
    """
    centers = get_hexagon_centers_sorted(layers)
    N = len(centers)
    
    # Initialize connectivity matrix
    connectivity = np.zeros((N, N), dtype=int)
    
    # Create a mapping from (q, r, s) coordinates to index
    coord_to_index = {}
    for i, (x, y, q, r, s, dist) in enumerate(centers):
        coord_to_index[(q, r, s)] = i
    
    # In hexagonal/cube coordinates, the 6 neighbors of (q, r, s) are:
    neighbor_offsets = [
        (1, -1, 0),   # East
        (1, 0, -1),   # Northeast
        (0, 1, -1),   # Northwest
        (-1, 1, 0),   # West
        (-1, 0, 1),   # Southwest
        (0, -1, 1)    # Southeast
    ]
    
    # For each hexagon, check its 6 potential neighbors
    for i, (x, y, q, r, s, dist) in enumerate(centers):
        for dq, dr, ds in neighbor_offsets:
            neighbor_coords = (q + dq, r + dr, s + ds)
            
            # Check if neighbor exists in our lattice
            if neighbor_coords in coord_to_index:
                j = coord_to_index[neighbor_coords]
                connectivity[i, j] = 1
    
    return connectivity, centers

# ==================== VISUALIZATION CODE ====================

class InteractiveHexagonalLattice:
    """
    Interactive hexagonal lattice with play/pause and slider controls
    Supports smooth interpolation between hexagons
    """
    def __init__(self, layers, data_array, time_array=None, cmap='viridis', vmin=None, vmax=None, fps=10, 
                 interpolate=False, show_edges=True, perturbed_osc=None):
        """
        Parameters:
        -----------
        layers : int
            Number of layers in the hexagonal lattice
        data_array : ndarray, shape (N, samples)
            Array with values for each hexagon at each time point
        time_array : ndarray, optional
            Time values corresponding to each sample
        cmap : str or Colormap
            Colormap to use (default: 'viridis')
        vmin, vmax : float, optional
            Min and max values for colormap normalization
        fps : int
            Frames per second for animation (default: 10)
        interpolate : bool
            If True, use smooth interpolation between hexagon centers (default: False)
        show_edges : bool
            If True, show hexagon edges (only visible when interpolate=False) (default: True)
        perturbed_osc : int or list, optional
            Index/indices of perturbed oscillator(s) to highlight
        """
        self.layers = layers
        self.data_array = data_array
        self.N, self.samples = data_array.shape
        self.time_array = time_array if time_array is not None else np.arange(self.samples)
        self.cmap = cmap
        self.current_frame = 0
        self.is_playing = False
        self.fps = fps
        self.anim = None
        self.interpolate = interpolate
        self.show_edges = show_edges
        self.perturbed_osc = perturbed_osc if perturbed_osc is not None else []
        if isinstance(self.perturbed_osc, int):
            self.perturbed_osc = [self.perturbed_osc]
        
        # Set colormap range
        self.vmin = vmin if vmin is not None else np.min(data_array)
        self.vmax = vmax if vmax is not None else np.max(data_array)
        
        # Get hexagon centers
        self.centers = get_hexagon_centers_sorted(layers)
        
        # Extract center coordinates
        self.center_coords = np.array([(x, y) for x, y, q, r, s, dist in self.centers])
        
        # Setup figure
        self.setup_figure()
        
    def setup_figure(self):
        """Initialize the figure and hexagons"""
        # Create figure with extra space for controls
        self.fig = plt.figure(figsize=(10, 10))
        
        # Main plot area (leave space at bottom for controls)
        self.ax = plt.axes([0.1, 0.25, 0.8, 0.7])
        
        self.center_spacing = 1.0
        self.angles = np.linspace(0, 2*np.pi, 7)
        
        # Normalize data to [0, 1] for colormap
        self.norm = plt.Normalize(vmin=self.vmin, vmax=self.vmax)
        self.colormap = plt.get_cmap(self.cmap)
        
        max_coord = self.layers * 2
        self.ax.set_xlim(-max_coord, max_coord)
        self.ax.set_ylim(-max_coord, max_coord)
        self.ax.set_aspect('equal')
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        
        if self.interpolate:
            # Create interpolated image
            self.setup_interpolated_display()
        else:
            # Create hexagon patches
            self.setup_hexagon_patches()
        
        # Add perturbed oscillator markers
        self.perturb_markers = []
        for osc_idx in self.perturbed_osc:
            if osc_idx < len(self.centers):
                x, y = self.center_coords[osc_idx]
                marker = self.ax.plot(x, y, 'r*', markersize=20, 
                                     markeredgecolor='white', markeredgewidth=1.5,
                                     label='Perturbed' if not self.perturb_markers else '')[0]
                self.perturb_markers.append(marker)
        
        if self.perturb_markers:
            self.ax.legend(loc='upper right')
        
        self.title = self.ax.set_title(f'Hexagonal Lattice - t={self.time_array[0]:.2f}s',
                                       fontsize=14, pad=15)
        
        # Add colorbar
        cbar_ax = plt.axes([0.92, 0.25, 0.02, 0.7])
        sm = plt.cm.ScalarMappable(cmap=self.colormap, norm=self.norm)
        sm.set_array([])
        self.cbar = plt.colorbar(sm, cax=cbar_ax)
        self.cbar.set_label('Excitatory Activity', rotation=270, labelpad=20)
        
        # Add controls
        self.setup_controls()
    
    def setup_hexagon_patches(self):
        """Create hexagon patches (non-interpolated mode)"""
        self.hexagons = []
        for i, (x, y, q, r, s, dist) in enumerate(self.centers):
            hex_x = x + np.cos(self.angles) * self.center_spacing
            hex_y = y + np.sin(self.angles) * self.center_spacing
            
            # Get initial color
            color = self.colormap(self.norm(self.data_array[i, 0]))
            
            edgecolor = 'black' if self.show_edges else 'none'
            linewidth = 1.5 if self.show_edges else 0
            patch = self.ax.fill(hex_x, hex_y, color=color, edgecolor=edgecolor, 
                               linewidth=linewidth)[0]
            self.hexagons.append(patch)
    
    def setup_interpolated_display(self):
        """Create interpolated image display"""
        # Create a regular grid for interpolation
        max_coord = self.layers * 2
        grid_resolution = 200  # Higher = smoother but slower
        
        x_grid = np.linspace(-max_coord, max_coord, grid_resolution)
        y_grid = np.linspace(-max_coord, max_coord, grid_resolution)
        self.X_grid, self.Y_grid = np.meshgrid(x_grid, y_grid)
        
        # Interpolate initial frame
        Z_grid = griddata(
            self.center_coords, 
            self.data_array[:, 0], 
            (self.X_grid, self.Y_grid), 
            method='cubic',
            fill_value=np.nan
        )
        
        # Create image
        self.image = self.ax.imshow(
            Z_grid, 
            extent=[-max_coord, max_coord, -max_coord, max_coord],
            origin='lower',
            cmap=self.colormap,
            vmin=self.vmin,
            vmax=self.vmax,
            interpolation='bilinear',
            aspect='equal'
        )
        
        # Optionally overlay hexagon edges
        if self.show_edges:
            for i, (x, y, q, r, s, dist) in enumerate(self.centers):
                hex_x = x + np.cos(self.angles) * self.center_spacing
                hex_y = y + np.sin(self.angles) * self.center_spacing
                self.ax.plot(hex_x, hex_y, 'k-', linewidth=0.5, alpha=0.3)
    
    def setup_controls(self):
        """Create interactive controls"""
        # Play/Pause button
        ax_play = plt.axes([0.15, 0.15, 0.1, 0.04])
        self.play_button = Button(ax_play, '▶ Play', color='lightgreen', hovercolor='green')
        self.play_button.on_clicked(self.toggle_play)
        
        # Reset button
        ax_reset = plt.axes([0.27, 0.15, 0.1, 0.04])
        self.reset_button = Button(ax_reset, '⟲ Reset', color='lightblue', hovercolor='blue')
        self.reset_button.on_clicked(self.reset)
        
        # Toggle interpolation button
        ax_toggle = plt.axes([0.39, 0.15, 0.15, 0.04])
        toggle_text = 'Smooth: ON' if self.interpolate else 'Smooth: OFF'
        self.toggle_button = Button(ax_toggle, toggle_text, color='lightyellow', hovercolor='yellow')
        self.toggle_button.on_clicked(self.toggle_interpolation)
        
        # Frame slider
        ax_slider = plt.axes([0.15, 0.09, 0.7, 0.03])
        self.slider = Slider(
            ax_slider, 
            'Frame', 
            0, 
            self.samples-1, 
            valinit=0, 
            valstep=1,
            color='steelblue'
        )
        self.slider.on_changed(self.on_slider_change)
        
        # FPS slider
        ax_fps = plt.axes([0.15, 0.03, 0.3, 0.03])
        self.fps_slider = Slider(
            ax_fps,
            'FPS',
            1,
            30,
            valinit=self.fps,
            valstep=1,
            color='coral'
        )
        self.fps_slider.on_changed(self.on_fps_change)
    
    def toggle_interpolation(self, event):
        """Toggle between interpolated and patch display"""
        self.interpolate = not self.interpolate
        
        # Clear current display - safely handle both modes
        if self.interpolate:
            # Switching TO interpolated mode - remove hexagon patches
            for patch in self.ax.patches[:]:
                patch.remove()
        else:
            # Switching TO hexagon mode - remove image if it exists
            if hasattr(self, 'image') and self.image in self.ax.images:
                self.image.remove()
        
        # Remove any lines that aren't perturb markers
        for line in self.ax.lines[:]:
            if line not in self.perturb_markers:
                line.remove()
        
        # Recreate display with new mode
        if self.interpolate:
            self.setup_interpolated_display()
            self.toggle_button.label.set_text('Smooth: ON')
        else:
            self.setup_hexagon_patches()
            self.toggle_button.label.set_text('Smooth: OFF')
        
        # Restore markers to front
        for marker in self.perturb_markers:
            marker.set_zorder(10)
        
        # Update to current frame
        self.update_frame(self.current_frame)
        plt.draw()
    
    def update_frame(self, frame):
        """Update hexagon colors for given frame"""
        if frame >= self.samples:
            frame = self.samples - 1
        if frame < 0:
            frame = 0
        
        if self.interpolate:
            # Update interpolated image
            Z_grid = griddata(
                self.center_coords, 
                self.data_array[:, frame], 
                (self.X_grid, self.Y_grid), 
                method='cubic',
                fill_value=np.nan
            )
            self.image.set_data(Z_grid)
        else:
            # Update hexagon patches
            for i, patch in enumerate(self.hexagons):
                color = self.colormap(self.norm(self.data_array[i, frame]))
                patch.set_facecolor(color)
        
        self.title.set_text(f'Hexagonal Lattice - t={self.time_array[frame]:.2f}s')
        self.current_frame = frame
        
        # Update slider without triggering callback
        self.slider.eventson = False
        self.slider.set_val(frame)
        self.slider.eventson = True
        
        if self.interpolate:
            return [self.image, self.title]
        else:
            return self.hexagons + [self.title]
    
    def toggle_play(self, event):
        """Toggle play/pause"""
        self.is_playing = not self.is_playing
        
        if self.is_playing:
            self.play_button.label.set_text('|| Pause')  # Using || instead of Unicode pause symbol
            self.play_button.color = 'yellow'
            self.start_animation()
        else:
            self.play_button.label.set_text('▶ Play')
            self.play_button.color = 'lightgreen'
            self.stop_animation()
    
    def start_animation(self):
        """Start the animation"""
        if self.anim is None:
            interval = 1000.0 / self.fps
            self.anim = animation.FuncAnimation(
                self.fig,
                self._animate_step,
                interval=interval,
                blit=False,
                cache_frame_data=False
            )
            plt.draw()
    
    def _animate_step(self, frame_number):
        """Internal animation step"""
        if self.is_playing:
            self.current_frame = (self.current_frame + 1) % self.samples
            return self.update_frame(self.current_frame)
        if self.interpolate:
            return [self.image, self.title]
        else:
            return self.hexagons + [self.title]
    
    def stop_animation(self):
        """Stop the animation"""
        if self.anim is not None:
            self.anim.event_source.stop()
            self.anim = None
    
    def on_slider_change(self, val):
        """Handle slider change"""
        if not self.is_playing:
            frame = int(val)
            self.update_frame(frame)
            plt.draw()
    
    def on_fps_change(self, val):
        """Handle FPS slider change"""
        self.fps = int(val)
        if self.is_playing:
            self.stop_animation()
            self.start_animation()
    
    def reset(self, event):
        """Reset to first frame"""
        was_playing = self.is_playing
        if self.is_playing:
            self.toggle_play(None)
        self.update_frame(0)
        plt.draw()
        if was_playing:
            self.toggle_play(None)
    
    def show(self):
        """Display the interactive plot"""
        plt.show()

# ==================== MAIN EXECUTION ====================

def run_simulation():
    """Run the hexagonal oscillator simulation"""
    print("=" * 60)
    print("HEXAGONAL OSCILLATOR NETWORK SIMULATION")
    print("=" * 60)
    
    # Number of hexagonal layers
    L = 4
    N = 3*L*(L-1)+1
    
    print(f"Number of layers: {L}")
    print(f"Number of oscillators: {N}")
    
    # Excitatory input parameter
    h_ex_0 = -6.56
    h_in_0 = -4.0
    
    eps = 0.01
    RANDOM_STATE = 11111
    np.random.seed(RANDOM_STATE)
    random_vals = eps*np.random.normal(0,1,size=N)
    random_vals_sorted = np.sort(random_vals)
    h_ex_rand = h_ex_0 - random_vals_sorted
    h_in_rand = h_in_0 - eps*np.random.normal(0,1,size=N)
    
    # Parameters
    pars = (1, 1.5, 10, 0)  # Homoclinic
    
    coupling_strength_EE, coupling_strength_EI = 5., 10.
    frac_EE, frac_EI = 0.1, 0.0
    
    # Coupling
    coupling_matrix_EE_ini, centers = create_connectivity_matrix(L)
    coupling_matrix_EI_ini = coupling_matrix_EE_ini.copy()
    
    coupling_matrix_EE = coupling_matrix_EE_ini*frac_EE
    coupling_matrix_EI = coupling_matrix_EI_ini*frac_EI
    
    np.fill_diagonal(coupling_matrix_EE, 1)
    np.fill_diagonal(coupling_matrix_EI, 1)
    
    # Time array
    time_stop = 100
    sr = 1000
    samples = time_stop*sr
    time = np.linspace(start=0, stop=time_stop, num=samples)
    
    pulse_fre = 0.23
    pulse_wid = 4.5
    pulse_amp = 10.0
    
    # Initial conditions
    y_ini = np.random.normal(size=2*N)
    y_ini = [-1.43744236, -12.9323564]*N
    
    pert = h_ex_0 + pulse_amp*ss.rect(np.mod(time, 1/pulse_fre)-(1/pulse_fre)/2+pulse_wid/2, pulse_wid)
    pert_osc = [0]
    
    print(f"\nSimulation parameters:")
    print(f"  Time: {time_stop}s")
    print(f"  Bifurcation: Homoclinic")
    print(f"  Coupling: E-E={coupling_strength_EE}, E-I={coupling_strength_EI}")
    print(f"  Perturbation: {pulse_fre} Hz, {pulse_wid}s width, amp={pulse_amp}")
    print(f"  Perturbed oscillator: {pert_osc[0]}")
    
    print(f"\nRunning simulation...")
    solution = solve_ivp(
        fun=N_oscillators, 
        t_span=(0, time_stop), 
        y0=y_ini, 
        t_eval=time,  
        args=(N, h_ex_rand, h_in_rand,  
              coupling_matrix_EE, 
              coupling_matrix_EI, 
              coupling_strength_EE, 
              coupling_strength_EI, 
              pars, sr, time_stop, pert, pert_osc), 
        method='BDF', 
        max_step=0.1
    )
    
    t, y_pert = solution.t, solution.y.T
    y_ex_only = y_pert[:,::2]
    
    print(f"Simulation complete!")
    print(f"Output shape: {y_ex_only.shape}")
    
    return L, y_ex_only, t, pert_osc

def main():
    """Main execution function"""
    # Run simulation
    L, y_ex_only, time_array, pert_osc = run_simulation()
    
    # Downsample for smoother interactive experience
    downsample = 10
    y_ex_downsampled = y_ex_only[::downsample, :].T  # Transpose to get (N_oscillators, N_samples)
    time_downsampled = time_array[::downsample]
    
    print(f"\nDownsampled data shape: {y_ex_downsampled.shape}")
    print(f"Creating interactive visualization...")
    print("\nControls:")
    print("  - Click '▶ Play' to start/pause animation")
    print("  - Click 'Smooth: ON/OFF' to toggle interpolation")
    print("  - Use 'Frame' slider to jump to specific time")
    print("  - Use 'FPS' slider to adjust animation speed")
    print("  - Click '⟲ Reset' to go back to t=0")
    
    # Create interactive visualization
    interactive = InteractiveHexagonalLattice(
        L, 
        y_ex_downsampled,
        time_array=time_downsampled,
        cmap='RdYlBu_r',
        vmin=np.percentile(y_ex_only, 1),
        vmax=np.percentile(y_ex_only, 99),
        fps=20,
        interpolate=True,   # Start with smooth interpolation
        show_edges=False,   # Hide edges for cleaner look
        perturbed_osc=pert_osc
    )
    
    interactive.show()

if __name__ == "__main__":
    main()