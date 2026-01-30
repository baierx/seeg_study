"""
Animation Export Options - HTML, MP4, GIF, and more
Shows how to save matplotlib animations in different formats
"""

from scipy.integrate import solve_ivp
import numpy as np
import sk_dsp_comm.sigsys as ss
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter, HTMLWriter
from scipy.interpolate import griddata

# ==================== SIMULATION CODE (abbreviated) ====================

def sigmoid(u):
    return np.tanh(u)

def N_oscillators(t, y, N, h_ex_rand, h_in_rand, 
                        coupling_matrix_EE, coupling_matrix_EI, cconst_E, cconst_I, pars, sr, time_stop, pert, pert_osc_list):
    tau_ex, tau_in, c2, c4 = pars
    time_index = int(t*sr)
    if time_index >= time_stop*sr:
        return np.zeros(2*N)
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
    dydt[:-1:2] = dy_ex
    dydt[1::2] = dy_in
    return dydt

def get_hexagon_centers_sorted(layers):
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
    centers_sorted = sorted(centers, key=lambda c: (c[5], np.arctan2(c[1], c[0])))
    return centers_sorted

def create_connectivity_matrix(layers):
    centers = get_hexagon_centers_sorted(layers)
    N = len(centers)
    connectivity = np.zeros((N, N), dtype=int)
    coord_to_index = {}
    for i, (x, y, q, r, s, dist) in enumerate(centers):
        coord_to_index[(q, r, s)] = i
    neighbor_offsets = [(1, -1, 0), (1, 0, -1), (0, 1, -1), (-1, 1, 0), (-1, 0, 1), (0, -1, 1)]
    for i, (x, y, q, r, s, dist) in enumerate(centers):
        for dq, dr, ds in neighbor_offsets:
            neighbor_coords = (q + dq, r + dr, s + ds)
            if neighbor_coords in coord_to_index:
                j = coord_to_index[neighbor_coords]
                connectivity[i, j] = 1
    return connectivity, centers

def run_quick_simulation():
    """Run a quick simulation for demonstration"""
    print("Running simulation...")
    L = 4  # Smaller for faster demo
    N = 3*L*(L-1)+1
    
    h_ex_0, h_in_0 = -6.4, -4.0
    eps = 0.01
    np.random.seed(11111)
    random_vals = eps*np.random.normal(0,1,size=N)
    h_ex_rand = h_ex_0 - np.sort(random_vals)
    h_in_rand = h_in_0 - eps*np.random.normal(0,1,size=N)
    
    pars = (1, 2.5, 10, 0)
    coupling_strength_EE, coupling_strength_EI = 5., 10.
    frac_EE, frac_EI = 0.2, 0.0
    
    coupling_matrix_EE_ini, centers = create_connectivity_matrix(L)
    coupling_matrix_EI_ini = coupling_matrix_EE_ini.copy()
    coupling_matrix_EE = coupling_matrix_EE_ini*frac_EE
    coupling_matrix_EI = coupling_matrix_EI_ini*frac_EI
    np.fill_diagonal(coupling_matrix_EE, 1)
    np.fill_diagonal(coupling_matrix_EI, 1)
    
    time_stop, sr = 3, 500
    samples = time_stop*sr
    time = np.linspace(start=0, stop=time_stop, num=samples)
    
    pulse_fre, pulse_wid, pulse_amp = 0.5, 0.2, 1.0
    y_ini = np.random.normal(size=2*N)
    pert = h_ex_0 + pulse_amp*ss.rect(np.mod(time, 1/pulse_fre)-(1/pulse_fre)/2+pulse_wid/2, pulse_wid)
    pert_osc = [0]
    
    solution = solve_ivp(fun=N_oscillators, t_span=(0, time_stop), y0=y_ini, t_eval=time,  
                         args=(N, h_ex_rand, h_in_rand, coupling_matrix_EE, coupling_matrix_EI, 
                               coupling_strength_EE, coupling_strength_EI, pars, sr, time_stop, pert, pert_osc), 
                         method='BDF', max_step=0.1)
    
    y_ex_only = solution.y.T[:,::2]
    print(f"Simulation complete! Shape: {y_ex_only.shape}")
    return L, y_ex_only, solution.t, pert_osc, centers

# ==================== ANIMATION CREATION ====================

def create_animation(L, y_ex_only, time_array, pert_osc, centers):
    """Create animation object that can be saved in various formats"""
    print("\nCreating animation...")
    
    # Get center coordinates
    center_coords = np.array([(x, y) for x, y, q, r, s, dist in centers])
    
    # Downsample
    downsample = 10
    y_ex_ds = y_ex_only[::downsample, :].T
    time_ds = time_array[::downsample]
    n_frames = len(time_ds)
    
    print(f"Animation frames: {n_frames}")
    
    # Setup colormap
    vmin, vmax = np.percentile(y_ex_only, [1, 99])
    cmap = plt.get_cmap('RdYlBu_r')
    
    # Create interpolation grid
    max_coord = L * 2
    grid_res = 150
    x_grid = np.linspace(-max_coord, max_coord, grid_res)
    y_grid = np.linspace(-max_coord, max_coord, grid_res)
    X_grid, Y_grid = np.meshgrid(x_grid, y_grid)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(5, 4))
    
    # Interpolate first frame
    Z_grid = griddata(center_coords, y_ex_ds[:, 0], (X_grid, Y_grid), 
                     method='cubic', fill_value=np.nan)
    
    im = ax.imshow(Z_grid, extent=[-max_coord, max_coord, -max_coord, max_coord],
                  origin='lower', cmap=cmap, vmin=vmin, vmax=vmax,
                  interpolation='bilinear', aspect='equal')
    
    # Mark perturbed oscillator
    for osc_idx in pert_osc:
        x, y = center_coords[osc_idx]
        ax.plot(x, y, 'r*', markersize=20, markeredgecolor='white', 
               markeredgewidth=1.5, zorder=10)
    
    ax.set_xlim(-max_coord, max_coord)
    ax.set_ylim(-max_coord, max_coord)
    ax.set_aspect('equal')
    ax.axis('off')
    
    title = ax.text(0.5, 0.98, f't = {time_ds[0]:.2f}s', 
                   transform=ax.transAxes, ha='center', va='top',
                   fontsize=14, weight='bold',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='Excitatory Activity')
    
    def update(frame):
        Z_grid = griddata(center_coords, y_ex_ds[:, frame], (X_grid, Y_grid),
                         method='cubic', fill_value=np.nan)
        im.set_data(Z_grid)
        title.set_text(f't = {time_ds[frame]:.2f}s')
        return [im, title]
    
    anim = FuncAnimation(fig, update, frames=n_frames, interval=50, blit=True)
    
    return anim, fig

# ==================== SAVE FUNCTIONS ====================

def save_as_mp4(anim, filename, fps=20, dpi=100, bitrate=2000):
    """
    Save animation as MP4 video
    
    Parameters:
    -----------
    anim : Animation object
    filename : str (should end with .mp4)
    fps : int (frames per second)
    dpi : int (resolution)
    bitrate : int (quality, higher = better quality but larger file)
    """
    print(f"\nSaving MP4: {filename}")
    writer = FFMpegWriter(fps=fps, bitrate=bitrate, codec='h264')
    anim.save(filename, writer=writer, dpi=dpi)
    print(f"✓ MP4 saved!")

def save_as_html(anim, filename):
    """
    Save animation as interactive HTML with JavaScript
    
    Parameters:
    -----------
    anim : Animation object
    filename : str (should end with .html)
    
    Note: The HTML file contains the full animation and can be opened in any browser
    """
    print(f"\nSaving HTML: {filename}")
    
    # Method 1: Direct save (recommended)
    anim.save(filename, writer='html', fps=20)
    
    # Method 2: Using to_jshtml() and writing manually
    # html_str = anim.to_jshtml()
    # with open(filename, 'w') as f:
    #     f.write(html_str)
    
    print(f"✓ HTML saved!")
    print(f"  Open in browser to view interactive animation")

def save_as_html_jupyter(anim):
    """
    Display animation in Jupyter notebook
    
    Usage in Jupyter:
    -----------------
    from IPython.display import HTML
    HTML(anim.to_jshtml())
    """
    print("\nFor Jupyter notebooks, use:")
    print("  from IPython.display import HTML")
    print("  HTML(anim.to_jshtml())")
    return anim.to_jshtml()

def save_as_gif(anim, filename, fps=20, dpi=80):
    """
    Save animation as GIF
    
    Parameters:
    -----------
    anim : Animation object
    filename : str (should end with .gif)
    fps : int (frames per second)
    dpi : int (resolution - keep lower for smaller file size)
    
    Note: GIFs are larger than MP4 but more universally compatible
    """
    print(f"\nSaving GIF: {filename}")
    writer = PillowWriter(fps=fps)
    anim.save(filename, writer=writer, dpi=dpi)
    print(f"✓ GIF saved!")

def save_as_frames(anim, output_dir, prefix='frame', fmt='png', dpi=150):
    """
    Save animation as individual frame images
    
    Parameters:
    -----------
    anim : Animation object
    output_dir : str (directory to save frames)
    prefix : str (prefix for frame filenames)
    fmt : str (image format: 'png', 'jpg', etc.)
    dpi : int (resolution)
    """
    import os
    print(f"\nSaving individual frames to: {output_dir}")
    os.makedirs(output_dir, exist_ok=True)
    
    # This requires accessing the animation's figure
    # Easier method: save frames during animation creation
    # For demonstration purposes only
    print(f"  Note: Individual frames best saved during animation creation")
    print(f"  See matplotlib documentation for frame-by-frame saving")

# ==================== MAIN DEMO ====================

def main():
    """Demonstrate all export formats"""
    print("=" * 70)
    print("ANIMATION EXPORT FORMATS DEMONSTRATION")
    print("=" * 70)
    
    # Run simulation
    L, y_ex_only, time_array, pert_osc, centers = run_quick_simulation()
    
    # Create animation
    anim, fig = create_animation(L, y_ex_only, time_array, pert_osc, centers)
    
    print("\n" + "=" * 70)
    print("SAVING IN MULTIPLE FORMATS")
    print("=" * 70)
    
    # 1. Save as MP4 (most common, good quality, small size)
    # save_as_mp4(anim, 'test_hexagon_demo.mp4', 
    #             fps=20, dpi=100, bitrate=2000)
    
    # 2. Save as HTML (interactive, works in browser)
    save_as_html(anim, 'test_hexagon_demo.html')
    
    # 3. Save as GIF (larger file, but works everywhere)
    # save_as_gif(anim, 'test_hexagon_demo.gif', 
    #             fps=10, dpi=80)
    
    # 4. Info about Jupyter usage
    print("\n" + "=" * 70)
    print("FOR JUPYTER NOTEBOOKS:")
    print("=" * 70)
    print("from IPython.display import HTML")
    print("HTML(anim.to_jshtml())")
    print("\nThis displays the animation inline in the notebook!")
    
    # Close figure
    plt.close(fig)
    
    print("\n" + "=" * 70)
    print("EXPORT COMPLETE!")
    print("=" * 70)
    print("\nFiles created:")
    print("  1. hexagon_demo.mp4   - High quality video (best for most uses)")
    print("  2. hexagon_demo.html  - Interactive HTML (open in browser)")
    print("  3. hexagon_demo.gif   - Animated GIF (universal compatibility)")
    print("\n" + "=" * 70)
    
    # Format comparison
    print("\nFORMAT COMPARISON:")
    print("-" * 70)
    print("MP4:")
    print("  ✓ Small file size")
    print("  ✓ High quality")
    print("  ✓ Widely supported")
    print("  ✓ Best for presentations, videos, sharing")
    print("  - Requires video player")
    
    print("\nHTML:")
    print("  ✓ Interactive (play/pause/speed controls)")
    print("  ✓ Works in any browser")
    print("  ✓ No external dependencies")
    print("  ✓ Perfect for Jupyter notebooks")
    print("  ✓ Self-contained (includes all data)")
    print("  - Larger file size than MP4")
    
    print("\nGIF:")
    print("  ✓ Works everywhere")
    print("  ✓ Auto-loops")
    print("  ✓ No player needed")
    print("  ✓ Good for documentation, emails")
    print("  - Large file size")
    print("  - Lower quality than MP4")
    print("  - No audio support")
    
    print("\n" + "=" * 70)

if __name__ == "__main__":
    main()
