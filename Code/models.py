import numpy as np

def sigmoid(u):

    return np.tanh(u)

def two_oscillators(t, y, h_ex, h_in, pars, frac_E, frac_I, sr, time_stop, pert):
               
    tau_ex, tau_in, c2, c4, c_EE, c_EI = pars

    time_index = int(t*sr)

    if time_index >= time_stop*sr:

        dydt = np.zeros(4)

        return dydt


    else:
        
        dydt = (

            (pert[time_index] - y[0] - c2*sigmoid(y[1]) + c_EE*sigmoid(y[0]+frac_E*y[2]))*tau_ex,
             
            (h_in             - y[1] - c4*sigmoid(y[1]) + c_EI*sigmoid(y[0]+frac_I*y[2]))*tau_in,

             
            (h_ex             - y[2] - c2*sigmoid(y[3]) + c_EE*sigmoid(y[2]+frac_E*y[0]))*tau_ex,
            
            (h_in             - y[3] - c4*sigmoid(y[3]) + c_EI*sigmoid(y[2]+frac_I*y[0]))*tau_in
           )

    return dydt

def N_oscillators(y, t, N, h_ex_rand, h_in_rand, 
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

        if osc in pert_osc_list:  # Changed from == to 'in'
                
            dy_ex[osc] = (pert[time_index] - y_ex[osc] - c2*sigmoid(y_in[osc]) + coup_EE)*tau_ex 
            dy_in[osc] = (h_in_rand[osc]   - y_in[osc] - c4*sigmoid(y_in[osc]) + coup_EI)*tau_in
            
        else:
                
            dy_ex[osc] = (h_ex_rand[osc]   - y_ex[osc] - c2*sigmoid(y_in[osc]) + coup_EE)*tau_ex 
            dy_in[osc] = (h_in_rand[osc]   - y_in[osc] - c4*sigmoid(y_in[osc]) + coup_EI)*tau_in
            
    # Combine Variables
    dydt[:-1:2] = dy_ex
    dydt[1: :2] = dy_in
    
    return dydt


# def N_oscillators(y, t, N, h_ex_rand, h_in_rand, 
#                         coupling_matrix_EE, coupling_matrix_EI, cconst_E, cconst_I, pars, sr, time_stop, pert, pert_osc_list):
        
#     tau_ex, tau_in, c2, c4 = pars
    
#     time_index = int(t*sr)
    
#     if time_index >= time_stop*sr:
#         dydt = zeros(2*N)
        
#         return dydt
#     # Separate Variables
#     y_ex = y[:-1:2]
#     y_in = y[1::2]
#     dy_ex, dy_in = zeros(N), zeros(N)
#     dydt = zeros(2*N)
#     for osc in arange(N):
#         coup_EE = cconst_E*sigmoid(sum(coupling_matrix_EE[:, osc] * y_ex))
#         coup_EI = cconst_I*sigmoid(sum(coupling_matrix_EI[:, osc] * y_ex))

#         if osc in pert_osc_list:  # Changed from == to 'in'
                
#             dy_ex[osc] = (pert[time_index] - y_ex[osc] - c2*sigmoid(y_in[osc]) + coup_EE)*tau_ex 
#             dy_in[osc] = (h_in_rand[osc]   - y_in[osc] - c4*sigmoid(y_in[osc]) + coup_EI)*tau_in
            
#         else:
                
#             dy_ex[osc] = (h_ex_rand[osc]   - y_ex[osc] - c2*sigmoid(y_in[osc]) + coup_EE)*tau_ex 
#             dy_in[osc] = (h_in_rand[osc]   - y_in[osc] - c4*sigmoid(y_in[osc]) + coup_EI)*tau_in
            
#     # Combine Variables
#     dydt[:-1:2] = dy_ex
#     dydt[1: :2] = dy_in
    
#     return dydt


def m01_2var_feedback_inhibition(t, variables, a1, b1, a2, b2, k_max, K_m, k_i, n, m, q):
    """Single reaction with feedback inhibition"""
    S, P = variables
    
    enzymatic_rate = (k_max * S**n) / (K_m**m + S**m) / (1 + k_i * P**q)
    
    dSdt = a1 - b1 * S - enzymatic_rate
    dPdt = a2 - b2 * P + enzymatic_rate
    
    return [dSdt, dPdt]

def m02_2var_forward_inhibition(t, variables, a1, b1, a2, b2, k_max, K_m, n, m):
    """Single reaction with forward inhibition"""
    S, P = variables
    
    enzymatic_rate = (k_max * S**n) / (K_m**n + S**m)
    
    dSdt = a1 - b1 * S - enzymatic_rate
    dPdt = a2 - b2 * P + enzymatic_rate
    
    return [dSdt, dPdt]
