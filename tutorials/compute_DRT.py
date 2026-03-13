import numpy as np
from math import pi, log, exp


def compute_A_complex(freq_vec, log_tau_vec):
    
    y = np.copy(log_tau_vec)
    N_f = freq_vec.size
    N_y = y.size

    out_A = np.zeros((N_f, N_y), dtype='cdouble')

    for m in range(0, N_f):
        for n in range(0, N_y):
            lhs = 0
            rhs = 0
            if n>0:
                y_lhs_mid = (y[n-1]+y[n])/2
                y_lhs_right = y[n]
                Delta_y = y[n] - y[n-1]
                I_lhs_mid = 2.0/(1.0+1j*2.0*pi*freq_vec[m]*exp(y_lhs_mid))
                I_lhs_right = 1.0/(1.0+1j*2.0*pi*freq_vec[m]*exp(y_lhs_right))
                lhs = Delta_y/6.0*(I_lhs_mid+I_lhs_right)

            if n<N_y-1:
                y_rhs_left = y[n]
                y_rhs_mid = (y[n]+y[n+1])/2
                Delta_y = y[n+1] - y[n]
                I_rhs_left = 1.0/(1.0+1j*2.0*pi*freq_vec[m]*exp(y_rhs_left))
                I_rhs_mid = 2.0/(1.0+1j*2.0*pi*freq_vec[m]*exp(y_rhs_mid))
                rhs = Delta_y/6.0*(I_rhs_left+I_rhs_mid)

            out_A[m, n] = lhs + rhs
            
    return out_A

def compute_L1(log_tau):
    
    N_tau = log_tau.size
    out_L1 = np.zeros((N_tau-1, N_tau))
    
    for n in range(0, N_tau-1):

        delta_loc = log_tau[n+1]-log_tau[n]
        
        out_L1[n,n] = 1./delta_loc
        out_L1[n,n+1] = -1./delta_loc

    return out_L1


# def compute_L2(log_tau):
#     '''
#     Create a differentiation matrix for second-order numerical derivatives.

#     :param log_tau: Vector of nodal points.
#     :return: Second-order numerical differentiation matrix.
#     '''

#     N_tau = len(log_tau)

#     # Initialize the differentiation matrix with zeros
#     out_L = np.zeros((N_tau, N_tau))

#     # Fill in the coefficients for each row in the matrix
#     for iter in range(1, N_tau-1):
#         h1 = log_tau[iter] - log_tau[iter-1]
#         h2 = log_tau[iter+1] - log_tau[iter]
#         out_L[iter, iter - 1] = (2 * h2**2) / (h1 * h2 * (h1**2 + h1 * h2 + h2**2))
#         out_L[iter, iter] = (-2 * (h2**2 - h1**2)) / (h1 * h2 * (h1**2 + h1 * h2 + h2**2))
#         out_L[iter, iter+1] = (2 * h1**2) / (h1 * h2 * (h1**2 + h1 * h2 + h2**2))

#     # Handle boundaries (forward/backward difference or other methods)
#     # These can be adjusted based on specific requirements
#     out_L[0, 0] = -2 / (log_tau[1] - log_tau[0])**2
#     out_L[0, 1] = 2 / (log_tau[1] - log_tau[0])**2
#     out_L[-1, -2] = 2 / (log_tau[-1] - log_tau[-2])**2
#     out_L[-1, -1] = -2 / (log_tau[-1] - log_tau[-2])**2

#     return out_L


def compute_L2_uniform(log_tau):
    
    N_tau = log_tau.size    
    out_L2 = np.zeros((N_tau-2, N_tau))
    
    for n in range(0, N_tau-2):

        delta_loc = log_tau[n+1]-log_tau[n]
        
        out_L2[n,n] = 1./(delta_loc**2)
        out_L2[n,n+1] = -2./(delta_loc**2)
        out_L2[n,n+2] = 1./(delta_loc**2)

    return out_L2


def compute_L2(log_tau):
    # https://pure.rug.nl/ws/files/3332271/1992JEngMathVeldman.pdf
    # Get the number of elements in log_tau
    N_tau = len(log_tau)
    
    # Initialize the output matrix with zeros, dimensions (N_tau - 2) x N_tau
    out_L2 = np.zeros((N_tau - 2, N_tau))

    # Loop through log_tau from the second to the second-to-last element
    for n in range(1, N_tau - 1):
        # Calculate the differences between consecutive elements
        h1 = log_tau[n] - log_tau[n-1]
        h2 = log_tau[n+1] - log_tau[n]
        
        # Calculate the sum of the differences
        h_sum = h1 + h2

        # Fill in the elements of the second derivative matrix
        out_L2[n-1, n-1] = 2 / (h1 * h_sum)
        out_L2[n-1, n]   = -2 * (1/h1 + 1/h2) / h_sum
        out_L2[n-1, n+1] = 2 / (h2 * h_sum)

    # Return the computed second derivative matrix
    return out_L2


def compute_L2_aristo(log_tau, sigma):
    '''
    Computes the augmented second derivative matrix with aristotelic boundary conditions.

    Parameters:
    - log_tau (numpy.ndarray): Input array of logarithmic tau values.
    - sigma (float): Standard deviation or scaling factor for the boundary conditions.

    Returns:
    - numpy.ndarray: Augmented second derivative matrix.
    '''
    # Compute the raw second derivative matrix using the compute_L2 function
    L2_standard = compute_L2(log_tau)
    
    # Get the number of elements in log_tau
    N_tau = len(log_tau)

    # Create a top vector with (1/sigma, 0, ..., 0)
    top_vector = np.zeros((1, N_tau))
    top_vector[0, 0] = 1 / sigma

    # Create a bottom vector with (0, 0, ..., 1/sigma)
    bottom_vector = np.zeros((1, N_tau))
    bottom_vector[0, -1] = 1 / sigma

    # Vertically stack the top vector, raw second derivative matrix, and bottom vector
    out_L2 = np.vstack([top_vector, L2_standard, bottom_vector])

    # Return the augmented second derivative matrix
    return out_L2


def compute_L2_neumann(log_tau, sigma):

    N_tau = log_tau.size
    out_L2 = np.zeros((N_tau, N_tau))

    # Interior points
    for n in range(1, N_tau-1):
        h1 = log_tau[n]-log_tau[n-1]
        h2 = log_tau[n+1]-log_tau[n]
        h_sum = h1+h2

        out_L2[n, n-1] = 2/(h1*h_sum)
        out_L2[n, n]   = -2*(1/h1+1/h2)/h_sum
        out_L2[n, n+1] = 2/(h2*h_sum)

    # Neumann boundary conditions
    # Left boundary (using forward difference)
    h1 = log_tau[1]-log_tau[0]
    out_L2[0, 0] = -1/(sigma*h1)
    out_L2[0, 1] = 1/(sigma*h1)

    # Right boundary (using backward difference)
    hN = log_tau[-1]-log_tau[-2]
    out_L2[-1, -2] = -1/(sigma*hN)
    out_L2[-1, -1] = 1/(sigma*hN)

    return out_L2

def compute_L3(log_tau):
    
    N_tau = log_tau.size    
    out_L3 = np.zeros((N_tau-3, N_tau))
    
    for n in range(0, N_tau-3):

        delta_loc = log_tau[n+1]-log_tau[n]
        
        out_L3[n,n] = -1./(delta_loc**3)
        out_L3[n,n+1] = 3./(delta_loc**3)
        out_L3[n,n+2] = -3./(delta_loc**3)
        out_L3[n,n+3] = 1./(delta_loc**3)

    return out_L3

def entropy(A, x):
    epsilon = 1E-10
    Ax = A@x
    p_temp = np.abs(Ax)
    p_sum = np.sum(p_temp)
    p = p_temp/p_sum + epsilon

    # Entropy calculation
    out_entropy_loc = p*np.log(p)
    out_entropy = -np.sum(out_entropy_loc)

    return out_entropy

def entropy_gradient(A, x):
    epsilon = 1E-10
    Ax = A@x
    p_temp = np.abs(Ax)
    p_sum = np.sum(p_temp)
    p = p_temp/p_sum + epsilon

    # Gradient calculation
    dp_dx = (np.sign(Ax)@A)/p_sum-(p_temp/p_sum**2)*np.sum(np.sign(Ax)@A, axis=1, keepdims=True)
    d_entropy_dp = -np.log(p)-1
    out_grad_entropy = dp_dx.T@d_entropy_dp

    return out_grad_entropy