def covariance(tau, tau_matrix, N, z, zbin):
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    N_matrix = np.array([tau_matrix[i] - tau for i in range(N)])

    tranpose = np.transpose(N_matrix)  # matrix math

    covariance = np.dot(tranpose, N_matrix)

    sigma = np.sqrt(np.diagonal(covariance) / (N - 1))
    
    fig, ax = plt.subplots(figsize=(10, 10))
    plt.errorbar(z, tau, yerr=sigma, fmt="o", label=zbin)
    plt.xlabel("Redshift (z)", fontsize = 15)
    plt.ylabel("Effective Opacity", fontsize = 15)
    plt.legend()
    plt.show()
    
    return sigma
    
    