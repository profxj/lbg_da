#the function to fit to
def tau_ev(z,a,b):
    return a*((1+z)/(1+2.22))**b

# simple chi^2
def chisquare(data,model,sigma):
    num = (data-model)**2
    den = sigma**2
    return sum(num/den)

def chi_grid(tau, red, error):
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    ind = 100
    N = np.linspace(0.001,1,ind)#scale
    M = np.linspace(1.5,3,ind)#power
    chisq = np.zeros((ind,ind))#empty array

    for i, N_i in enumerate(N):
        for j, M_j in enumerate(M):
        
            model = tau_ev(red,N_i,M_j)
    
            chisq[i,j] = np.sum(chisquare(tau,model,error))
            
    data_len = len(tau)
    reduce = len(tau) + 2 #data points plus free parameters
    
    chisq = chisq/reduce #reduced chi^2
    chi_min = np.min(chisq)
    
    print("Mininum value in the grid: " + str(chi_min))

    plt.figure(figsize=(10,10))
    plt.contourf(N,M, chisq, cmap=plt.cm.magma_r)
    plt.xlabel("Scale Parameter (A)",fontsize=15)
    plt.ylabel("Power Parameter (B)",fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    cb = plt.colorbar()
    cb.set_label(label="$\chi^2$",fontsize=15)
    plt.show()

    Nb, Mb = np.where(chisq==chi_min)
    best_scale = float(N[Nb])
    best_power = float(M[Mb])

    print("Best scale parameter:" + str(best_scale))
    print("Best power parameter:" + str(best_power))
    
    #make the function
    
    z_space = np.linspace(1,4,data_len)
    evolution = tau_ev(z_space, best_scale, best_power)
    
    #plot it!
    
    plt.figure(figsize=(10,5))
    plt.errorbar(red, tau, yerr=error,fmt="o",
                 color = "gray", label = "Effective opacity")
    plt.plot(z_space, evolution, label="Best-Fit Power Law",
             color = "#f03b20", ls = "dashdot")
    plt.xlim(1.95, 2.5)
    plt.ylim(-.05,.65)
    plt.legend()
    plt.xlabel("$z$",fontsize=15)
    plt.ylabel("$τ_{eff}$",fontsize=15)
    #plt.savefig("/Users/jsmonzon/IGM-UCSC/figures/tauscatter.pdf")
    plt.show()
    
    return best_scale, best_power, z_space
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    