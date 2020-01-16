import numpy as np
#simple Chi^2
def chisquare(data,model,sigma):
    
    num = (data-model)**2
    
    den = sigma**2
    
    return sum(num/den)

#Faucher-Giguere  (2 free parameters)
def F_G(z):
    
    τ_eff = 0.0018*(1+z)**3.92
        
    return τ_eff

#Becker  (3 free parameters)
def Becker(z):
    
    τ_eff = .751*((1+z)/(1+3.5))**2.9 - .132
        
    return τ_eff

#Thomas  (2 free parameters)
def Thomas(z):
    
    return .0089*(1+z)**2.55

#Schaye et al

z_vals_s = np.array([1.756,2.000,1.799,2.034,1.843,2.080,1.986,2.224,2.003,2.239,1.998,
            2.242,2.010,2.256,2.103,2.366,2.217,2.496,2.243,2.506,2.308,2.572,
            2.447,2.715,2.509,2.805,2.626,2.920,2.670,3.009,2.752,3.078,3.058,
            3.382,3.088,3.411,3.708,3.912,3.517,3.862,3.862,4.287])


τ_eff_s = np.array([0.099,0.093,0.120,0.162,0.224,0.124,0.128,0.164,0.180,0.113,0.139,
           0.149,0.115,0.156,0.131,0.175,0.137,0.214,0.180,0.205,0.234,0.283,
           0.177,0.308,0.182,0.273,0.279,0.343,0.232,0.360,0.329,0.271,0.423,
           0.496,0.366,0.445,0.705,0.811,0.644,0.843,0.839,0.827])


τ_err_s = np.array([0.016,0.016,0.026,0.026,0.033,0.019,0.019,0.026,0.027,0.018,0.021,
           0.018,0.014,0.021,0.019,0.022,0.020,0.029,0.025,0.022,0.029,0.032,
           0.025,0.041,0.021,0.030,0.033,0.037,0.023,0.040,0.043,0.030,0.034,
           0.048,0.039,0.041,0.074,0.065,0.056,0.075,0.069,0.060])
    

#Kirkman 
def D_a_to_τ_eff(D_a):
    
    τ_eff = -np.log(1-D_a)
    
    return τ_eff

z_vals_k = [1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2]

D_a_k = np.array([0.098,0.099,0.128,0.143,0.202,0.213,0.250,0.252,0.338])

D_a_k_sig = np.array([0.015,0.006,0.006,0.013,0.014,0.014,0.016,0.015,0.024])
                 
τ_eff_k = D_a_to_τ_eff(D_a_k)

τ_err_k =  D_a_to_τ_eff(D_a_k_sig)

#now for the comparison

def literature(tau,red,error,z_space):
    
    import matplotlib.pyplot as plt

    FG_array = F_G(z_space)
    FG_reduce = len(FG_array)+2 
    FG_chi = chisquare(tau, FG_array, error)/FG_reduce


    Becker_array = Becker(z_space)
    Becker_reduce = len(Becker_array)+3 
    Becker_chi = chisquare(tau, Becker_array, error)/Becker_reduce


    Thomas_array = Thomas(z_space)
    Thomas_reduce = len(Thomas_array)+2
    Thomas_chi = chisquare(tau, Thomas_array, error)/Thomas_reduce
    
    print("Thomas Chi^2:"+str(Thomas_chi))
    print("Faucher-Giguere Chi^2:"+str(FG_chi))
    print("Becker Chi^2:"+str(Becker_chi))
    
    
    plt.figure(figsize=(10,5))

    plt.errorbar(red, tau, yerr=error,fmt="o", label="This Work",color="gray")

    plt.errorbar(z_vals_k,τ_eff_k, yerr = τ_err_k,
                 capsize=5, marker="*",label="Kirkman et al. 2005",
                 ls="none",color="orange", markersize=10)
                 
    plt.errorbar(z_vals_s,τ_eff_s, yerr = τ_err_s,
                 capsize=5, marker="d",label="Schaye et al. 2003",
                 ls="none",color="darkblue", markersize=10)    
             
    plt.plot(z_space, Thomas_array, color = "#f03b20",label = "Thomas et al. 2017")
    
    plt.ylim(-.1, .7)
    plt.xlim(1.95, 2.5)
    plt.xlabel("$z$",fontsize=15)
    plt.ylabel("$τ_{eff}$",fontsize=15)
    plt.legend()
    plt.savefig("/Users/jsmonzon/IGM-UCSC/figures/comptau.pdf")
    plt.show()
    




