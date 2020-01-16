import tau_eff
import tau_bootstrap
import tau_errors
import line_fitting
import compare
import numpy as np

lowz_tau, lowz_red = tau_eff.HI_opacity("lowz") #the actual stack

lowz_tau_matrix = tau_bootstrap.HI_opacity(5000,"lowz") #the bootstrap

lowz_sigma = tau_errors.covariance(lowz_tau, lowz_tau_matrix, 5000, lowz_red, "lowz" )
    

hiz_tau, hiz_red = tau_eff.HI_opacity("hiz")

hiz_tau_matrix = tau_bootstrap.HI_opacity(5000,"hiz")

hiz_sigma = tau_errors.covariance(hiz_tau, hiz_tau_matrix, 5000, hiz_red, "hiz")

#finally we use both bins!
tau = np.append(lowz_tau, hiz_tau)
red = np.append(lowz_red, hiz_red)
error = np.append(lowz_sigma, hiz_sigma)

mask = [~np.isnan(tau)] #to ensure all data points are valid
tau = tau[mask]
red = red[mask]
error = error[mask]

print(min(red), max(red))

print(len(tau))


#my own line
scale, power, z_space = line_fitting.chi_grid(tau, red, error)

#comparison plots
compare.literature(tau,red,error,z_space)







    
    
