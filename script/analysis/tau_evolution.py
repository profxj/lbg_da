import numpy as np
from astropy.table import Table
import tau_eff

print("---------calculating tau_eff--------")

low_tau, low_z, low_mask, low_zmed = tau_eff.composite("lowz")

print("the lowz bin found",len(low_tau), "valid indicies from which to measure tau_eff")

hi_tau, hi_z, hi_mask, hi_zmed = tau_eff.composite("hiz")

print("the hiz bin found",len(hi_tau), "valid indicies from which to measure tau_eff")

print("---------bootstrap error calculation--------")

N = 5000

low_tmat = tau_eff.bootstrap("lowz", N, low_mask, low_zmed)

hi_tmat= tau_eff.bootstrap("hiz", N, hi_mask, hi_zmed)


tau_boot = np.concatenate((low_tmat, hi_tmat), axis = 1)
tau = np.append(low_tau, hi_tau)
red = np.append(low_z, hi_z)

print(tau.shape, tau_boot.shape)

#covar = tau_eff.covariance(tau, tau_boot)

tau_eff = np.array([red,tau])
np.save("/Users/jsmonzon/lbg_da/tau_data/test/tau.npy", tau_eff)
#np.save("/Users/jsmonzon/lbg_da/tau_data/test/covar.npy", covar)
np.save("/Users/jsmonzon/lbg_da/tau_data/test/bootstrap.npy", tau_boot)





#low_tsig = tau_eff.covariance(low_tau, low_tmat)

#hi_tsig = tau_eff.covariance(hi_tau, hi_tmat)




#print("---------writing out the results--------")

#low_info = Table([np.round(low_z,4), np.round(low_tau,4), np.round(low_tsig,4)], names=("z", "$τ_eff$", "$τ_sig$"))
#low_info.write("/Users/jsmonzon/lbg_da/figures/low.txt", format="latex", overwrite="True")

#hi_info = Table([np.round(hi_z,4), np.round(hi_tau,4), np.round(hi_tsig,4)], names=("z", "$τ_eff$", "$τ_sig$"))
#hi_info.write("/Users/jsmonzon/lbg_da/figures/hi.txt", format="latex", overwrite="True")

#tau = np.append(low_tau, hi_tau)
#red = np.append(low_z, hi_z)
#tau_sig = np.append(low_tsig, hi_tsig)

#tau_eff = np.array([red,tau,tau_sig])
#np.save("/Users/jsmonzon/lbg_da/tau_data/tau.npy", tau_eff)