import numpy as np
import tau_eff

print("---------calculating tau_eff--------")

low_tau, low_z, low_mask, low_zmed = tau_eff.composite("lowz")

print("the lowz bin found",len(low_tau), "valid indicies from which to measure tau_eff")

hi_tau, hi_z, hi_mask, hi_zmed = tau_eff.composite("hiz")

print("the hiz bin found",len(hi_tau), "valid indicies from which to measure tau_eff")

print("---------bootstrap error calculation--------")

N = 5000

low_tmat, low_zmat = tau_eff.bootstrap("lowz", N, low_mask, low_zmed)

hi_tmat, hi_zmat = tau_eff.bootstrap("hiz", N, hi_mask, hi_zmed)

low_tsig = tau_eff.covariance(low_tau, low_tmat)

hi_tsig = tau_eff.covariance(hi_tau, hi_tmat)

print("---------writing out the results--------")

tau = np.append(low_tau, hi_tau)
red = np.append(low_z, hi_z)
tau_sig = np.append(low_tsig, hi_tsig)

tau_eff = np.array([red,tau,tau_sig])
np.save("/Users/jsmonzon/lbg_da/tau_data/tau_data.npy", tau_eff)