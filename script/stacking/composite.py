import clamato_read
import cuts
import milkyway
import prep
import bootstrap
import numpy as np

#grab the spectra
spec, red, coord = clamato_read.xspec([2.25,2.75])

print("---------CLAMATO data read in--------")

#cut the sample
spec, red, coord = cuts.lyman_limit(spec, red, coord)
spec, red, coord = cuts.s2n(spec, red, coord)

print("--------data sample cut--------")

#deredden
spec = milkyway.extinction(spec, red, coord)

print("--------spectra de-reddened--------")

#creating the two bins
lowz_cut = [red < 2.5]
hiz_cut = [red > 2.5]

low_spec, low_red = spec[lowz_cut], red[lowz_cut]
hi_spec, hi_red = spec[hiz_cut], red[hiz_cut]

print("lowz")
print("z median =", np.median(low_red), "z std = ",np.std(low_red), "nspec =",len(low_spec))

print("hiz")
print("z median =", np.median(hi_red), "z std = ",np.std(hi_red), "nspec =",len(hi_spec))

"""
low_rest = prep.clean(low_spec, low_red)
hi_rest = prep.clean(hi_spec, hi_red)

print("--------spectra seperated into bins and cleaned for stacking--------")

#the stacks!
print("--------low_z bootstrap beginning--------")
bootstrap.stacks(low_rest,low_red,5000,"lowz")

print("--------hi_z bootstrap beginning--------")
bootstrap.stacks(hi_rest,hi_red,5000,"hiz")
"""