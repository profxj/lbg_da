import astropy.units as u
from astropy.coordinates import SkyCoord
from dustmaps.bayestar import BayestarQuery
from astropy.cosmology import WMAP9 as cosmo
import stack
import clamato_read
import bootstrap

#the full composite creation process

z_min = float(input("what is the minimum redshift (z) bin value?"))

z_max = float(input("what is the maximum redshift (z) bin value?"))

#the read in

spec, red, radec = clamato_read.xspec([z_min, z_max]) 

# unreddening

Mpc = cosmo.comoving_distance(red)

bayestar = BayestarQuery()

coords = [SkyCoord(radec[i][0] * u.deg, radec[i][1] * u.deg, distance=Mpc[i], frame='fk5') for i in range(len(spec))]

ebv = [bayestar(i) for i in coords]

# the stack

s2n = float(input("What is the minimum S/N value you want to select from the sample?"))

dv = float(input("What velocity dispersion (km/s) do you want to rebin to?"))

stack, rest_spec = stack.lbg_smash(spec, red, ebv, s2n, dv)

N = int(input("How many times do you want to iterate over the bootstrap?"))

#bootstrap

composite = bootstrap.cov_err(rest_spec, N, stack)

location = input("Where do you want to save the composite?")

composite.write_to_fits(location)

