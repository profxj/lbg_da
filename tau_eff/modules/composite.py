import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from dustmaps.bayestar import BayestarQuery
from astropy.cosmology import WMAP9 as cosmo
import numpy as np
import stack
import clamato_read
import bootstrap

z_min = float(input("what is the minimum redshift (z) bin value?"))

z_max = float(input("what is the maximum redshift (z) bin value?"))

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

# the bootstrap

N = int(input("How many times do you want to iterate over the bootstrap?"))

composite = bootstrap.cov_err(rest_spec, N, stack)

# formatting the composite for SB99

mask = [1000*u.AA < i < 1500*u.AA for i in composite.wavelength]  # only the data I want 1000-1500A

wave = composite.wavelength[mask]

flux = composite.flux[mask]

sig = composite.sig[mask]

c_wave = fits.Column(name='WAVELENGTH', array=wave, unit="angstroms", format="E")

c_flux = fits.Column(name='FLUX', array=flux, unit="relative flux", format="E")

c_noise = fits.Column(name='NOISE', array=sig, unit="relative flux", format="E")

t = fits.BinTableHDU.from_columns([c_wave, c_flux, c_noise])

hdr = t.header

hdr.set("REDSHIFT", z)

location2 = input("where do you want to save the formatted composite?")

t.writeto(location2)





