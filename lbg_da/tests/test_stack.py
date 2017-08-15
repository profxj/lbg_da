import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from linetools.spectra import utils as ltsu
from linetools.spectra.xspectrum1d import XSpectrum1D
import astropy.units as u
from array import array

def test_on_5():



#reading in the five spec

    gal_1 = XSpectrum1D.from_file('/home/sebastian/Spectra/spec/c16_24_zph2.3_15045_F.fits')

    gal_2 = XSpectrum1D.from_file('/home/sebastian/Spectra/spec/c16_18_zph2.4_17949_F.fits')

    gal_3 = XSpectrum1D.from_file('/home/sebastian/Spectra/spec/c16_20_zph2.4_17955_F.fits')

    gal_4 = XSpectrum1D.from_file('/home/sebastian/Spectra/spec/c16_22_zph2.4_15040_F.fits')

    gal_5 = XSpectrum1D.from_file('/home/sebastian/Spectra/spec/c16_11_zph2.3_12434_F.fits')


#the in initial stack without normalizing the flux array


    spec_list = [gal_1, gal_2, gal_3, gal_4, gal_5]

    z_values = ([2.128, 1.94, 2.51, 2.54, 2.42])

    collate = ltsu.collate(spec_list)

    restspec = ltsu.rebin_to_rest(collate, z_values, 200 * u.km / u.s)

    stack_1 = ltsu.smash_spectra(restspec)

    return stack_1

#a new xspec object scaled to the highest median


    flux_av = []

    for i in range(restspec.nspec):
        flux_av.append(np.median(restspec[i].flux))

    factor = max(flux_av)

    scale = []

    for number in flux_av:
        scale.append(factor / (number))

    scaled_spec = []

    for i in range(restspec.nspec):
        scaled_spec.append(XSpectrum1D(restspec[i].wavelength, restspec[i].flux * scale[i], restspec[i].sig))

    scaled_col = ltsu.collate(scaled_spec)

    stack_2 = ltsu.smash_spectra(scaled_col)

    return stack_2

test_on_5()

