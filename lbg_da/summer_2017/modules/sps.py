import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from linetools.spectra import utils as ltsu
from linetools.spectra.xspectrum1d import XSpectrum1D
import astropy.units as u
from array import array
import fsps

#Stellar population object.
model = fsps.StellarPopulation(imf_type=0, zcontinuous=2, logzsol=.9, zred=2.5, dust_type=2, dust2 = .9)

wave, flux = model.get_spectrum(tage= 6, peraa = True)


#the stack to scale the model
stack = XSpectrum1D.from_file("/home/sebastian/lbg_da/fits/fullstack.fits")


#converting the model to Xspec and scaling it
xmodel = XSpectrum1D(wave= wave, flux= flux)

rest_xmodel = ltsu.rebin(xmodel, stack.wavelength)

factor_2 = np.median(stack.flux)/np.median(rest_xmodel.flux)

xfsps =  XSpectrum1D(wave= rest_xmodel.wavelength, flux= rest_xmodel.flux*factor_2)

xfsps.write_to_fits("sps.fits")