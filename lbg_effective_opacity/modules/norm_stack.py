def stack(spec , red , wv_norm , s2n , dv):

    """

    :param spec: array of XSpectrum1D objects: use clamato_read.py
    :param red: array of their redshift values (z)
    :param wv_norm: array of the wavelength range (angstroms) in the rest frame in which to normalize, ex:[1450,1500]
    :param s2n: float lowest value of S/N. will be computed over wv_norm
    :param dv: velocity width of new pixels for rebin (will be measured in km/s)
    :return: a single composite XSpectrum1D object only flux & wavelength values.


    """

    # Note: this routine normalizes spectra individually before stacking them,
    # see stack_norm.py for the opposite

    import numpy as np
    from astropy.io import fits
    from astropy.table import Table as Table
    import matplotlib.pyplot as plt
    from linetools.spectra import utils as ltsu
    from linetools.spectra.xspectrum1d import XSpectrum1D
    import astropy.units as u
    from astropy import constants as const

    print("Number of spectra provided =", len(spec))


# initial S/N cut

    r = range(len(spec))

    temp = [np.asarray(spec[i].wavelength / (1 + red[i])) for i in r]  # rest-frame

    norm = [(np.min(wv_norm) < entry) & (entry < np.max(wv_norm)) for entry in temp]

    signal = np.asarray([spec[i].flux[norm[i]] for i in r])

    noise = np.asarray([spec[i].sig[norm[i]] for i in r])

    s2n_array = np.asarray([np.median(signal[i] / noise[i]) for i in r])

    s2n_cut = [i > s2n for i in s2n_array]

    spec = spec[s2n_cut]

    red = red[s2n_cut]

    print("S/N cut left", len(spec), "spectra")


# the normalization

    r = range(len(spec))

    temp = [np.asarray(spec[i].wavelength / (1 + red[i])) for i in r]

    wv_norm = [(1260 < entry) & (entry < 1304) for entry in temp]  # just between the SiII lines.

    flux_range = np.asarray([spec[i].flux[wv_norm[i]] for i in r])

    medians = np.asarray([np.median(flux_range[i]) for i in r])

    norm_flux = np.asarray([(spec[i].flux / medians[i]) for i in r])

    norm_spec = [XSpectrum1D(spec[i].wavelength, norm_flux[i],spec[i].sig) for i in r] # the normalized spec array


# truncating each spectra to only include wavelength values from 950A to 1500A

    r = range(len(norm_spec))

    temp = np.array([norm_spec[i].wavelength / (1 + red[i]) for i in r])

    trim = []

    for i in r:
        trim.append([950 * u.AA <= temp[i][j] <= 1500 * u.AA for j in range(len(temp[i]))])

    t_spec = []

    for i in r:
        t_spec.append(XSpectrum1D(norm_spec[i].wavelength[trim[i]][0:935],
                                  norm_spec[i].flux[trim[i]][0:935],
                                  norm_spec[i].flux[trim[i]][0:935]))

    collate = ltsu.collate(t_spec) # to get a single XSpec object

    print("Number of spectra to be stacked =", collate.nspec)

# rebin to rest frame with dv dispersion

    rest_spec = ltsu.rebin_to_rest(collate, red, dv*u.km/u.s, grow_bad_sig=True)

# the stack! Note: the error array was not normalized; do not trust it

    stack = ltsu.smash_spectra(rest_spec)

    stack.plot()

    return stack, rest_spec




