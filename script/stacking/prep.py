def clean(spec, red):
    
    """
    :param spec: [array] of XSpectrum1D objects
    :param red: [array] of their redshift values (z)
    :returns: [array] of truncatted, normalized and restframed XSpectrum1D objects
    """
    
    import numpy as np
    import astropy.units as u
    from linetools.spectra import utils as ltsu
    from linetools.spectra.xspectrum1d import XSpectrum1D

    r = range(len(spec))

    temp = [np.asarray(spec[i].wavelength / (1 + red[i])) for i in r]
    # truncating each spectra to only include wavelength values from 1000 to 1450
    wv_coverage = [(1000 < entry) & (entry < 1450) for entry in temp]
    
    wave = np.asarray([spec[i].wavelength[wv_coverage[i]] for i in r])
    flux = np.asarray([spec[i].flux[wv_coverage[i]] for i in r])
    error = np.asarray([spec[i].sig[wv_coverage[i]] for i in r])

    temp = [np.asarray(wave[i] / (1 + red[i])) for i in r]
    # normalizing just between the SiII lines
    wv_norm = [(1260 < entry) & (entry < 1304) for entry in temp]  # just between the SiII lines

    flux_range = np.asarray([flux[i][wv_norm[i]] for i in r])
    medians = np.asarray([np.median(flux_range[i]) for i in r])
    norm_flux = np.asarray([(flux[i] / medians[i]) for i in r])

    # getting a single XSpec object to feed into the stack
    spec = [XSpectrum1D(wave[i], norm_flux[i], error[i]) for i in r]
    collate = ltsu.collate(spec)
    rest_spec = ltsu.rebin_to_rest(collate, red, 300 * u.km / u.s, grow_bad_sig=True)

    return rest_spec