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
    
    # truncating each spectra to only include wavelength values from 950A to 1500A
    
    r = range(len(spec))

    temp = [np.asarray(spec[i].wavelength / (1 + red[i])) for i in r]
    
    wv_coverage = [(1050 < entry) & (entry < 1400) for entry in temp]
    
    wave = np.asarray([spec[i].wavelength[wv_coverage[i]] for i in r])
    
    flux = np.asarray([spec[i].flux[wv_coverage[i]] for i in r])
    
    error = np.asarray([spec[i].sig[wv_coverage[i]] for i in r])
    
    #the normalization what if I dont normailize?
    
    temp = [np.asarray(wave[i] / (1 + red[i])) for i in r]

    wv_norm = [(1260 < entry) & (entry < 1304) for entry in temp]  # just between the SiII lines

    flux_range = np.asarray([flux[i][wv_norm[i]] for i in r])

    medians = np.asarray([np.median(flux_range[i]) for i in r])

    norm_flux = np.asarray([(flux[i] / medians[i]) for i in r])

    t_spec = []

    for i in r: #maybe try to not messing with indecies?
        
        t_spec.append(XSpectrum1D(wave[i][0:935], norm_flux[i][0:935], error[i][0:935])) #[0:935]

    collate = ltsu.collate(t_spec)  # to get a single XSpec object

    fix = []  # to make sure none of the values included will mess up the stack

    a = collate.data['wave']

    for i in range(collate.nspec):

        if np.min(a[i]) == 0.0:
            
            fix.append(i)

    if len(fix) > 0:

        print("deleting these spectra", fix)

        good_spec = np.delete(t_spec,fix) # exclude all spec[fix]

        red = np.delete(red,fix)

        collate = ltsu.collate(good_spec)
        
        rest_spec = ltsu.rebin_to_rest(collate, red, 300 * u.km / u.s, grow_bad_sig=True)
        
    else:
        
        rest_spec = ltsu.rebin_to_rest(collate, red, 300 * u.km / u.s, grow_bad_sig=True)
        
    return rest_spec