def HI_opacity(zbin):

    """
    :param stack: the stack formatted in the same way as for SB99 analysis
    :param fit: the fitted model to that stack (also from SB99)
    :param filename: what the file will be saved as
    :return: fits file with following data:

             z_emm: redshift values per pixel
             t_eff: the corresponding effective opacity

    """

    import numpy as np
    from astropy.io import fits
    
    
    stack = fits.open("/Users/jsmonzon/IGM-UCSC/fits/composites/"+zbin+"/stack.fits")

    model = fits.open("/Users/jsmonzon/IGM-UCSC/fits/composites/"+zbin+"/model.fits")
    
    wave = stack[1].data["wavelength"]
    
    flux = stack[1].data["flux"]
    
    sigma = stack[1].data["noise"]
    
    model = model[0].data
    
    z_med = stack[1].header["redshift"]
    
    # initial mask as we only want the opacity in the lya forest

    lya_mask = [1070 < i < 1170 for i in wave]

    wave = wave[lya_mask]

    flux = flux[lya_mask]

    model = model[lya_mask]

    sigma = sigma[lya_mask]

    # second mask to ensure we don't include any ISM features

    ism = np.array([1083.99, 1117.97, 1122.52, 1128.01, 1144.93, 1152.81])

    buffer = 5  # angstroms

    masks = []

    for j in ism:
        masks.append([(j - buffer < i) & (i < j + buffer) for i in wave])

    ism_mask = np.invert([any(t) for t in zip(masks[0], masks[1], masks[2], masks[3], masks[4], masks[5])])

    wave = wave[ism_mask]

    flux = flux[ism_mask]

    model = model[ism_mask]

    sigma = sigma[ism_mask]

    npix = len(wave)

    # final mask to ensure that we get positive values of tau
        
    tau = np.zeros(npix)

    for j in range(npix):
            
        if flux[j] < model[j]:
                                
            t_eff = -np.log(flux[j] / model[j])  # an effective opacity
                
            tau[j] = t_eff
                
        else:
                
            tau[j] = np.nan
                        
    z_abs = (wave/ 1216) * (1 + z_med) - 1 # values scaled by lya and median z value 

    return tau, z_abs