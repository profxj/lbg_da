def lya_opacity(stack, fit, filename):

    """
    :param stack: the stack formatted in the same way as for SB99 analysis
    :param fit: the fitted model to that stack (also from SB99)
    :param filename: what the file will be saved as
    :return: fits file with following data:

             z_emm: redshift values per pixel
             t_eff: the corresponding effective opacity
             t_sig: errors on t_eff

    """

    import numpy as np
    from astropy.io import fits
    import matplotlib.pyplot as plt

    wave = stack[1].data["wavelength"]

    flux = stack[1].data["flux"]

    sigma = stack[1].data["noise"]

    z_med = stack[1].header["redshift"]

    model = fit[0].data

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

    t_mask = [flux[i] < model[i] for i in range(npix)]

    z_emm = ((wave[t_mask] / 1216) * (1 + z_med) - 1)  # values scaled by lya and median z value of the composite

    t_eff = -np.log(flux[t_mask] / model[t_mask])  # an effective opacity

    # errors! #teff = -ln(Fobs/Fmod)

    t_sig = np.sqrt( (-sigma[t_mask]/flux[t_mask] )**2 )

    # a simple plot

    plt.errorbar(z_emm, t_eff, yerr=t_sig, fmt='o')
    plt.xlabel("$z$", fontsize=18)
    plt.ylabel("$τ_{eff}$", fontsize=18)
    plt.show()

    # writing out the data

    c_z = fits.Column(name='redshift', array=z_emm, format="E")

    c_tau = fits.Column(name='tau', array=t_eff, format="E")

    c_sig = fits.Column(name='noise', array=t_sig, format="E")

    t = fits.BinTableHDU.from_columns([c_z, c_tau, c_sig])

    path = "/home/jsm/PycharmProjects/tau_eff/fits/Tau_eff/"

    t.writeto(path+filename)

    return z_emm, t_eff, np.array(t_sig)