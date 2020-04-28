def composite(zbin):

    """
    :param zbin: (str) "lowz" or "hiz"

    :return:
        tau_eff: (numpy array) calculated tau values
        z_eff: (numpy array) corresponding redshift values
        lya_mask: (numpy array) mask to select valid lya values, used in the bootstrapping step
        z_med: (float) median redshift value for zbin
    """
    import numpy as np
    from astropy.io import fits

    stack = fits.open("/Users/jsmonzon/lbg_da/fits_data/composites/observed/"+zbin+"/composite.fits")
    model = fits.open("/Users/jsmonzon/lbg_da/fits_data/composites/modeled/"+zbin+"/sed.fits")

    wave = stack[1].data["wavelength"]
    flux = stack[1].data["flux"]
    z_med = stack[0].header["redshift"]
    model = model[0].data

    ism = np.array([1083.99, 1117.97, 1122.52, 1128.01, 1144.93, 1152.81])
    buffer = 3 # +- angstroms

    forest_mask = [] # between 1070-1170
    line_mask = [] # skipping the ism lines
    wave_ind = np.arange(0, len(wave))

    for i, lam in enumerate(wave):

        if 1070 < lam < 1170:

            if flux[i] < model[i]:

                forest_mask.append(i)

                for line in ism:

                    if (line - buffer) < lam < (line + buffer):
                        line_mask.append(i)

    lya_mask = np.isin(wave_ind, np.array(forest_mask)[~np.isin(forest_mask, line_mask)])

    tau_eff = -np.log(flux[lya_mask] / model[lya_mask])
    z_eff = (wave[lya_mask]/ 1216) * (1 + z_med) - 1

    tau_eff = kirkman_correction(tau_eff, z_eff, plot=False, p=False)

    return tau_eff, z_eff, lya_mask, z_med

def bootstrap(zbin, N, lya_mask, z_med):

    """
    :param zbin: (str) "lowz" or "hiz"
    :param N: (int) number of bootstrap iterations
    :param lya_mask: (numpy array) mask to select correct wavelengths
    :param z_med: (float) median redshift value for zbin

    :return:
        tau_matrix: (numpy array) tau values across bootstrap
        z_matrix: (numpy array) corresponding redshift values
    """
    import numpy as np
    from astropy.io import fits

    stacks = fits.open("/Users/jsmonzon/lbg_da/fits_data/bootstrap/observed/"+zbin+"/stacks.fits")
    wave = stacks[1].data["wavelength"]

    tau_matrix = []
    z_matrix = []

    for i in range(N):

        model = fits.open("/Users/jsmonzon/lbg_da/fits_data/bootstrap/modeled/"+zbin+"/sed_"+str(i)+".fits")

        model = model[0].data
        flux = stacks[2].data[i]

        tau_eff = -np.log(flux[lya_mask] / model[lya_mask])
        z_eff = (wave[lya_mask] / 1216) * (1 + z_med) - 1

        tau_eff = kirkman_correction(tau_eff, z_eff)

        tau_matrix.append(tau_eff)
        #z_matrix.append(z_eff)

    return np.array(tau_matrix) #, np.array(z_matrix)


def covariance(tau_eff, tau_matrix):
    """
    :param tau_eff: (numpy array) tau_eff values from actual composite
    :param tau_matrix: (numpy array) tau_eff vales from bootstrap

    :return:
        sigma: (numpy array) 1D errors for tau_eff
    """
    import numpy as np

    N = tau_matrix.shape[0]
    diff_matrix = np.array([tau_matrix[i] - tau_eff for i in range(N)])

    tranpose = np.transpose(diff_matrix)
    covariance = (1/(N-1))*(np.dot(tranpose, diff_matrix))
    #sigma = np.sqrt(np.diagonal(covariance))

    return covariance


def kirkman_correction(tau, red, plot=False, p=False):
    """
    :param tau: (numpy array) the tau_eff values from both bins
    :param red: (numpy array) the redshift values from both bins

    :return:
        (numpy array) the corrected tau values
    """
    import numpy as np
    import matplotlib.pyplot as plt

    lam = 1216 * ((1 + red) / (1 + np.median(red)))  # using the pivot to return to angstroms

    tau_m = np.log(1 + (0.01564 - (4.646 * 10 ** (-5) * (lam - 1360))))

    if p == True:
        print("the average percentage contributed by metals:", 100 * np.average(tau_m / tau))

    if plot == True:
        plt.scatter(red, tau, label="original")
        plt.scatter(red, tau - tau_m, label="corrected")
        plt.legend()
        plt.show()

    return tau - tau_m
