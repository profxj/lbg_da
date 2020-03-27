def stacks(spec, red, N, zbin, plot=False):
    """
    :param spec: (numpy array) XSpectrum1D objects
    :param red: (numpy array) their redshift values (z)
    :param N: (int) number of bootstrap iterations
    :param zbin: (str) redshift bin "lowz" or "hiz"
    """

    import numpy as np
    from astropy.io import fits
    from numpy import random as ran
    from linetools.spectra import utils as ltsu
    from linetools.spectra.xspectrum1d import XSpectrum1D

    print("Number of spectra provided =", spec.nspec)

    # the actual stack
    stack = ltsu.smash_spectra(spec)
    z_med = np.median(red)

    # the bootstrap
    R = (spec.nspec)  # restraints on the random variable
    N_stack = []
    N_z_med = []

    for i in np.arange(N):  # this might take a while

        if i % 10 == 0:
            print(i)
        choice = np.asarray(ran.randint(0, R, R))  # the shuffle

        rand_spec = np.asarray([spec[index] for index in choice])
        rand_red = np.asarray([red[index] for index in choice])

        rand_collate = ltsu.collate(rand_spec)  # collate and stack
        N_stack.append(ltsu.smash_spectra(rand_collate))
        N_z_med.append(np.median(rand_red))

    # matrix math to create the error array
    N_flux = np.array([entry.flux for entry in N_stack])
    N_matrix = np.array([N_flux[i] - stack.flux for i in range(N)])

    tranpose = np.transpose(N_matrix)  # matrix math
    covariance = np.dot(tranpose, N_matrix)
    sigma = np.sqrt(np.diagonal(covariance) / (N - 1))  # the final error array

    composite = XSpectrum1D(stack.wavelength, stack.flux, sig=sigma)
    # plot with bootstrap generated error
    if plot == True:
        composite.plot()

    # writing out the stack
    hdr = fits.Header()
    hdr["REDSHIFT"] = z_med
    hdr["NSPEC"] = spec.nspec
    header = fits.PrimaryHDU(header=hdr)

    c_wave = fits.Column(name='WAVELENGTH', array=stack.wavelength, unit="angstroms", format="E")
    c_flux = fits.Column(name='FLUX', array=stack.flux, unit="relative flux", format="E")
    c_noise = fits.Column(name='FLUX_ERR', array=sigma, unit="relative flux", format="E")
    table = fits.BinTableHDU.from_columns([c_wave, c_flux, c_noise])

    full_hdu = fits.HDUList([header, table])
    full_hdu.writeto("../../fits/composites/" + zbin + "/composite.fits")

    # writing out the boostrap
    hdr = fits.Header()
    hdr["NSPEC"] = spec.nspec
    hdr["NBOOT"] = N
    header = fits.PrimaryHDU(header=hdr)

    c_wave = fits.Column(name='WAVELENGTH', array=stack.wavelength, unit="angstroms", format="E")
    c_noise = fits.Column(name='FLUX_ERR', array=sigma, unit="relative flux", format="E")
    table = fits.BinTableHDU.from_columns([c_wave, c_noise])

    flux_data = fits.ImageHDU(N_flux, name="FLUX_ITERATIONS")

    full_hdu = fits.HDUList([header, table, flux_data])
    full_hdu.writeto("../../fits/bootstrap/" + zbin + "/stacks.fits")