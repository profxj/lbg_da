def covariance_error(spec, N, stack):

    """

    :param spec: XSpectrum1D collated object :rest_spec from the stacking modules
    :param N: number of iterations for the bootstrap
    :param stack: XSpectrum1D stack that you are calculating the error of
    :return: error array to be combined with stacked spectra


    """

    import numpy as np
    from linetools.spectra import utils as ltsu
    from linetools.spectra.xspectrum1d import XSpectrum1D
    import astropy.units as u
    from numpy import random as ran

    print("Number of spectra provided =", spec.nspec)

    ## the bootstrap

    R = (spec.nspec)  # restraints on the random variable

    N_stack = []

    for i in np.arange(N): #this should take a while!

        choice = np.asarray(ran.randint(0, R, R))  # the shuffle

        rand_spec = np.asarray([spec[index] for index in choice])

        rand_collate = ltsu.collate(rand_spec)  # collate and stack

        N_stack.append(ltsu.smash_spectra(rand_collate))


    N_flux = np.array([entry.flux for entry in N_stack])  # getting into a np.array for matrix math

    N_matrix = np.array([N_flux[i] - stack.flux for i in range(N)])

    tranpose = np.transpose(N_matrix)  # matrix math

    covariance = np.dot(tranpose, N_matrix)

    sigma = np.sqrt(np.diagonal(covariance) / (N - 1))  # the final error array

    composite = XSpectrum1D(stack.wavelength,stack.flux, sig = sigma)

    return composite
