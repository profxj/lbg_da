def stack_and_strap(spec , zem , wv_norm, dv, N, ):

    """

    :param spec: array of XSpectrum1D objects
    :param zem: array of their redshift values (z)
    :param wv_norm: array of the wavelength range (angstroms) in the rest frame in which to normalize, ex:[1450,1500]
    :param dv: velocity width of new pixels for rebin (will be measured in km/s)
    :param N: number of iterations fro the bootstrap
    :return: a single composite XSpectrum1D object with bootstrap generated errors


    """

    import numpy as np
    from linetools.spectra import utils as ltsu
    from linetools.spectra.xspectrum1d import XSpectrum1D
    import astropy.units as u


    print("Number of spectra provided =", len(spec))

    ## the normalization
    r_2 = range(len(spec))

    temp_2 = [np.asarray(spec[i].wavelength / (1 + zem[i])) for i in r_2]

    wave_norm = [(wv_norm[0] < entry) & (entry < wv_norm[1]) for entry in temp_2]

    rough_med = np.asarray([np.median(spec[i].flux[wave_norm[i]]) for i in r_2])

    g_med = [median > 0.0 for median in rough_med] #to ensure a solid rebin

    medians = rough_med[g_med]

    g_spec = spec[g_med]

    g_red = zem[g_med]

    r_3 = range(len(g_spec))

    norm_flux = np.asarray([(g_spec[i].flux / medians[i]) for i in r_3])

    ## the new scaled Xspec objects

    scaled_spec = []

    for i in r_3:

        scaled_spec.append(XSpectrum1D(g_spec[i].wavelength, norm_flux[i], sig=g_spec[i].sig))

    ## the first trim on the wavelegth array to ensure a solid stack

    trim_spec = []

    for i in range(len(scaled_spec)):

        trim_spec.append(XSpectrum1D(scaled_spec[i].data["wave"][0][440:1500],
                                     scaled_spec[i].data["flux"][0][440:1500],
                                     scaled_spec[i].data["sig"][0][440:1500]))

    collate = ltsu.collate(trim_spec)

    print("Number of spectra to be stacked =", collate.nspec)

    ## rest-frame wave values and smash

    rest_spec = ltsu.rebin_to_rest(collate, g_red, dv*u.km / u.s, grow_bad_sig=True)

    stack = ltsu.smash_spectra(rest_spec)

    ## the bootstrap

    R = (len(trim_spec))  # restraints on the random variable

    M = len(trim_spec[0].flux)  # columns

    N_stack = []

    for i in np.arange(N): #this should take a while!

        choice = np.asarray(ran.randint(0, R, R))  # the shuffle

        rand_spec = np.asarray([rest_spec[index] for index in choice])

        rand_collate = ltsu.collate(rand_spec)  # collate and stack

        N_stack.append(ltsu.smash_spectra(rand_collate))


    N_flux = np.array([entry.flux for entry in N_stack])  # getting into a np.array for matrix math

    N_matrix = np.array([N_flux[i] - stack.flux for i in range(N)])

    tranpose = np.transpose(N_matrix)  # matrix math

    covariance = np.dot(tranpose, N_matrix)

    sigma = np.sqrt(np.diagonal(covariance) / (N - 1))  # the final error array

    composite = XSpectrum1D(stack.wavelength, stack.flux, sig=sigma)  # the full composite

    return composite
