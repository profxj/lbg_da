def lbg_smash(spec, red, ebv, s2n, dv):
    """
    :param spec: array of XSpectrum1D objects: use clamato_read.py
    :param red: array of their redshift values (z)
    :param ebv: array of E(B-V) values from attenuation caused by the Milky Way
    :param s2n: float lowest value of S/N. will be computed over wv_norm
    :param dv: velocity width of new pixels for rebin (will be measured in km/s)
    :return: a single composite XSpectrum1D object only flux & wavelength values.

    """

    import numpy as np
    from numpy.lib.polynomial import poly1d

    import matplotlib.pyplot as plt

    import astropy.units as u

    from linetools.spectra import utils as ltsu
    from linetools.spectra.xspectrum1d import XSpectrum1D



    # S/N try 2.5
    # dv try 300

    # Note: this routine normalizes spectra individually before stacking them,

    print("Number of spectra provided =", len(spec))

    # initial S/N cut

    r = range(len(spec))

    temp = [np.asarray(spec[i].wavelength / (1 + red[i])) for i in r]  # rest-frame

    norm = [(1260 < entry) & (entry < 1304) for entry in temp]

    signal = np.asarray([spec[i].flux[norm[i]] for i in r])

    noise = np.asarray([spec[i].sig[norm[i]] for i in r])

    s2n_array = np.asarray([np.median(signal[i] / noise[i]) for i in r])

    s2n_cut = [i > s2n for i in s2n_array]

    spec = spec[s2n_cut]

    plt.scatter(red, s2n_array, label="Full sample", color="#f03b20")
    plt.scatter(red[s2n_cut], s2n_array[s2n_cut], label="Reduced sample", color="#feb24c")
    plt.ylabel("S/N")
    plt.xlabel("Spectroscopic Redshift $z$")
    plt.legend()
    plt.show()

    red = red[s2n_cut]

    print("S/N cut left", len(spec), "spectra")

    # milky way extinction

    r = range(len(spec))

    flux = []

    for i in r:

        x = 10000. / np.array(spec[i].wavelength)  # Convert to inverse microns
        npts = x.size
        a = np.zeros(npts, dtype=np.float)
        b = np.zeros(npts, dtype=np.float)
        r_v = 3.1

        good = np.where((x >= 0.3) & (x < 1.1))
        if len(good[0]) > 0:
            a[good] = 0.574 * x[good] ** (1.61)
            b[good] = -0.527 * x[good] ** (1.61)

        good = np.where((x >= 1.1) & (x < 3.3))
        if len(good[0]) > 0:  # Use new constants from O'Donnell (1994)
            y = x[good] - 1.82

            c1 = np.array([1., 0.104, -0.609, 0.701, 1.137, -1.718, -0.827, 1.647, -0.505])  # from O'Donnell
            c2 = np.array([0., 1.952, 2.908, -3.989, -7.985, 11.102, 5.491, -10.805, 3.347])

            a[good] = poly1d(c1[::-1])(y)
            b[good] = poly1d(c2[::-1])(y)

        good = np.where((x >= 3.3) & (x < 8))
        if len(good[0]) > 0:
            y = x[good]

            a[good] = 1.752 - 0.316 * y - (0.104 / ((y - 4.67) ** 2 + 0.341))  # + f_a
            b[good] = -3.090 + 1.825 * y + (1.206 / ((y - 4.62) ** 2 + 0.263))  # + f_b

        good = np.where((x >= 8) & (x <= 11))
        if len(good[0]) > 0:
            y = x[good] - 8.

            c1 = np.array([-1.073, -0.628, 0.137, -0.070])
            c2 = np.array([13.670, 4.257, -0.420, 0.374])
            a[good] = poly1d(c1[::-1])(y)
            b[good] = poly1d(c2[::-1])(y)

        # Now apply extinction correction to input flux vector

        a_v = r_v * ebv[i]

        a_lambda = a_v * (a + b / r_v)

        funred = spec[i].flux * 10. ** (0.4 * a_lambda)  # Derive unreddened flux

        flux.append(funred)

    # the normalization

    temp = [np.asarray(spec[i].wavelength / (1 + red[i])) for i in r]

    wv_norm = [(1260 < entry) & (entry < 1304) for entry in temp]  # just between the SiII lines.

    flux_range = np.asarray([flux[i][wv_norm[i]] for i in r])

    medians = np.asarray([np.median(flux_range[i]) for i in r])

    norm_flux = np.asarray([(flux[i] / medians[i]) for i in r])

    norm_spec = [XSpectrum1D(spec[i].wavelength, norm_flux[i], spec[i].sig) for i in r]  # the normalized spec array

    # truncating each spectra to only include wavelength values from 950A to 1500A

    r = range(len(norm_spec))

    temp = np.array([norm_spec[i].wavelength / (1 + red[i]) for i in r])

    trim = []

    for i in r:
        trim.append([950 * u.AA <= temp[i][j] <= 1500 * u.AA for j in range(len(temp[i]))])
        # to include only enough to extrapolate from

    t_spec = []

    for i in r:
        t_spec.append(XSpectrum1D(norm_spec[i].wavelength[trim[i]][0:935],
                                  norm_spec[i].flux[trim[i]][0:935],
                                  norm_spec[i].flux[trim[i]][0:935]))

    collate = ltsu.collate(t_spec)  # to get a single XSpec object

    print("Number of spectra to be stacked =", collate.nspec)

    fix = []  # to make sure none of the values included will mess up the stack

    a = collate.data['wave']

    for i in range(collate.nspec):

        if np.min(a[i]) == 0.0:
            fix.append(i)

    if len(fix) > 0:
        print(fix)  # exclude all spec[fix]

    # rebin to rest frame with dv dispersion

    rest_spec = ltsu.rebin_to_rest(collate, red, dv * u.km / u.s, grow_bad_sig=True)

    # the stack! Note: the error array was not normalized; do not trust it

    stack = ltsu.smash_spectra(rest_spec)

    stack.plot()

    return stack, rest_spec