def extinction(spec, red, coord):
    
    """
    :param spec: [array] of XSpectrum1D objects: use clamato_read.py
    :param red: [array] of their redshift values (z)
    :param coord: [array] of coordinates
    :return: [array] of unreddened Xspectrum1D objects
    """
    
    import numpy as np
    from numpy.lib.polynomial import poly1d
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from dustmaps.bayestar import BayestarQuery
    from astropy.cosmology import WMAP9 as cosmo
    from linetools.spectra.xspectrum1d import XSpectrum1D

    
    r = range(len(spec))

    Mpc = cosmo.comoving_distance(red)
    bayestar = BayestarQuery()
    coords = [SkyCoord(coord[i][0] * u.deg, coord[i][1] * u.deg, distance=Mpc[i], frame='fk5') for i in r]
    ebv = [bayestar(i) for i in coords] #to get the ebv values for each galaxy
    
    unred_spec = []

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
        funred = np.asarray(funred)
        unred_spec.append(XSpectrum1D(spec[i].wavelength, funred, spec[i].sig))
        
    return np.asarray(unred_spec)