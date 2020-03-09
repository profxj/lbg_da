def lyman_limit(spec, red, coord, plot=False):
    
    """
    :param spec: [array] of XSpectrum1D objects: use clamato_read.py
    :param red: [array] of their redshift values (z)
    :param coord: [array] of coordinates
    """
    
    import numpy as np
    import matplotlib.pyplot as plt
    from astropy import units as u

    print("Number of spectra provided =", len(spec))

    r = range(len(spec))

    ll_meds = []
    ll_cut = []

    for i in r:

        rest_wave = spec[i].wavelength / (1 + red[i])  # rest-frame
        sig = np.array(spec[i].sig[rest_wave < 912 * u.AA])
        flux = np.array(spec[i].flux[rest_wave < 912 * u.AA])

        try:  # take the median

            ll_flux = np.sort(flux[sig > 0.0])
            n = len(ll_flux)

            if n % 2 == 0:
                med1 = ll_flux[n // 2]
                med2 = ll_flux[n // 2 - 1]
                median = (med1 + med2) / 2
            else:
                median = ll_flux[n // 2]

            ll_meds.append(median)

            if (-.2 < median) & (median < .2):
                ll_cut.append(i)

        except IndexError:
            ll_cut.append(i)

    ll_meds = np.asarray(ll_meds)
    ll_cut = np.asarray(ll_cut)

    # the plot!

    if plot == True:

        plt.figure(figsize=(8, 8))
        plt.hist(ll_meds[ll_meds < 1.5], bins=30, edgecolor="black",
             linewidth=1.2, label="Full Sample", color="grey")

        plt.hist(ll_meds[(-.2 < ll_meds) & (ll_meds < .2)], bins=6, edgecolor="black",
             linewidth=1.2, label="Reduced Sample", color="#fd8d3c")

        plt.xlabel("median flux value past 912A", fontsize=15)
        plt.ylabel("number of spectra", fontsize=15)

        plt.legend(fontsize=12)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.show()

    #the cut!
    spec = spec[ll_cut]
    red = red[ll_cut]
    coord = coord[ll_cut]

    print("Number of spectra left after the lyman-limit cut = ", len(spec))

    return spec, red, coord
    
def s2n(spec, red, coord, plot=False):
    
    """
    :param spec: [array] of XSpectrum1D objects: use clamato_read.py
    :param red: [array] of their redshift values (z)
    :param coord: [array] of coordinates
    :param s2n: [float] lowest value of S/N to include in the cut
    """
    
    import numpy as np
    import matplotlib.pyplot as plt

    # now we exlcude any galaxies below the s2n

    # now we exlcude any galaxies below the s2n
    print("Number of spectra provided =", len(spec))

    r = range(len(spec))

    s2n = 1.5

    temp = [np.asarray(spec[i].wavelength / (1 + red[i])) for i in r]  # rest-frame
    good = [(1260 < entry) & (entry < 1304) for entry in temp]  # masking to find a good region
    signal = np.asarray([spec[i].flux[good[i]] for i in r])
    noise = np.asarray([spec[i].sig[good[i]] for i in r])
    s2n_array = np.asarray([np.median(signal[i] / noise[i]) for i in r])
    s2n_cut = [i > s2n for i in s2n_array]
    
    #plot!
    if plot == True:

        plt.figure(figsize=(10, 10))

        plt.scatter(red, s2n_array, label="Excluded spectra", color="#f03b20")
        plt.scatter(red[s2n_cut], s2n_array[s2n_cut], label="Reduced sample", color="#fd8d3c")

        plt.hlines(s2n, 2.2, 2.8, linestyle="--", color="black", label="S/N cutoff")
        plt.ylabel("median S/N value ($1260 - 1304\AA$)", fontsize=15)
        plt.xlabel("$z$", fontsize=15)
        plt.ylim(0, 9)

        plt.xlim(2.2, 2.8)
        plt.legend(fontsize=12)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)

        plt.show()

    # cut!
    spec = spec[s2n_cut]
    red = red[s2n_cut]
    coord = coord[s2n_cut, :]

    print("Number of spectra left after the S/N cut = ", len(spec))
    
    return spec, red, coord