def xspec(z_bin):

    """
    :param z_bin: array of min and max redshift(z) values of the desired bin ex:[2.0,2.5]
    :return: array of XSpectrum1D objects and corresponding redshifts within the specified z_bin

    """
    import numpy as np
    from astropy.io import fits
    from astropy.table import Table as Table
    import matplotlib.pyplot as plt
    from linetools.spectra import utils as ltsu
    from linetools.spectra.xspectrum1d import XSpectrum1D
    import astropy.units as u
    from astropy import constants as const

    path_16 = "/home/jsm/PycharmProjects/tau_eff/Spectra/2016/"

    path_17 = "/home/jsm/PycharmProjects/tau_eff/Spectra/2017/"


    spec_atr_16 = Table.read(path_16 + "spec_atr.txt", format = 'ascii')

    spec_atr_17 = Table.read(path_17 + "spec_atr.txt", format = 'ascii')

#read in

    init_spec_16 = []

    init_z_16 = []

    for entry in spec_atr_16: # from the CLAMATO 2016 survey

        if np.min(z_bin) < entry["zspec"] < np.max(z_bin):  # creating the bin size

            if entry["Conf"] < 10.0:  # excluding any QSOs

                temp = XSpectrum1D.from_file(path_16 + entry["Filename"])

                if temp.wvmin < (1216 * u.AA) * (1 + entry["zspec"]) < temp.wvmax:

                    init_z_16.append(entry["zspec"])

                    init_spec_16.append(XSpectrum1D.from_file(path_16 + entry["Filename"]))


    init_spec_17 = []

    init_z_17 = []

    for entry in spec_atr_17:

        if np.min(z_bin) < entry["col5"] < np.max(z_bin):

            if entry["col4"] < 10.0:

                temp = XSpectrum1D.from_file(path_17 + entry["col1"])

                if temp.wvmin < (1120 * u.AA) * (1 + entry["col5"]) < temp.wvmax:

                    init_z_17.append(entry["col5"])

                    init_spec_17.append(XSpectrum1D.from_file(path_17 + entry["col1"]))

#array creation

    spec = np.asarray(init_spec_16 + init_spec_17)

    print("Number of spectra (Nspec) in the redshift bin", len(spec))

    red = np.asarray(init_z_16 + init_z_17)

#a simple histogram to illustrate the sample selected


    z_16 = [entry["zspec"] for entry in spec_atr_16]  # from the CLAMATO 2016 survey

    z_17 = [entry["col5"] for entry in spec_atr_17]  # from the CLAMATO 2017 survey

    z_tot = z_16 + z_17

#the plot

    plt.hist(z_tot, bins=30, edgecolor='white', linewidth=1.2, label="Full Sample", color="#f03b20")
    plt.hist(red, bins=10, edgecolor='white', linewidth=1.2, label="Reduced Sample", color="#feb24c")
    plt.ylabel("Number of Galaxies")
    plt.xlabel("Spectroscopic Redshift $z$")
    plt.xlim(1.5,3.5)
    plt.legend()
    plt.show()

    return spec, red

