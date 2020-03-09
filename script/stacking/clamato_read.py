def xspec(z_bin, plot=False):
    
    """
    :param z_bin: array of min and max redshift(z) values of the desired bin ex:[2.0,2.5]
    :return: array of XSpectrum1D objects, corresponding redshifts and coordiantes

    """

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.table import Table
    import astropy.units as u
    from linetools.spectra.xspectrum1d import XSpectrum1D

    path_16 = "../../spectra/2016/"
    path_17 = "../../spectra/2017/"

    spec_atr_16 = Table.read(path_16 + "spec_atr.txt", format='ascii')
    spec_atr_17 = Table.read(path_17 + "spec_atr.txt", format='ascii')

    # read in

    spec_16 = []
    z_16 = []
    coord_16 = []

    for entry in spec_atr_16:  # from the CLAMATO 2016 survey

        if np.min(z_bin) <= entry["zspec"] <= np.max(z_bin):  # creating the bin size

            if entry["Conf"] < 10.0:  # excluding any QSOs

                temp = XSpectrum1D.from_file(path_16 + entry["Filename"])

                if temp.wvmin < (1216 * u.AA) * (1 + entry["zspec"]) < temp.wvmax:

                    coord_16.append([entry["RA"], entry["Dec"]])  # coordinates in deg

                    z_16.append(entry["zspec"])

                    spec_16.append(XSpectrum1D.from_file(path_16 + entry["Filename"]))

    spec_17 = []
    z_17 = []
    coord_17 = []

    for entry in spec_atr_17:

        if np.min(z_bin) <= entry["col5"] <= np.max(z_bin):

            if entry["col4"] < 10.0:

                temp = XSpectrum1D.from_file(path_17 + entry["col1"])

                if temp.wvmin < (1216 * u.AA) * (1 + entry["col5"]) < temp.wvmax:

                    coord_17.append([entry["col7"], entry["col8"]])

                    z_17.append(entry["col5"])

                    spec_17.append(XSpectrum1D.from_file(path_17 + entry["col1"]))

    # array creation

    spec = np.asarray(spec_16 + spec_17)
    red = np.asarray(z_16 + z_17)
    coord = np.asarray(coord_16 + coord_17)

    print("Number of spectra in the redshift bin:", len(spec))

    z_16_tot = [entry["zspec"] for entry in spec_atr_16]  # from the CLAMATO 2016 survey
    z_17_tot = [entry["col5"] for entry in spec_atr_17]  # from the CLAMATO 2017 survey
    z_tot = z_16_tot + z_17_tot

    # a histogram to illustrate the sample selected
    
    if plot == True:

        plt.figure(figsize=(8, 8))

        plt.hist(z_tot, bins=70, edgecolor='white', linewidth=1.2, label="Full Sample", color="grey")
        plt.hist(red, bins=12, edgecolor='black', linewidth=1.2, label="Reduced Sample", color="#fd8d3c")
    
        plt.ylabel("Number of Galaxies",fontsize=15)
        plt.xlabel("$z$",fontsize=15)
    
        plt.xlim(1.5, 3.5)
    
        plt.legend(fontsize=12)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.show()
    
    return spec, red, coord
