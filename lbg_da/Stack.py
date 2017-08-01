
def stack(spec_atr,basestring, s2n,dv):
    """Stack XSpec files from a table of attributes

    Parameters
    ----------
    spec_atr : table(ascii), includes spectral files and their Zem values
    basestring : str, the path to the spectral files
    s2n : float, the minimum value of s/n
    dv : float, velocity width of new pixels for rebin will be measured in km/s
    """

    import numpy as np
    from astropy.io import fits
    from astropy.table import Table
    import matplotlib.pyplot as plt
    from linetools.spectra import utils as ltsu
    from linetools.spectra.xspectrum1d import XSpectrum1D
    import astropy.units as u
    from array import array
    import pdb

    # initial filter on speclist to ensure a successful rebin

    good_spec = []

    good_z = []

    for entry in spec_atr:

        if entry["zspec"] > 0.00:  # for non zero redshifts

            temp = XSpectrum1D.from_file(basestring + entry['Filename'])

            if temp.wvmin < (1216*u.AA) * (1+entry["zspec"]) < temp.wvmax:  # for lya coverage

                good_z.append(entry['zspec'])

                good_spec.append(XSpectrum1D.from_file(basestring + entry['Filename']))


    # first rebin to calculate s/n


    coll_1 = ltsu.collate(good_spec)  # to create a single xspec object

    reb_1 = ltsu.rebin_to_rest(coll_1, good_z, dv * u.km / u.s, grow_bad_sig=True)

    wv = np.median(reb_1.wavelength)  # median near lya

    s2n_list = [XSpectrum1D.get_local_s2n(i, wv)[0] for i in reb_1]  # index zero so that we exclude the deviation

    # 2nd rebin on s/n filter

    trim = np.asarray([good_spec, s2n_list]).T  # an array to match the spectra with their S/N values

    viable = trim[trim[:, 1] > s2n][:, 0]  # to filter the spectra according to s/n

    coll_2 = ltsu.collate(viable)  # for the new Xspec object

    new_z = np.asarray(good_z)[trim[:, 1] > s2n].tolist()  # to match the new list with Zem values

    reb_2 = ltsu.rebin_to_rest(coll_2, new_z, dv * u.km / u.s, grow_bad_sig=True)  # another rebin

    #the stack

    smash = ltsu.smash_spectra(reb_2)  # the stack

    outfile = str(input("enter the name of the file to be saved to disk:" ))

    smash.write_to_fits(outfile)






