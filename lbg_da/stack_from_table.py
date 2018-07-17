def stack_from_table(spec_atr, basestring, s2n_min, dv):


    """Stack XSpec files from a table of attributes
    Parameters
    ----------
    spec_atr : table(ascii), includes spectra filenames and their Zem values
    basestring : str, the path to the spectra filenames
    s2n : float, the minimum value of s/n
    dv : float, velocity width of new pixels for rebin will be measured in km/s
    """

    import numpy as np
    from astropy.io import fits
    from astropy.table import Table as Table
    import matplotlib.pyplot as plt
    from linetools.spectra import utils as ltsu
    from linetools.spectra.xspectrum1d import XSpectrum1D
    import astropy.units as u
    from array import array
    from astropy import constants as const


    # an initial filter to ensure the spectra included cover lya

    init_spec = []

    init_z = []

    for entry in spec_atr:

        if 2.00 < entry["zspec"] < 2.9:

            temp = XSpectrum1D.from_file(basestring + entry["Filename"])

            if temp.wvmin < (1216 * u.AA) * (1 + entry["zspec"]) < temp.wvmax:
                init_z.append(entry["zspec"])

                init_spec.append(XSpectrum1D.from_file(basestring + entry["Filename"]))

    init_spec = np.asarray(init_spec)

    init_z = np.asarray(init_z)

    print("Number of spectra (Nspec) after 1st filter =", len(init_spec))


    # second filter to ensure that I have full coverage of the values I want

    wave_range_mask = [(3600 * u.AA > spec.wvmin) and (4600 * u.AA < spec.wvmax) for spec in init_spec]

    speclist = init_spec[wave_range_mask]

    z_val = init_z[wave_range_mask]

    print("Nspec after 2nd filter=", len(speclist))


    # now to calculate the local S/N and use it as the third filter

    wave = np.linspace(3600, 4600)

    s2n = []

    for i in range(len(speclist)):

        s2n.append([XSpectrum1D.get_local_s2n(speclist[i], wavelength * u.AA)[0] for wavelength in wave])

    s2n_mean = []

    s2n_med = []

    for i in s2n:

        s2n_mean.append(np.mean(i))

        s2n_med.append(np.median(i))

    # an array to match the spectra with their s2n values
    trim = np.asarray([speclist, s2n_mean]).T

    # grabbing only those above the min value
    trimmed_spec = trim[trim[:, 1] > s2n_min][:, 0]

    trimmed_z = z_val[trim[:, 1] > s2n_min]

    print("Nspec after 3rd filter=", len(trimmed_spec))


    # normalizing the flux
    flux_med = [np.median(trimmed_spec[i].flux) for i in range(len(trimmed_spec))]

    scaled_flux = np.asarray([(trimmed_spec[i].flux / flux_med[i]) for i in range(len(trimmed_spec))])

    # print(len(scaled_flux[0]))

    scaled_spec = []

    # the new scaled Xspec objects
    for i in range(len(trimmed_spec)):

        scaled_spec.append(XSpectrum1D(trimmed_spec[i].wavelength, scaled_flux[i], sig=trimmed_spec[i].sig))


    # a final trim on the wavelegth array to ensure a solid stack
    new_spec = []

    for i in range(len(scaled_spec)):
        new_spec.append(XSpectrum1D(scaled_spec[i].data["wave"][0][640:1460],
                                    scaled_spec[i].data["flux"][0][640:1460],
                                    scaled_spec[i].data["sig"][0][640:1460]))

    collate = ltsu.collate(new_spec)

    # rest frame wave values
    rest_spec = ltsu.rebin_to_rest(collate, trimmed_z, dv * u.km/u.s, grow_bad_sig=True)

    # the stack!
    stack = ltsu.smash_spectra(rest_spec)

    return stack
