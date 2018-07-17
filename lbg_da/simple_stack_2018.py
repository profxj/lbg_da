def simple_stack(spec , zem , wave_norm, dv):

    """

    :param spec: array of XSpectrum1D objects
    :param zem: array of their redshift values (z), they cannot exceed a redshift value of 3.00
    :param wave_norm: list of the wavelength range in the rest frame from which to normalize, ex: [1450,1500]
    :param dv: velocity width of new pixels for rebin measured in km/s
    :return: a single composite XSpectrum1D object


    """

    print("Number of spectra provided =", len(spec))


    wave_range_mask = [(4000*u.AA >= spec.wvmin) & (5600*u.AA <= spec.wvmax) for spec in spec]
    # his masks ensures that viableparts of the spectrum are included

    speclist = spec[wave_range_mask]

    z_val = zem[wave_range_mask]


    print("Number of spectra left after initial filter =", len(speclist))


    temp = [np.asarray(speclist[i].wavelength / (1 + z_val[i])) for i in range(len(speclist))]
    # Rebining so that waverange can be used in the rest frame

    good_wave = [(wave_norm[0] < entry) & (entry < wave_norm[1]) for entry in temp]

    raw_med = np.asarray([np.median(speclist[i].flux[good_wave[i]]) for i in range(len(speclist))])


    good_med = [median > 0.0 for median in raw_med]
    # to ensure no Nan values are created

    medians = raw_med[good_med]

    short = speclist[good_med]

    z_short = z_val[good_med]


    scaled_flux = np.asarray([(short[i].flux / medians[i]) for i in range(len(short))])

    scaled_spec = []

    # the new scaled Xspec objects
    for i in range(len(short)):

        scaled_spec.append(XSpectrum1D(short[i].wavelength, scaled_flux[i], sig=short[i].sig))

    # the final trim on the wavelength array to ensure a solid stack
    new_spec = []

    for i in range(len(scaled_spec)):
        new_spec.append(XSpectrum1D(scaled_spec[i].data["wave"][0][240:1460],
                                    scaled_spec[i].data["flux"][0][240:1460],
                                    scaled_spec[i].data["sig"][0][240:1460]))

    collate = ltsu.collate(new_spec)

    print ("The final number of spectra stacked =", collate.nspec)

    # rest frame wave values
    rest_spec = ltsu.rebin_to_rest(collate, z_short, dv*u.km / u.s, grow_bad_sig=True)

    stack = ltsu.smash_spectra(rest_spec)

    return stack

