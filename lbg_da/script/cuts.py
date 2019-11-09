def lyman_limit(spec, red, coord, thres):
    
    """
    :param spec: [array] of XSpectrum1D objects: use clamato_read.py
    :param red: [array] of their redshift values (z)
    :param coord: [array] of coordinates
    :param thres: [int] value for the maximum std to include in the lyman limit cut
    """
    
    import numpy as np
    import matplotlib.pyplot as plt
    from astropy import units as u


    print("Number of spectra provided =", len(spec))

    r = range(len(spec))
    
    lylim_med = []
    
    for i in r:
        
        rest_wave = spec[i].wavelength / (1 + red[i]) # rest-frame
                
        sig = np.array(spec[i].sig[rest_wave < 912*u.AA])
        
        flux = np.array(spec[i].flux[rest_wave < 912*u.AA])

        if len(sig) > 0:
                
            flux = flux[sig > 0.0]
            
            med = np.median(flux)
            
            lylim_med.append(med)
            
        else:
            
            lylim_med.append(np.nan)
            
            
    lylim_med = np.asarray(lylim_med)
    
    keep = np.isnan(lylim_med)
    cut = np.invert(keep)
            
             
    nan_spec = spec[keep] #to keep the rest of the spectra that dont reach
    nan_red = red[keep]
    nan_coord = coord[keep]
    
    lylim_spec = spec[cut] #these are the spectra to cut    
    lylim_red = red[cut]
    lylim_coord = coord[cut]
    
    lylim_dist = lylim_med[cut] #the stats used to cut
    mean = np.mean(lylim_dist)
    std = np.std(lylim_dist)
    
    upper = (mean + (thres*std))
    lower = (mean - (thres*std))
    
    lylim_cut = [(lower < i) & (i < upper) for i in lylim_dist]
    
    #plot!
    
    plt.figure(figsize=(10,10))
    
    plt.hist(lylim_dist, bins = 60,edgecolor='white', linewidth=1.2, color="grey")
    
    plt.vlines((mean - (thres*std)),0,100, linestyle="--", color="green", label = "1 $\sigma$")
    plt.vlines((mean + (thres*std)),0,100, linestyle="--", color="green")
    
    #plt.vlines((mean - (thres*std*2)),0,80, linestyle="--", color="orange", label = "2 $\sigma$")
    #plt.vlines((mean + (thres*std*2)),0,80, linestyle="--", color="orange")
    
    #plt.vlines((mean - (thres*std*3)),0,60, linestyle="--", color="red", label = "3 $\sigma$")
    #plt.vlines((mean + (thres*std*3)),0,60, linestyle="--", color="red")
    
    plt.vlines(mean,0,120, linestyle="--", color="black", label = "mean")
        
    plt.ylabel("Number of Spectra",fontsize=15)
    plt.xlabel("Median flux values blueward of 912$\AA$",fontsize=15)
    
    plt.legend(fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlim(-1,1)
    
    #plt.savefig("/Users/jsmonzon/IGM-UCSC/figures/ll_cut.pdf")
    
    plt.show()

    
    #cut!
    
    spec = np.append(nan_spec, lylim_spec[lylim_cut])

    red = np.append(nan_red, lylim_red[lylim_cut])
    
    coord = (list(nan_coord)) + (list(lylim_coord[lylim_cut,:])) # a little patch to add the two arrays
    
    coord = np.asarray(coord)
        
    print("Number of spectra left after the lyman-limit cut = ", len(spec))
    
    return spec, red, coord
    
    
def s2n(spec, red, coord, s2n):
    
    """
    :param spec: [array] of XSpectrum1D objects: use clamato_read.py
    :param red: [array] of their redshift values (z)
    :param coord: [array] of coordinates
    :param s2n: [float] lowest value of S/N to include in the cut
    """
    
    import numpy as np
    import matplotlib.pyplot as plt

    # now we exlcude any galaxies below the s2n

    r = range(len(spec))

    temp = [np.asarray(spec[i].wavelength / (1 + red[i])) for i in r]  # rest-frame

    good = [(1260 < entry) & (entry < 1304) for entry in temp] # masking to find a good region

    signal = np.asarray([spec[i].flux[good[i]] for i in r])

    noise = np.asarray([spec[i].sig[good[i]] for i in r])

    s2n_array = np.asarray([np.median(signal[i] / noise[i]) for i in r])

    s2n_cut = [i > s2n for i in s2n_array]
    
    #plot!
    
    plt.figure(figsize=(10,10))
    
    plt.scatter(red, s2n_array, label="Excluded galaxies", color="#f03b20")
    plt.scatter(red[s2n_cut], s2n_array[s2n_cut], label="Reduced sample", color="#fd8d3c")
                
    plt.hlines(s2n,min(red),max(red), linestyle="--", color="grey",label="S/N cutoff")
    plt.vlines(2.5,0,9,linestyle="--", color="grey",)
    plt.ylabel("S/N",fontsize=15)
    plt.xlabel("Spectroscopic Redshift $z$",fontsize=15)
    plt.ylim(0,9)
    
    plt.legend(fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    
    #plt.savefig("/Users/jsmonzon/IGM-UCSC/figures/noise_scatter.pdf")
    plt.show()


    
    #cut!
    
    spec = spec[s2n_cut]
    
    red = red[s2n_cut]
    
    coord = coord[s2n_cut,:]
    
    print("Number of spectra left after the S/N cut = ", len(spec))
    
    return spec, red, coord