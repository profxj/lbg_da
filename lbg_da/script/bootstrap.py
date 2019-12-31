def stacks(spec, red, N, zbin):
    """
    :param spec: XSpectrum1D collated object :rest_spec from the stacking modules
    :red: redshift array
    :param N: number of iterations for the bootstrap
    :param zbinn: which bin is being used (str) "lowz" or "hiz"
    """

    import numpy as np
    from astropy.io import fits
    from numpy import random as ran
    from linetools.spectra import utils as ltsu
    from linetools.spectra.xspectrum1d import XSpectrum1D

    print("Number of spectra provided =", spec.nspec)

    # the bootstrap

    R = (spec.nspec)  # restraints on the random variable

    N_stack = []
    
    N_z_med = []

    for i in np.arange(N):  # this might take a while
        
        print(i)

        choice = np.asarray(ran.randint(0, R, R))  # the shuffle

        rand_spec = np.asarray([spec[index] for index in choice])
        
        rand_red = np.asarray([red[index] for index in choice])

        rand_collate = ltsu.collate(rand_spec)  # collate and stack

        N_stack.append(ltsu.smash_spectra(rand_collate))
        
        N_z_med.append(np.median(rand_red))
        
    # getting into a array for matrix math
        
    N_flux = np.array([entry.flux for entry in N_stack])
    
    #the actual stack
    
    collate  = ltsu.collate(spec)
    
    stack = ltsu.smash_spectra(collate)
    
    z_med = np.median(red)
    
    #to create the error array

    N_matrix = np.array([N_flux[i] - stack.flux for i in range(N)])

    tranpose = np.transpose(N_matrix)  # matrix math

    covariance = np.dot(tranpose, N_matrix)

    sigma = np.sqrt(np.diagonal(covariance) / (N - 1))  # the final error array
    
    #write out of the actual stack
    
    c_wave = fits.Column(name='WAVELENGTH', array=stack.wavelength, unit="angstroms", format="E") 
    c_flux = fits.Column(name='FLUX', array=stack.flux, unit="relative flux", format="E")
    c_noise = fits.Column(name='NOISE', array=sigma, unit="relative flux", format="E")
    t = fits.BinTableHDU.from_columns([c_wave, c_flux, c_noise])
    hdr = t.header
    hdr.set("REDSHIFT", z_med)
    t.writeto("/Users/jsmonzon/IGM-UCSC/fits/composites/"+zbin+"/stack.fits")

    #plot with error
    composite = XSpectrum1D(stack.wavelength, stack.flux, sig=sigma)
    composite.plot()
    
    #now to write out the randomized stacks with the same error array
    
    for j in range(N):
        
        print(j)
        
        c_wave = fits.Column(name='WAVELENGTH', array=N_stack[j].wavelength, unit="angstroms", format="E") 
        c_flux = fits.Column(name='FLUX', array=N_stack[j].flux, unit="relative flux", format="E")
        c_noise = fits.Column(name='NOISE', array=sigma, unit="relative flux", format="E")
        t = fits.BinTableHDU.from_columns([c_wave, c_flux, c_noise])
        hdr = t.header
        hdr.set("REDSHIFT", N_z_med[j])
        t.writeto("/Users/jsmonzon/IGM-UCSC/fits/bootstrap/"+zbin+"/composite_"+str(j)+".fits")
        
        