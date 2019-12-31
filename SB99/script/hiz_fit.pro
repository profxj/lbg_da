;This code fits the stellar continuum of restframe UV observations between 1230-1600 Angstroms (A) with stellar continuum models.
;I include an example spectrum to illustrate how to use the code.
;To use the code you will need ancillary programs including in the zip file:
;       (1) MPFIT.pro: a non-linear least-squares fitter (https://www.physics.wisc.edu/~craigm/idl/fitting.html)
;       (2) gconv.pro: a Gaussian convolution function
;       (3) calz_unred.pro: A program that reddens the stellar continuum models according to Calzetti et al. 2000. 
;To fit your own data you will need to:
;        (A) Put in your own data (line 29)
;        (B) Mask any intervening absorbers or reduction issues (line 76)
;        (C) Put in the spectral resolution (FWHM) of you data (line 108)



;This just tells IDL which directory you are in.
;Make sure you open IDL in the directory that you have all the files in. This includes MPFIT, gconv, and calz_unred. 
;Or you can make sure that all of these programs are in your IDL path
;path = '/Users/jsmonzon/IGM-UCSC/SB99/data/'


cspeed = 2.99E5 ;This is the speed of light. This constant will be used many places

;First we have to read in the data. I don't know how your data is organized.
;I would say that the best way to store the data is in an IDL ’structure’ where 
; you can keep things like wavelength, flux and error all together. 
; see here: http://slugidl.pbworks.com/w/page/29675376/IDL%20Structures
;But if your data is in .txt files we can change it.
;As part of this I am attaching a fits file with the spectrum of a galaxy
;from the MEGaSaURA sample (Rigby et al. 2018 AJ, 155, 104R)
s = mrdfits('/Users/jsmonzon/IGM-UCSC/fits/composites/hiz/stack.fits', 1)
;NOTE: the wavelength is the rest wavelength and the flux/err are in F_\lambda (ergs/s/cm^-2/A^-1)


;We are going to fit the observations with models.
;To keep things simple, we will normalize the flux of the data and models
;in the same wavelength. The normalized variable will be called norm.
;norm = median(s.flux[where(s.wavelength gt 1260 and s.wavelength lt 1304)])

;Now divide BOTH the error and the flux by this normalization constant.
;This keeps the relative SNR constant. 
;s.flux = s.flux/norm
;s.noise = s.noise/norm


;Now we have to read in the stellar continuum models. These models are
; a model of a single star at a single age and a single stellar 
; metallicity. There are ten ages that we include in the fits: 1, 2, 3, 4, 5, 8, 10, 20, 40 Myr.
; Further, there are 5 stellar metallicities included in the fits: 0.05, 0.2, 0.4, 1.0, 2.0 Z_\odot. 
models = mrdfits('models.fits', 1)
;The models are in a structure and contain a flux and a wavelength array. 
; You should get used to the models, and look for the differences
; in models as you increase the stellar age and metallicity (number of metals in the star).
; Plotting the models is done like:
; plot, models.wave, models.flux[*, 0], xr = [1500, 1575]
; IDL is 0 indexed, so the above selects the first model (1 Myr old and 0.05 Z_odot)
; The command xr = [1500, 1575] selects the x-range to be between 1500 and
; 1576 angstroms. This is near the VERY important C IV stellar wind feature.
; This line changes in depth and width as the stars get older and more metal rich.
; You can then over-plot the 2 Myr as:
; djs_oplot, models.wave, models.flux[*, 1], color = 'red'
;Try to plot all of the different models to get a feel of how the stellar
; wind line changes with age and metallicity.
; (hint: metallicity changes every 10 models so there are 50 models!!!
; djs_oplot, models.wave, models.flux[*, 10], color = 'red' 
; overplots a 1 Myr 0.2 Z_\odot model)


;When we do the fit, we only want to include points that provide
;information on the STARS, not the gas.
; This means that we want to mask out (not include) absorption lines from interstellar gas.
ISM_lines = [1031.92, 1036.33, 1048.21, 1083.99, 1117.97, 1122.52, 1128.00, 1144.93, 1152.81,$
             1190.41, 1193.28, 1197.18, 1199.39, 1199.54, 1206.50, 1215.67, 1238.82, 1260.42,$
             1302.17, 1304.37, 1335.71, 1393.76, 1402.77];, 1526.71, 1548.19, 1550.77]
;These lines will easily confuse and contaminate the stellar fitting.

;Higher redshift galaxies also have absorbtion from gas within the Inter galactic medium.
; You also want to remove that. This galaxy has an absorption feature at 1499A
;ISM_lines = [ISM_lines, 1499.] 

;At these spectral resolutions, the ISM featues are usually less than 1000 km/s in width.
;Therefore we are going to mask out 500 km/s on either side of the lines marked in ISM_lines.
dv = 500. ;velocity interval 

;Now we need to know how many lines we are masking
n_lines = n_elements(ISM_lines)

print, n_elements(ISM_lines)

;And then we need to create a mask array so we know which wavelengths should be included
mask = fltarr(n_elements(s.wavelength))+1.


;Now we will go through and mark the mask array with a 0 if we do not want to include it.
for maskn = 0, n_lines-1 do begin & $
  tempmask = where(s.wavelength gt -dv/cspeed*ISM_lines[maskn]+ISM_lines[maskn] and s.wavelength lt dv/cspeed*ISM_lines[maskn]+ISM_lines[maskn])   & $ ;This excludes +/-500km/s around each line
  if ISM_lines[maskn] eq 1550.77 then tempmask = where(s.wavelength gt -dv/cspeed*ISM_lines[maskn]+ISM_lines[maskn] and s.wavelength lt 50./cspeed*ISM_lines[maskn]+ISM_lines[maskn])   & $ ;CIV has a special stellar feature that we want redward of line center
  mask[tempmask] = 0 & $ ;This marks the bad wavelength regions in the mask array.
endfor

;Now we want to select the region of the spectrum that we want to fit.
;We do not want to go lower than 1225A because Ly-alpha absorption and emission 
; are strong in that region. Plus the Calzetti attenuation curve isn't defined down there.
good_reg = where(s.wavelength gt 1250 and s.wavelength lt 1450 and mask eq 1 and s.flux/s.noise gt 1, comp = bad_reg)
;bad_reg is the compliment, i.e. the stuff you don't want to include in the fit.

;Plot the wavelength regime to make sure that there is no unmasked narrow absorption lines 
;plot, s.wavelength, s.flux, xr = [1000, 1500], yr = [0,2], xtit = 'Rest Wavelength', ytit = 'Normalized flux'
;djs_oplot, s.wavelength[bad_reg], s.flux[bad_reg], color = 'gray'
;If you see narrow absorption lines in black and not in gray go back and add them to the ISM_lines call on line 76.


vdisp = 300. ; FWHM of spectral resolution of the observations in km/s


;This is where all of the work is done. The program entitled s99_continuum
;takes the models, the wavelength flux and error in the good regions (good_reg),
; and the velocity dispersion that we set and finds a linear combination of 
; the stellar models. It returns
; the best-fit parameters in a structure called cont_params and the best
; fit continuum model as cont_fit.
; 
cont_params = s99_continuum(models, s.wavelength, s.flux,$
  s.noise, good_reg, vdisp, yfit = cont_fit)
help, cont_params, /st
;Within cont_params has all the information that came out of the fit.
;Including: the fraction of light multiplied by each stellar continuum model (light_frac)
;           the best-fit attenuation parameter (E(B-V) or ebv)
;           the age of each stellar continuum model (model_age)
;           the errors on each of these parameters are also included as the name +_err

;The array cont_fit is the best-fit linear combination of the stellar continuum.

;Now we plot the observed flux and the fit
!p.font = 0 ;This makes the fonts look better
plot, s.wavelength, s.flux, yr = [0, 2], xr = [1050, 1400], /xs, xtitle = 'Wavelength', $
ytitle = 'Normalized Flux', charthick = 4
;djs_oplot, s.wavelength[bad_reg], s.flux[bad_reg], color = 'gray'
oplot, s.wavelength, s.noise, thick = 3 ;this is the error on the flux
oplot, s.wavelength, cont_fit, thick = 3; this is the stellar continuum fit

;We can calculate the average age and metallicity of the stellar continuum using a weighted average 
; of the light fractions.
exage = total(cont_params.light_frac*cont_params.model_age)/total(cont_params.light_frac)
exmet = total(cont_params.light_frac*cont_params.model_z)/total(cont_params.light_frac)

print,"Age",exage

print, "Metalicity", exmet

;hdr_1 = "COMMENT The continuum fit",$ "MKEY = "FIT"
;mwrfits, cont_fit, '/Users/jsmonzon/IGM-UCSC/fits/composites/hiz/model.fits'
;mwrfits, cont_params, 'high_z_param.fits'


;plot, s.wavelength[good_reg], s.flux[good_reg], xr = [1000, 1500], yr = [0,2], xtit = 'Rest Wavelength', ytit = 'Normalized flux'
;djs_oplot, s.wavelength[bad_reg], s.flux[bad_reg], color = 'gray'

;And we are finished with the initial pass through! There are a lot of ways we can improve this,
;But lets make sure you understand every step first..

end
