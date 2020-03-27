pro sedfitter,zbin

;This routine fits SEDs to the 5000 bootstrap generated composites from the "zbin" redshift interval
;as well as the "actual" composite. It returns the models individually as fits files. The routine reads in two fits
;files that are included in the 'spec_data' subdirectory: 'zbin_stacks.fits' and 'zbin_composite.fits'

;select from the following redshift bins
;zbin = "lowz"  2.25<z<2.5
;zbin = "hiz"   2.5<z<2.75

;John extended the dv, masked some emission lines, and made the fitting regime was equal between the bootstrap and composite fits. 

;---------------------------------------
;intializing some things
;---------------------------------------

models = mrdfits('/Users/johnchisholm/Downloads/for_john 2/S99_models.fits', 1) ; the SB99 models

cspeed = 2.99E5 
dv = 600. ;velocity interval for the masking of the lines
vdisp = 300.; FWHM of spectral resolution of the observations in km/s

;---------------------------------------
;reading in the bootstrap data
;---------------------------------------

spec_table = mrdfits('/Users/johnchisholm/Downloads/for_john 2/spec_data/'+zbin+'_stacks.fits',1,data) ; the wavelength and flux_err arrays (372)
fluxes = mrdfits('/Users/johnchisholm/Downloads/for_john 2/spec_data/'+zbin+'_stacks.fits',2,data) ; the flux iteration matrix (5000x372) 

wave = spec_table.wavelength ;the same for every bootstrap iteration
flux_err = spec_table.flux_err ;the same for every bootstrap iteration

;---------------------------------------
;masking out the ISM lines
;---------------------------------------

ISM_lines = [1031.92, 1036.33, 1048.21, 1083.99, 1117.97, 1122.52, 1128.00, 1144.93, 1152.81, 1190.41,$
            1193.28, 1197.18, 1199.39, 1199.54, 1206.50, 1215.67, 1238.82, 1260.42, 1262., 1264., 1266., 1299., 1302.17, 1304.37, 1307., 1310., 1333.,$
            1335.71, 1393.76, 1402.77, 1499., 1526.71, 1548.19, 1550.77]

n_lines = n_elements(ISM_lines)
mask = fltarr(n_elements(wave))+1. ;for the ISM_lines

for maskn = 0, n_lines-1 do begin & $
   tempmask = where(wave gt -dv/cspeed*ISM_lines[maskn]+ISM_lines[maskn] and wave lt dv/cspeed*ISM_lines[maskn]+ISM_lines[maskn]) & $
   if ISM_lines[maskn] eq 1550.77 then tempmask = where(wave gt -dv/cspeed*ISM_lines[maskn]+ISM_lines[maskn] and $
                          wave lt 50./cspeed*ISM_lines[maskn]+ISM_lines[maskn]) & $
   mask[tempmask] = 0 & $ ;This marks the bad wavelength regions in the mask array.
endfor

;---------------------------------------
;the bootstrap fitting loop!
;---------------------------------------

;i am hoping that fluxes[ii] takes the ii'th ROW in the 5000x327 table
for ii = 0,4999 do begin  
    flux_t = fluxes[*, ii]
    good_reg = where(wave gt 1225 and wave lt 1440 and mask eq 1 and fluxes[*,ii]/flux_err gt 1, comp = bad_reg)
    cont_params = s99_continuum(models, wave, reform(fluxes[*,ii]), flux_err, good_reg, vdisp, yfit = cont_fit)
    mwrfits, cont_fit, '/Users/johnchisholm/Downloads/for_john 2/fitting_results/'+zbin+'_bootstrap/sed_'+strn(ii)+'.fits'
endfor

;---------------------------------------
;fitting the actual composite
;---------------------------------------

s = mrdfits('/Users/johnchisholm/Downloads/for_john 2/spec_data/'+zbin+'_composite.fits', 1)

good_reg = where(s.wavelength gt 1225 and s.wavelength lt 1440 and mask eq 1 and s.flux/s.flux_err gt 1, comp = bad_reg)
cont_params = s99_continuum(models, s.wavelength, s.flux, s.flux_err, good_reg, vdisp, yfit = cont_fit)

;Now we plot the observed flux and the fit
!p.font = 0
plot, s.wavelength, s.flux, yr = [0, 2], xr = [1050, 1500], /xs, xtitle = 'Wavelength', $
ytitle = 'Normalized Flux', charthick = 4
oplot, s.wavelength, s.flux_err, thick = 3 ;this is the error on the flux
oplot, s.wavelength, cont_fit, thick = 3 ;this is the stellar continuum fit

;We can calculate the average age and metallicity of the stellar continuum using a weighted average of the light fractions.
exage = total(cont_params.light_frac*cont_params.model_age)/total(cont_params.light_frac)
exmet = total(cont_params.light_frac*cont_params.model_z)/total(cont_params.light_frac)

;---------------------------------------
;writing out
;---------------------------------------

;it would be nice if you could print our the chi^2 values from cont_params or cont_fit
;print, "Chi^2 value in the 1225-1450A fitting region for the "+zbin+"redshift bin", chisq
chisq = total((s[good_reg].flux-cont_fit[good_reg])^2/s[good_reg].flux_err^2)/double(n_elements(good_reg))
print,"Age of the SED model for the "+zbin+" redshift bin", exage
print, "Metalicity of the SED for the "+zbin+" redshift bin", exmet
print, "Chi^2 value in the 1225-1450A fitting region for the "+zbin+"redshift bin", chisq

mwrfits, cont_fit, '/Users/johnchisholm/Downloads/for_john 2/fitting_results/'+zbin+'_composites/sed.fits'

end
