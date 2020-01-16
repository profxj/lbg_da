function gconv, x, sigma, edge_wrap=edge_wrap, fwhm = fwhm
  ;; Convolve the x-vector with a Gaussian kernel - the kernel size
  ;; is set to 4 times the sigma of the Gaussian.

  ; x should be a float!!!

  if (n_elements(fwhm) gt 0) then $
    sigma = fwhm/2.3548

  ;; special case for no smoothing
  if sigma eq 0 then return, x

  binfactor=1
  ksize=round(4.0*sigma+1.0)*2
  xx = findgen(ksize)-ksize/2

  kernel=exp(-xx^2/(2*sigma^2))
  kernel=kernel/total(kernel)

  sm = convol(x, kernel, edge_wrap=edge_wrap)

  return, sm
end



function s99_continuum, model, restwl, flux, err, good_reg, vdisp, $
         yfit = yfit

nmodels = n_elements(model.age) ;This is the total number of continuum models

; Define the output structure. This includes light fractions, model ages, model metallicites,
; and the errors on each.
coefs = {light_frac: fltarr(nmodels), light_frac_err: fltarr(nmodels), $
         model_age: model.age, model_z: model.z, ebv: 0., ebv_err:0.}


;-----------------------------------------------------------------------------
;Now we need to establish the initial and boundary conditions for the fitted parameters
; MPFIT requires that you make a structure with a starting value, a fixed flag, a limited flag,
; a tied flag, and the actual limits.
; There are currently 51 free parameters: 50 coefficients for continuum models and 1 for attenuation.
; See https://www.physics.wisc.edu/~craigm/idl/fitting.html for details on MPFIT


;The way the structure is set up has indicies 1 through nmodels-1 for stellar continuum. The stellar attenuation is index 0.
parinfo = replicate({value:0.D, fixed:0, limited:[0,0], tied:'', $
                    limits:[0.D,0.D]}, nmodels+1) ;nmodels = number of stellar continuum, 3 for LyB 
 
;I changed the nmodels+1 to nmodels+2

parinfo[0].limited = [1,1] ;This determines if the zeroth parameter (E(B-V)) has limits
parinfo[0].limits = [0.0,5.0]   ; This gives the range of possible E(B-V) values
parinfo[0].value = 0.3     ;This is the initial E(B-V) value 

parinfo[1:*].limited = [1,0]    ; This limits the light fraction to not be negative, but there is no upper limit
parinfo[1:*].limits = [0.0,0.0] ; This puts the lower-limit on the light fraction at 0, but no upper limit
parinfo[1:*].value = 0.1     ;This is the initial guess for each light fraction.   

;I added this parameter to describe the effects of the IGM

;parinfo[nmodels+1].limited = [1,1] ;This determines if the last parameter (T_IGM) has limits
;parinfo[nmodels+1].limits = [0.0,1]   ; This gives the range of possible T_IGM values
;parinfo[nmodels+1].value = 0.5    ;This is the initial T_IGM value       
                    
                  

;-----------------------------------------------------------------------------
; Convolve models to velocity dispersion of data and interpolate to
; match data
s99_pix = 146.107 ; size of 1 model pixel in km/s 
s99_vdisp = 92.0 ; approximate velocity dispersion of S99 models
npix = n_elements(restwl) ;number of wavelength pixels

;Deconvolve template instrumental resolution, 
if vdisp lt s99_vdisp then vdisp_add = 0 $
else vdisp_add = sqrt(vdisp^2 - s99_vdisp^2)  
sigma_pix = vdisp_add / s99_pix

custom_lib = dblarr(npix, nmodels)  ;creates the array for the new models
for ii = 0, nmodels - 1 do custom_lib[*,ii] = gconv(interpol(model.flux[*,ii], model.wave, restwl), sigma_pix, fwhm=sigma_pix)& $ 
;This puts the models at the same resolution and on the same wavelength grid as the observations.

;-------------------------------------------------------------------------------
;Now we are actually going to run MPFIT
fitcoefs = mpfitfun('s99_mcombine', restwl[good_reg], flux[good_reg], err[good_reg], $
                     parinfo = parinfo, functargs = {mlib: custom_lib[good_reg, *], vdisp: vdisp}, $
                     perror=perror, niter=niter, status=status, /quiet, $
                     maxiter = 500, errmsg = errmsg)
print, 'CONTINUUM FIT ITERATIONS: ', strtrim(niter, 2)
print, 'CONTINUUM_FIT EXIT STATUS: ', strtrim(status, 2)
;We want to make sure that the exit status is always 1


;MPFIT puts the best fit models into fitcoefs and the errors on those parameters on perror.

;We then put this parameters into the output structure
coefs.light_frac = fitcoefs[1:nmodels-1]
coefs.light_frac_err = perror[1:nmodels-1]
coefs.ebv = fitcoefs[0]
coefs.ebv_err = perror[0]
; fit to full spectrum including masked pixels
yfit = s99_mcombine(restwl, fitcoefs, mlib=custom_lib)

;Now we can return the fitted coefficients
return, coefs

end
