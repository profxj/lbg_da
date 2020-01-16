pro reddy_unred, wave, flux, ebv, funred, R_V = R_V
;This program is uses the attenuation curve from Reddy et al. 2016 ApJ, 828, 107R

  if N_elements(R_V) EQ 0 then R_V = 4.05
  w1 = where((wave GE 6300) AND (wave LE 22000), c1)
  w2 = where((wave GE  912) AND (wave LT  6300), c2)
  x  = wave*0.0001                     ;Wavelength in microns

  ;IF (c1 + c2) NE N_elements(wave) THEN message,/INF, $
   ; 'Warning - some elements of wavelength vector outside valid domain'
    
  klam = 2.191 +.974/x

  funred = flux*10.0^(0.4*klam*ebv)
  if N_params() EQ 3 then flux = funred

end

pro calz_unred, wave, flux, ebv, funred, R_V = R_V
  On_error, 2
;This reddens the stellar continuum according to a Calzetti attenuation curve from Calzetti et al. 2000 ApJ, 533, 682C
  if N_params() LT 3 then begin
    print,'Syntax: CALZ_UNRED, wave, flux, ebv, [ funred, R_V=]'
    return
  endif

  if N_elements(R_V) EQ 0 then R_V = 4.05
  w1 = where((wave GE 6300) AND (wave LE 22000), c1)
  w2 = where((wave GE  912) AND (wave LT  6300), c2)
  x  = 10000.0/wave                      ;Wavelength in inverse microns

  IF (c1 + c2) NE N_elements(wave) THEN message,/INF, $
    'Warning - some elements of wavelength vector outside valid domain'

  klam = 0.0*flux

  IF c1 GT 0 THEN $
    klam[w1] = 2.659*(-1.857 + 1.040*x[w1]) + R_V

  IF c2 GT 0 THEN $
    klam[w2] = 2.659*(poly(x[w2], [-2.156, 1.509d0, -0.198d0, 0.011d0])) + R_V

  funred = flux*10.0^(0.4*klam*ebv)
  if N_params() EQ 3 then flux = funred

end

function s99_mcombine, x, a, mlib=mlib

; Create a linear combination of the templates.
; mlib are the templates.
; a is the array with the coefficients to create the linear combination.
; The first element of a is the attenuation (E(B-V))
; The rest of the array is the linear combination coefficients.
y = double(mlib) # a[1:*] ;This is a linear combination of the templates 

;Pick which reddening law you use. Calzetti is the most "standard"
; Calzetti reddening
calz_unred, x, y, -1.*a[0], ynew

; Reddy reddening
;reddy_unred, x, y, -1.*a[0], ynew ;This reddens the spectra with the given attenuation. 

return, ynew ;This returns the reddened linear combination. 

end
