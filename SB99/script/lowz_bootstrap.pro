pro lowsedfitting,N

cspeed = 2.99E5 

models = mrdfits('models.fits', 1) ; the SB99 models

ISM_lines = [1031.92, 1036.33, 1048.21, 1083.99, 1117.97, 1122.52, 1128.00, 1144.93, 1152.81, 1190.41,$
            1193.28, 1197.18, 1199.39, 1199.54, 1206.50, 1215.67, 1238.82, 1260.42, 1302.17, 1304.37,$
            1335.71, 1393.76, 1402.77, 1499., 1526.71, 1548.19, 1550.77]
            
n_lines = n_elements(ISM_lines)

dv = 500. ;velocity interval for the masking of the lines

vdisp = 300. ; FWHM of spectral resolution of the observations in km/s

  for I = 0,N-1 do begin
  
    model_id = strn(I)
    
    print, model_id
  
    s = mrdfits('/Users/jsmonzon/IGM-UCSC/fits/bootstrap/lowz/composite_'+model_id+'.fits', 1) ;the read in
  
    ;norm = median(s.flux[where(s.wavelength gt 1260 and s.wavelength lt 1304)]) ; normalizing
    ;s.flux = s.flux/norm
    ;s.noise = s.noise/norm

    mask = fltarr(n_elements(s.wavelength))+1. ;for the ISM_lines
  
    for maskn = 0, n_lines-1 do begin & $
      tempmask = where(s.wavelength gt -dv/cspeed*ISM_lines[maskn]+ISM_lines[maskn] and s.wavelength lt dv/cspeed*ISM_lines[maskn]+ISM_lines[maskn])   & $ ;This excludes +/-500km/s around each line
      if ISM_lines[maskn] eq 1550.77 then tempmask = where(s.wavelength gt -dv/cspeed*ISM_lines[maskn]+ISM_lines[maskn] and s.wavelength lt 50./cspeed*ISM_lines[maskn]+ISM_lines[maskn])   & $ ; CIV has a special stellar feature that we want redward of line center
      mask[tempmask] = 0 & $ ;This marks the bad wavelength regions in the mask array.
    endfor
  
    good_reg = where(s.wavelength gt 1225 and s.wavelength lt 1500 and mask eq 1 and s.flux/s.noise gt 1, comp = bad_reg)
  
    cont_params = s99_continuum(models, s.wavelength, s.flux,$ ;the actual fitting!
    s.noise, good_reg, vdisp, yfit = cont_fit)
    help, cont_params, /st
  
    ;Now we plot the observed flux and the fit
    ;!p.font = 0 ;This makes the fonts look better
    ;plot, s.wavelength, s.flux, yr = [0, 2], xr = [1000, 1500], /xs, xtitle = 'Wavelength', $
    ;ytitle = 'Normalized Flux', charthick = 4
    ;oplot, s.wavelength[bad_reg], s.flux[bad_reg], color = 'gray'
    ;oplot, s.wavelength, s.noise, thick = 2 ;this is the error on the flux
    ;oplot, s.wavelength, cont_fit, thick = 3; this is the stellar continuum fit
  
  
    ;the write out
    mwrfits, cont_fit, '/Users/jsmonzon/IGM-UCSC/fits/bootstrap/lowz/model_'+model_id+'.fits'

  endfor

end

