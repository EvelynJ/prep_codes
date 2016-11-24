; quick code to prepare datacubes for BUDDI.
; The code will:
;   - convert flux back into counts in a really dirty way (/counts)
;   - log-rebin the datacube in the wavelength direction if 
;     it's linearly binned (/log_rebin)
; 

;===============================================================

FUNCTION cts2maggies, cts, expt, zp
;zp: zeropoint is a positive number
   return, cts / expt * 10^(-0.4*zp)
END
FUNCTION maggies2fnu, maggies
   return, 3631e-23*maggies     ;[erg s-1 Hz-1 cm-2]
END
FUNCTION fnu2flam, fnu, lam     ;[erg s-1 Hz-1 cm-2], [angstroem]
   return, 299792458.*1e10/lam^2.*fnu ;[erg s-1 cm-2 A-1]
END


FUNCTION flam2fnu, flam, lam    ;[erg s-1 cm-2 A-1], [angstroem]
   return, flam/299792458./1e10*lam^2. ;[erg s-1 Hz-1 cm-2]
END
FUNCTION fnu2maggies, fnu
   return, 3631e23*fnu
END
FUNCTION maggies2cts, maggies, expt, zp
;zp: zeropoint is a positive number
   return, maggies * expt / 10^(-0.4*zp)
END


FUNCTION mags2cts, mags, expt, zp
;zp: zeropoint is a positive number
   return, maggies2cts(mags2maggies(mags), expt, zp)
END
FUNCTION mags2maggies, mags
   return, 10^(-0.4*mags)
END


FUNCTION cts2mags, cts, expt, zp
;zp: zeropoint is a positive number
   return, maggies2mags(cts2maggies(cts, expt, zp))
END
FUNCTION maggies2mags, maggies
   return, -2.5*alog10(maggies)
END

;===============================================================

pro BUDDI_prep,input_file,COUNTS=counts,LOG_REBIN=log_rebin,MANGA=manga

read_input, input_file, setup



;*** Set up directory structure for all output files
root=setup.root
galaxy_ref=setup.galaxy_ref
file=setup.file
kinematics=setup.kinematics
decomp=setup.decomp
median_dir=setup.median_dir
binned_dir=setup.binned_dir
slices_dir=setup.slices_dir
psf_file=setup.psf_file
stellib_dir=setup.stellib_dir

directory=root

fits_read,directory+file+'.fits',input,h
h = headfits(directory+file+'.fits')
h2 = headfits(directory+file+'.fits',exten=1)
x=sxpar(h2,'NAXIS1')
y=sxpar(h2,'NAXIS2')
z=sxpar(h2,'NAXIS3')


;write out wavelength array
if keyword_set(MANGA) then begin
    wavelength0_lin=sxpar(h2,'CRVAL3')
    wavelength0=alog10(wavelength0_lin)
    step_lin=sxpar(h2,'CD3_3')
    step=alog10(wavelength0_lin+step_lin)-wavelength0
    wavelength=fltarr(sxpar(h2,'NAXIS3'))
    for m=0,sxpar(h2,'NAXIS3')-1,1 do wavelength[m]=wavelength0+(m*step)
    wave=10^wavelength
  
    sxaddpar,h2,'CRVAL3',wavelength[0]
    sxaddpar,h2,'CDELT3',wavelength[1]-wavelength[0]
    sxaddpar,h2,'CD3_3',wavelength[1]-wavelength[0]
    
endif else begin
    wavelength0=sxpar(h2,'CRVAL3')
    step=sxpar(h2,'CD3_3')
    wavelength=fltarr(sxpar(h2,'NAXIS3'))
    for m=0,sxpar(h2,'NAXIS3')-1,1 do wavelength[m]=wavelength0+(m*step)
  
endelse







input_counts=fltarr(x,y,z)
if keyword_set(COUNTS) then begin
  exptime=sxpar(h,'exptime')/sxpar(h,'nexp')
  wave_bands=[3543,4770,6231,7625,9134]
  zp_bands=[23.2680,24.2500,23.9270,23.6130,21.8650]
  
  zp_temp=linear_interpolate(wave,wave_bands,zp_bands)
  
;  zp=fltarr(z)
;  
;  result=poly_fit(wave_bands,zp_bands,2)
;  for j=0,z-1,1 do zp[j]=result[0]+result[1]*(wave[j])+result[2]*(wave[j])^2
  zp=linear_interpolate(wavelength,wave_bands,zp_bands)
;    
  
  ;now use corrected zeropoints to convert flux to counts in the datacube
  datacube=fltarr(x,y,z)
  IFU_fnu=fltarr(x,y,z)
  IFU_maggies=fltarr(x,y,z)
  IFU_mags=fltarr(x,y,z)
  for column=0,x-1,1 do begin
    for row=0,y-1,1 do begin
      spec_xy=input[column,row,*]*1e-17
      IFU_fnu[column,row,*]=flam2fnu(spec_xy,(wavelength))
      IFU_maggies[column,row,*]=fnu2maggies(IFU_fnu[column,row,*])
      datacube[column,row,*]=maggies2cts(IFU_maggies[column,row,*],exptime,zp)
        

      
    endfor
  endfor
  
  input_counts=datacube
endif else input_counts=input


output=fltarr(x,y,z)
if keyword_set(LOG_REBIN) then begin
  lamRange2=[wavelength[0],wavelength[-1]]
  for column=0,x-1,1 do begin
    for row=0,y-1,1 do begin
      temp_spec=fltarr(z)
      temp_spec[*]=input_counts[column,row,*]
      log10_rebin, lamRange2, temp_spec, new_spec, logLam2, VELSCALE=velScale
      output[column,row,*]=new_spec
    endfor
  endfor
  sxaddpar,h2,'CRVAL3',logLam2[0]
  sxaddpar,h2,'CDELT3',logLam2[1]-logLam2[0]
  sxaddpar,h2,'CD3_3',logLam2[1]-logLam2[0]
endif else output=input_counts

fits_write,directory+file+'_counts.fits',output,h2



;mkhdr, hdr,var,/IMAGE
;sxaddhist,"Variance image",hdr,/COMMENT
mwrfits,0,directory+file+'_counts.fits',h
;mkhdr, hdr,badpix,/image
;sxaddhist,"Bad pixel mask",hdr,/COMMENT
;mwrfits,badpix,'/home/ejohnsto/NGC3311/data/NGC3311.fits',hdr

end