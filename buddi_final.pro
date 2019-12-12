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

pro BUDDI_final,ERGS=ergs

;read_input, input_file, setup



;*** Set up directory structure for all output files
root='/home/ejohnsto/Manga/7443-9102_ss/'
decomp='IFU_decomp/'
decomp_dir='decomposed_data/'
median_dir='median_image/'
binned_dir='binned_images/'
slices_dir='image_slices/'

directory=root

;fits_read,directory+file+'.fits',input,h
;
;;h = headfits(directory+file+'.fits')
;;h2 = headfits(directory+file+'.fits',exten=1)
;x=sxpar(h,'NAXIS1')
;y=sxpar(h,'NAXIS2')
;z=sxpar(h,'NAXIS3')




if keyword_set(ERGS) then begin
  fits_read,root+decomp+decomp_dir+'bulge_1D.fits',bulge_in,h_b
  fits_read,root+decomp+decomp_dir+'disk_1D.fits',disk_in,h_d
  z=sxpar(h_b,'NAXIS1')
output_ergs_b=fltarr(z)
output_ergs_d=fltarr(z)
  
  exptime=sxpar(h_b,'exptime')/sxpar(h_b,'nexp')
  wave_bands=[3543,4770,6231,7625,9134]
  zp_bands=[23.2680,24.2500,23.9270,23.6130,21.8650]
  
  ;zp_temp=linear_interpolate(wave,wave_bands,zp_bands)
  
  wavelength0=3.5598;sxpar(h_b,'CRVAL1')
  step=0.0001;sxpar(h_b,'CDELT1')
  wavelength=fltarr(sxpar(h_b,'NAXIS1'))
  for m=0,sxpar(h_b,'NAXIS1')-1,1 do wavelength[m]=wavelength0+(m*step)
  zp=linear_interpolate(10^wavelength,wave_bands,zp_bands)+15
  

  ;now use corrected zeropoints to convert flux to counts in the datacube
  datacube_b=fltarr(z)
  IFU_fnu_b=fltarr(z)
  IFU_maggies_b=fltarr(z)
  IFU_flam_b=fltarr(z)
  datacube_d=fltarr(z)
  IFU_fnu_d=fltarr(z)
  IFU_maggies_d=fltarr(z)
  IFU_flam_d=fltarr(z)
  
  IFU_maggies_b=cts2maggies(bulge_in,exptime,zp)
  IFU_maggies_d=cts2maggies(disk_in,exptime,zp)
  
  IFU_fnu_b=maggies2fnu(IFU_maggies_b)
  IFU_fnu_d=maggies2fnu(IFU_maggies_d)
  
  IFU_flam_b=fnu2flam(IFU_fnu_b,wavelength)/1e-17
  IFU_flam_d=fnu2flam(IFU_fnu_d,wavelength)/1e-17



  
  output_ergs_b=IFU_flam_b
  output_ergs_d=IFU_flam_d
endif 




fits_write,root+decomp+decomp_dir+'bulge_1D_ergs',output_ergs_b,h_b
fits_write,root+decomp+decomp_dir+'disk_1D_ergs',output_ergs_d,h_d
;mwrfits,0,root+decomp+decomp_dir+'bulge_1D_ergs.fits',h_b
;mwrfits,0,root+decomp+decomp_dir+'disk_1D_ergs.fits',h_d



fits_read,root+decomp+decomp_dir+'bulge.fits',bulge_in,h_b
fits_read,root+decomp+decomp_dir+'disk.fits',disk_in,h_d
x=sxpar(h_b,'NAXIS1')
y=sxpar(h_b,'NAXIS2')
z=sxpar(h_b,'NAXIS3')


datacube_b=fltarr(x,y,z)
IFU_fnu_b=fltarr(x,y,z)
IFU_maggies_b=fltarr(x,y,z)
IFU_flam_b=fltarr(x,y,z)
datacube_d=fltarr(x,y,z)
IFU_fnu_d=fltarr(x,y,z)
IFU_maggies_d=fltarr(x,y,z)
IFU_flam_d=fltarr(x,y,z)
for column=0,x-1,1 do begin
  for row=0,y-1,1 do begin
    IFU_maggies_b[column,row,*]=cts2maggies(bulge_in[column,row,*],exptime,zp)
    IFU_maggies_d[column,row,*]=cts2maggies(disk_in[column,row,*],exptime,zp)
    
    IFU_fnu_b[column,row,*]=maggies2fnu(IFU_maggies_b[column,row,*])
    IFU_fnu_d[column,row,*]=maggies2fnu(IFU_maggies_d[column,row,*])
    
    IFU_flam_b[column,row,*]=fnu2flam(IFU_fnu_b[column,row,*],wavelength)/1e-17
    IFU_flam_d[column,row,*]=fnu2flam(IFU_fnu_d[column,row,*],wavelength)/1e-17
    
  
  endfor
endfor
fits_write,root+decomp+decomp_dir+'bulge_ergs',IFU_flam_b,h_b
fits_write,root+decomp+decomp_dir+'disk_ergs',IFU_flam_d,h_d
;mwrfits,0,root+decomp+decomp_dir+'bulge_ergs.fits',h_b
;mwrfits,0,root+decomp+decomp_dir+'disk_ergs.fits',h_d


print,'#############################'
print,'## CRVAL3='+string(wavelength[0])+' ##'
print,'## CDELT3='+string(wavelength[1]-wavelength[0])+' ##'
print,'#############################'


stop
end