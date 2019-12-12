; code to create PSF profiles from standard stars for MUSE galaxies


pro MUSE_PSF_Es, input_file

read_input, input_file, setup
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
side=20.0


directory=root

readcol,directory+'star_centers.txt', F='F,F',x_in,y_in,comment='#',/silent
fits_read,directory+file+'.fits',input,h
psf_cube=fltarr(2*side+1,2*side+1,sxpar(h,'NAXIS3'))
psf_fwhm_wave=fltarr(sxpar(h,'NAXIS3'))

;if galaxy_ref eq '2MIG1814' then begin
;  fits_read,directory+'psf.fits',input
;  x_in=20
;  y_in=20
;endif

;x_in=[218.000,447.000,457.000,507.000,484.000,262.000]
;y_in=[1204.00,928.000,835.000,570.000,571.000,313.000]

;;combine the datcubes to create PSF profile cube
;for m=0,sxpar(h,'NAXIS3')-1,1 do begin
;;for m=100,100,1 do begin
;  in=input[*,*,m]
;  PSF_EXTRACT,x_in,y_in,0,0,in,2*side+1,psf_new,Psf_fwhm,Background,/SKY_MEDIAN
;  psf_cube[*,*,m]=psf_new
;  psf_fwhm_wave[m]=Psf_fwhm
;endfor
;
;psf_cube_final=psf_cube[2:*,2:*,*]
;
;nan=where(psf_cube_final[side,side,*] eq 0 or finite(psf_cube_final[side,side,*]) eq 0)
;no_nan=where(psf_cube_final[side,side,*] ne 0 or finite(psf_cube_final[side,side,*]) eq 1)
;
;;fix any slices where the values are NAN, likely first or final slices
;if nan[0] ne -1 then begin
;  for n=0,n_elements(nan)-1,1 do psf_cube_final[*,*,nan[n]]=psf_cube_final[*,*,no_nan[0]]
;endif
;
;fits_write,root+'psf.fits',psf_cube_final,h


openw,01,root+'psf_FWHM.txt'
psf_cube=fltarr(2*side+1,2*side+1,fix(sxpar(h,'NAXIS3')/100))
psf_cube2=fltarr(2*side+1,2*side+1,sxpar(h,'NAXIS3'))
;  for m=500,600,1 do begin
j=0
for m=51,sxpar(h,'NAXIS3')-51,100 do begin
  print,m
  ;    if m mod 10.0 eq 0 then print,string(m)+'/'+string(sxpar(h,'NAXIS3'))
  ;    if m lt 3 then low=0 else low=2
  ;    if m gt sxpar(h,'NAXIS3')-6 then high=0 else high=2
  PSF_EXTRACT,x_in,y_in,0,0,median(input[*,*,m-50:m+50],dimension=3),2*side+1,psf_new,Psf_fwhm,Background,/SKY_MEDIAN
  psf_new=PSF_GAUSSIAN( Npixel=2*side+1, FWHM=Psf_fwhm, /NORMAL )
  
  psf_cube[*,*,j]=psf_new
  psf_fwhm_wave[j]=Psf_fwhm
  printf,01,Psf_fwhm
  j+=1
endfor
psf_cube_final=psf_cube;[2:*,2:*,*]
close,01


for m=0,sxpar(h,'NAXIS3')-1,1 do begin
  if m lt 51 then psf_cube2[*,*,m]=psf_cube[*,*,0] $
  else if m ge 3500 then psf_cube2[*,*,m]=psf_cube[*,*,-1] $
  else begin
    n=fix(m/100)
    frac_1=(m-(n*100.))/100.
    frac_2=1.-frac_1
    min_val=min(frac_2*psf_cube[*,*,n] + frac_1*psf_cube[*,*,n+1])
    max_val=max(frac_2*psf_cube[*,*,n] + frac_1*psf_cube[*,*,n+1])
    im=(frac_2*psf_cube[*,*,n] + frac_1*psf_cube[*,*,n+1])
    ;if min_val lt 0 then im=im+min_val
    psf_cube2[*,*,m]=im;/total(im)
  endelse

endfor
psf_cube_final=psf_cube2;[2:*,2:*,*]
fits_write,root+'psf_smooth.fits',psf_cube_final,h

end