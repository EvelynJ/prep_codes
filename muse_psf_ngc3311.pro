; code to create PSF profiles from standard stars for MUSE galaxies


pro MUSE_PSF_NGC3311, input_file

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
side=30.0


directory=root
  
;readcol,directory+'star_centres.txt', F='F,F',x_in,y_in,comment='#',/silent
fits_read,directory+file+'.fits',input,h
psf_cube=fltarr(2*side+1,2*side+1,sxpar(h,'NAXIS3'))
psf_fwhm_wave=fltarr(sxpar(h,'NAXIS3'))

x_in=[218.000,447.000,457.000,507.000,484.000,262.000]
y_in=[1204.00,928.000,835.000,570.000,571.000,313.000]

;combine the datcubes to create PSF profile cube
for m=0,sxpar(h,'NAXIS3')-1,1 do begin
  in=input[*,*,m]
  PSF_EXTRACT,x_in,y_in,0,0,in,2*side+1,psf_new,Psf_fwhm,Background,/SKY_MEDIAN
  psf_cube[*,*,m]=psf_new
  psf_fwhm_wave[m]=Psf_fwhm
endfor
psf_cube_final=psf_cube[2:*,2:*,*]

fits_write,root+'psf.fits',psf_cube_final,h



end