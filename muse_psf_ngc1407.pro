; code to create PSF profiles from standard stars for MUSE galaxies


pro MUSE_psf_NGC1407, input_file
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
  
  
  ;read in a text file containing the list of styandard star datacubes
  ;readcol,root+'std_input.txt',format='a',files,/SILENT,comment='#'
  
  fits_read,'/home/ejohnsto/NGC1407_raw/products_D/DATACUBE_STD_0001.fits',input,h
  psf_cube=fltarr(2*side+1,2*side+1,sxpar(h,'NAXIS3'))
  psf_fwhm_wave=fltarr(sxpar(h,'NAXIS3'))
  
  x_cen=195+1
  y_cen=156+1
  x_in=[194,260];,63,50,42]
  y_in=[161,88];,133,99,51]
  ;for each file, read in
;  for n=0,n_elements(files)-1,1 do begin
;    fits_read,files[n],input,h
;    all_stars[(2*side+1)*n:(2*side+1)*(n+1)-1,*,*]=input[x_cen-side:x_cen+side,y_cen-side:y_cen+side,*]
;    x_in=[x_in,(n*side*2)+side]
;    y_in=[y_in,side]
;  endfor
  
  ;combine the datcubes to create PSF profile cube
  openw,01,root+'psf_FWHM.txt'
  for m=0,sxpar(h,'NAXIS3')-1,1 do begin
    PSF_EXTRACT,x_in,y_in,0,0,input[*,*,m],2*side+1,psf_new,Psf_fwhm,Background,/SKY_MEDIAN
    psf_cube[*,*,m]=psf_new
    psf_fwhm_wave[m]=Psf_fwhm
    printf,01,Psf_fwhm
  endfor
  psf_cube_final=psf_cube[2:*,2:*,*]
  close,01
  

  ;fits_write,root+'psf.fits',psf_cube_final,h
  stop
  
end