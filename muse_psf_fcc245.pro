; code to create PSF profiles from standard stars for MUSE galaxies


pro MUSE_psf_FCC245;, input_file
;  read_input, input_file, setup
;  root=setup.root
;  galaxy_ref=setup.galaxy_ref
;  file=setup.file
;  kinematics=setup.kinematics
;  decomp=setup.decomp
;  median_dir=setup.median_dir
;  binned_dir=setup.binned_dir
;  slices_dir=setup.slices_dir
;  psf_file=setup.psf_file
;  stellib_dir=setup.stellib_dir
root='/raid/ejohnston/FCC245/BUDDI/'
  side=25.0
  
  
  ;read in a text file containing the list of styandard star datacubes
  ;readcol,root+'std_input.txt',format='a',files,/SILENT,comment='#'
  
;  fits_read,root+'DATACUBE_FINAL_ZAPPED2.fits',input,h
  fits_read,'/raid/ejohnston/FCC245/products_std/DATACUBE_FINAL_ZAPPED.fits',input,h
  psf_cube=fltarr(2*side+1,2*side+1,sxpar(h,'NAXIS3'))
  psf_fwhm_wave=fltarr(sxpar(h,'NAXIS3'))
  
  x_in=[253, 138,  43, 230]   ;using std datacube
  y_in=[192, 208, 143, 248]
  ;for each file, read in
;  for n=0,n_elements(files)-1,1 do begin
;    fits_read,files[n],input,h
;    all_stars[(2*side+1)*n:(2*side+1)*(n+1)-1,*,*]=input[x_cen-side:x_cen+side,y_cen-side:y_cen+side,*]
;    x_in=[x_in,(n*side*2)+side]
;    y_in=[y_in,side]
;  endfor
  
  ;combine the datcubes to create PSF profile cube
  openw,01,root+'psf_FWHM.txt'
;  for m=500,600,1 do begin
  for m=0,sxpar(h,'NAXIS3')-1,1 do begin
    if m mod 10.0 eq 0 then print,string(m)+'/'+string(sxpar(h,'NAXIS3'))
    if m lt 3 then low=0 else low=2
    if m gt sxpar(h,'NAXIS3')-6 then high=0 else high=2
    PSF_EXTRACT,x_in,y_in,0,0,median(input[*,*,m-low:m+high],dimension=3),2*side+1,psf_new,Psf_fwhm,Background,/SKY_MEDIAN
    psf_cube[*,*,m]=psf_new
    psf_fwhm_wave[m]=Psf_fwhm
    printf,01,Psf_fwhm
  endfor
  psf_cube_final=psf_cube[2:*,2:*,*]
  close,01
  

  fits_write,root+'psf.fits',psf_cube_final,h
  stop
  
end