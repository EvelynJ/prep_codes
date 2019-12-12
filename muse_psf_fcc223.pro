; code to create PSF profiles from standard stars for MUSE galaxies


pro MUSE_psf_FCC223;, input_file
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
root='/raid/ejohnston/FCC223/BUDDI/'
  side=20.0
  
  
  ;read in a text file containing the list of styandard star datacubes
  ;readcol,root+'std_input.txt',format='a',files,/SILENT,comment='#'
  
;  fits_read,root+'DATACUBE_FINAL_ZAPPED2.fits',input,h
  fits_read,'/raid/ejohnston/FCC223/products_std/DATACUBE_FINAL_ZAPPED.fits',input,h
;  fits_read,'/raid/ejohnston/FCC223/BUDDI/DATACUBE_FINAL_ZAPPED2_FLUX.fits',input,h
  psf_cube=fltarr(41,41s,fix(sxpar(h,'NAXIS3')/100))
  psf_cube2=fltarr(41,41,sxpar(h,'NAXIS3'))
  psf_fwhm_wave=fltarr(sxpar(h,'NAXIS3'))
  
;  x_in=[160, 245,  96, 288, 41,  63, 49]   ;using std datacube
;  y_in=[260, 275, 146,  52, 48, 130, 96]
  x_in=[209, 193, 225,  19, 137, 163, 262, 281,  33]   ;using std datacube
  y_in=[285, 230, 251, 252, 237,  59, 116, 208, 296]
  ;for each file, read in
;  for n=0,n_elements(files)-1,1 do begin
;    fits_read,files[n],input,h
;    all_stars[(2*side+1)*n:(2*side+1)*(n+1)-1,*,*]=input[x_cen-side:x_cen+side,y_cen-side:y_cen+side,*]
;    x_in=[x_in,(n*side*2)+side]
;    y_in=[y_in,side]
;  endfors
  
  ;combine the datcubes to create PSF profile cube
  openw,01,root+'psf_FWHM.txt'
;  for m=500,600,1 do begin
  j=0
  for m=51,sxpar(h,'NAXIS3')-51,100 do begin
    print,m
;    if m mod 10.0 eq 0 then print,string(m)+'/'+string(sxpar(h,'NAXIS3'))
;    if m lt 3 then low=0 else low=2
;    if m gt sxpar(h,'NAXIS3')-6 then high=0 else high=2
    PSF_EXTRACT,x_in-1,y_in-1,0,0,median(input[*,*,m-50:m+50],dimension=3),2*side+1,psf_new,Psf_fwhm,Background,/SKY_MEDIAN
;    Psf_fwhm*=1.102    ;there seems to be some offset between the FWHM of the input stars and the gaussian model
;    psf_new=PSF_GAUSSIAN( Npixel=101, FWHM=Psf_fwhm, /NORMAL )
    psf_cube[*,*,j]=psf_new
    psf_fwhm_wave[j]=Psf_fwhm
    printf,01,Psf_fwhm
    j+=1
  endfor
  psf_cube_final=psf_cube;[2:*,2:*,*]
  close,01
  
  ;interpolate between the binned PSF profiles
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


  fits_write,root+'psf.fits',psf_cube_final,h
  stop
  
end