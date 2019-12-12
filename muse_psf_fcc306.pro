; code to create PSF profiles from standard stars for MUSE galaxies


pro MUSE_psf_FCC306


root='/data2/ejohnston/Fornax_new_sample/FCC306/'
  side=25.0
  
  
  ;read in a text file containing the list of styandard star datacubes
  ;readcol,root+'std_input.txt',format='a',files,/SILENT,comment='#'
  
;  fits_read,'/raid/ejohnston/FCC215/products_std/DATACUBE_FINAL_ZAPPED.fits',input,h
  fits_read,root+'products/DATACUBE_FINAL_ZAPPED.fits',input,h
  psf_cube=fltarr(2*side+1,2*side+1,fix(sxpar(h,'NAXIS3')/100))
  psf_cube2=fltarr(2*side+1,2*side+1,sxpar(h,'NAXIS3'))
  psf_fwhm_wave=fltarr(sxpar(h,'NAXIS3'))
  
  x_in=[302,234]
  y_in=[219,115]


  openw,01,root+'BUDDI/psf_FWHM.txt'
  psf_cube=fltarr(101,101,fix(sxpar(h,'NAXIS3')/100))
  psf_cube2=fltarr(101,101,sxpar(h,'NAXIS3'))
  ;  for m=500,600,1 do begin
  j=0
  for m=51,sxpar(h,'NAXIS3')-51,100 do begin
    print,m
    ;    if m mod 10.0 eq 0 then print,string(m)+'/'+string(sxpar(h,'NAXIS3'))
    ;    if m lt 3 then low=0 else low=2
    ;    if m gt sxpar(h,'NAXIS3')-6 then high=0 else high=2
    PSF_EXTRACT,x_in,y_in,0,0,median(input[*,*,m-50:m+50],dimension=3),2*side+1,psf_new,Psf_fwhm,Background,/SKY_MEDIAN
    psf_new=PSF_GAUSSIAN( Npixel=101, FWHM=Psf_fwhm, /NORMAL )
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


  fits_write,root+'BUDDI/psf.fits',psf_cube_final,h
  stop
  
end