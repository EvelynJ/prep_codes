; quick code to prepare datacubes for BUDDI.
; The code will:
;   - convert flux back into counts in a really dirty way (/counts)
;   - log-rebin the datacube in the wavelength direction if 
;     it's linearly binned (/log_rebin)
; 

;===============================================================


;===============================================================

pro BUDDI_prep,input_file,COUNTS=counts,LOG_REBIN=log_rebin,MANGA=manga,MUSE=muse,JY=jy,BADPIX=badpix,PSF=psf
;buddi_prep,'BUDDI_input.txt',/MANGA,/JY,/BADPIX
;buddi_prep,'BUDDI_input2_r2.txt',/LOG_REBIN,/JY,/MUSE
;buddi_prep,'BUDDI_input.txt',/LOG_REBIN,/JY,/MUSE,/BADPIX,/PSF
;buddi_prep,'BUDDI_input.txt',/LOG_REBIN  <= Victors cubes

same_dir='y'  ;yes or no- is the original file int he same directory?
file_temp='/raid/ejohnston/Luis/RCS0327_16mc_zap.fits.fits'

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
PSF_cube=setup.psf_file

directory=root


if same_dir eq 'y' then fits_read,directory+file+'.fits',input,h
if same_dir eq 'n' then fits_read,file_temp,input,h

;h = headfits(directory+file+'.fits')
;h2 = headfits(directory+file+'.fits',exten=1)
x=sxpar(h,'NAXIS1')
y=sxpar(h,'NAXIS2')
z=sxpar(h,'NAXIS3')


;write out wavelength array
if keyword_set(MANGA) then begin
    newfile = directory+file+'.fits'
    wave = MRDFITS(newfile, 4, hdr4)    ; wavelength in Angstrom
    wavelength=alog10(wave)
        
    
    sxaddpar,h,'CRVAL3',wavelength[0]
    sxaddpar,h,'CDELT3',wavelength[1]-wavelength[0]
    sxaddpar,h,'CD3_3',wavelength[1]-wavelength[0]
    

endif else begin
    wavelength0=sxpar(h,'CRVAL3')
    step=sxpar(h,'CD3_3')
    wavelength=fltarr(sxpar(h,'NAXIS3'))
    for m=0,sxpar(h,'NAXIS3')-1,1 do wavelength[m]=wavelength0+(m*step)
    
endelse


;print,'New wavelength solution for BUDDI input file:'
;print,'D00) wavelength start= '+string(wavelength[0])
;print,'D01) wavelength step=  '+string(wavelength[1]-wavelength[0])




input_counts=fltarr(x,y,z)
if keyword_set(Jy) then begin
  print,'*** Now converting units to Jy'
  
  for zz=0,z-1,1 do begin
    ;convert from 10^-17 erg/s/cm2/AA to Jy
;    input_counts[*,*,zz]=((input[*,*,zz]*1e-17)*wave[zz]*wave[zz]/(3.e5))/1e-23
    if keyword_set(MANGA) then input_counts[*,*,zz]=(input[*,*,zz]*1e-17)*wave[zz]*wave[zz]*3.34e4
    if keyword_set(MUSE) then input_counts[*,*,zz]=(input[*,*,zz]*1e-20)*wavelength[zz]*wavelength[zz]*3.34e4
  endfor
  sxaddpar,h,'BUNIT','Jy', 'Specific intensity (per spaxel)'

endif else input_counts=input







output=fltarr(x,y,z)
if keyword_set(LOG_REBIN) then begin
  print,'*** Now rebinning logarithmically'
  lamRange2=[wavelength[0],wavelength[-1]]
  
  for column=0,x-1,1 do begin
    for row=0,y-1,1 do begin
      temp_spec=fltarr(z)
      temp_spec[*]=input_counts[column,row,*]
      log10_rebin, lamRange2, temp_spec, new_spec, logLam2, VELSCALE=velScale
      zz=n_elements(new_spec)-1
      output[column,row,0:zz]=new_spec
    endfor
  endfor
  sxaddpar,h,'CRVAL3',logLam2[0]
  sxaddpar,h,'CDELT3',logLam2[1]-logLam2[0]
  sxaddpar,h,'CD3_3',logLam2[1]-logLam2[0]
  
  if keyword_set(PSF) then begin
    fits_read,directory+PSF_file,psf_input,h_psf
    x_psf=sxpar(h_psf,'NAXIS1')
    y_psf=sxpar(h_psf,'NAXIS2')
    z_psf=sxpar(h_psf,'NAXIS3')
    output_psf=fltarr(x_psf,y_psf,z_psf)
    for column=0,x_psf-1,1 do begin
      for row=0,y_psf-1,1 do begin
        temp_spec=fltarr(z_psf)
        temp_spec[*]=psf_input[column,row,*]
        log10_rebin, lamRange2, temp_spec, new_spec, logLam2, VELSCALE=velScale
        zz=n_elements(new_spec)-1
        output_psf[column,row,0:zz]=new_spec
      endfor
    endfor
    sxaddpar,h_psf,'CRVAL3',logLam2[0]
    sxaddpar,h_psf,'CDELT3',logLam2[1]-logLam2[0]
    sxaddpar,h_psf,'CD3_3',logLam2[1]-logLam2[0]
    fits_write,directory+'psf_LOG.fits',output_psf,h_psf;,extname='FLUX'
  endif
endif else output=input_counts


print,'New wavelength solution for BUDDI input file:'
print,'D00) wavelength start= '+string(logLam2[0])
print,'D01) wavelength step=  '+string(logLam2[1]-logLam2[0])


fits_write,directory+file+'_FLUX.fits',output,h;,extname='FLUX'
stop
if same_dir eq 'y' then temp=mrdfits(directory+file+'.fits',0,h_temp)
;;if same_dir eq 'n' then h_temp=headfits(file_temp,exten=0)
;
if keyword_set(Jy) then sxaddpar,h_temp,'BUNIT','Jy', 'Specific intensity (per spaxel)'
sxaddpar,h_temp,'CD1_1',sxpar(h,'CD1_1')
sxaddpar,h_temp,'CD2_2',sxpar(h,'CD2_2')



if keyword_set(MANGA) then begin
  ;extract IVAR datacube and convert to sigma images
  IVAR=MRDFITS(directory+file+'.fits', 2, hdr2)
  s=size(IVAR)
  sigma=fltarr(s[1],s[2],s[3])
  sigma[*,*,*]=1/sqrt(IVAR[*,*,*])
  

  if same_dir eq 'y' then fits_read,directory+file+'.fits',input,h2,exten_no=2
  if same_dir eq 'n' then fits_read,file_temp,input,h2,exten_no=2
  if keyword_set(Jy) then begin
    ;for zz=0,z-1,1 do sigma[*,*,zz]=((sigma[*,*,zz]*1e-17)*wave[zz]*wave[zz]/(3.e5))/1e-23
      for zz=0,z-1,1 do sigma[*,*,zz]=(sigma[*,*,zz]*1e-17)*wave[zz]*wave[zz]*3.34e4

  endif
  fits_write,directory+file+'_SIGMA.fits',sigma,extname='SIGMA'
  sxaddpar,h2,'EXTNAME','SIGMA'

  if same_dir eq 'y' then temp=mrdfits(directory+file+'.fits',0,h_temp)
;  if same_dir eq 'n' then h_temp=headfits(file_temp,exten=0)
  
  if keyword_set(Jy) then sxaddpar,h_temp,'BUNIT','Jy', 'Specific intensity (per spaxel)'
  modfits,directory+file+'_SIGMA.fits',0,h_temp
  modfits,directory+file+'_SIGMA.fits',1,h2,extname='SIGMA'
endif






;==================================
;sigma datacube



if keyword_set(MUSE) then begin
  ;extract IVAR datacube and convert to sigma images
  ;for MUSE data, the STAT extension shows sigma^2
  
  if same_dir eq 'y' then IVAR=MRDFITS(directory+file+'.fits', 2, hdr2)
  if same_dir eq 'n' then IVAR=MRDFITS(file_temp, 2, hdr2)
  s=size(IVAR)
  sigma=fltarr(s[1],s[2],s[3])
  sigma[*,*,*]=sqrt(IVAR[*,*,*])
  
  
  if same_dir eq 'y' then fits_read,directory+file+'.fits',input,h2,exten_no=2
  if same_dir eq 'n' then fits_read,file_temp,input,h2,exten_no=2
  if keyword_set(Jy) then begin
    ;for zz=0,z-1,1 do sigma[*,*,zz]=((sigma[*,*,zz]*1e-17)*wave[zz]*wave[zz]/(3.e5))/1e-23
    for zz=0,z-1,1 do sigma[*,*,zz]=(sigma[*,*,zz]*1e-20)*wavelength[zz]*wavelength[zz]*3.34e4
    
  endif
  print,'log-rebinning the sigma cube'
  output_sigma=fltarr(x,y,z)
  for column=0,x-1,1 do begin
    for row=0,y-1,1 do begin
      temp_spec=fltarr(z)
      temp_spec[*]=sigma[column,row,*]
      log10_rebin, lamRange2, temp_spec, new_spec, logLam2, VELSCALE=velScale
      zz=n_elements(new_spec)-1
      output_sigma[column,row,0:zz]=new_spec
    endfor
  endfor
 
  fits_write,directory+file+'_SIGMA.fits',output_sigma,extname='SIGMA'
  sxaddpar,h2,'EXTNAME','SIGMA'
  
  if same_dir eq 'y' then temp=mrdfits(directory+file+'.fits',0,h_temp)
;  if same_dir eq 'n' then h_temp=headfits(file_temp,exten=0)
  
  if keyword_set(Jy) then sxaddpar,h_temp,'BUNIT','Jy', 'Specific intensity (per spaxel)'
  modfits,directory+file+'_SIGMA.fits',0,h_temp
  modfits,directory+file+'_SIGMA.fits',1,h2,extname='SIGMA'
endif







;==================================
;bad pixel datacube



badpix=fltarr(x,y,z)
if keyword_set(BADPIX) then begin
  print,'*** Now creating bad pixel datacube'
  if keyword_set(MANGA) then begin
    print,'creating bad pixel mask'
    badpix=MRDFITS(directory+file+'.fits', 3, hdr3)
    sxaddpar,hdr3,'EXTNAME','BADPIX'
    fits_write,directory+file+'_BADPIX.fits',badpix,extname='BADPIX'
    ;mwrfits,0,directory+file+'_BADPIX.fits',hdr3
    temp=mrdfits(directory+file+'.fits',0,h_temp)
    modfits,directory+file+'_BADPIX.fits',0,h_temp
    modfits,directory+file+'_BADPIX.fits',1,hdr3,extname='BADPIX'
  endif else begin
    ;where no bad pixel data cube provided, create one by 
    ;identifying 0-value or NAN pixels
    fits_read,directory+file+'_FLUX.fits',output
    index = WHERE(output EQ 0 or finite(output) eq 0)
    s = SIZE(output)
    ncol = s[1]
    nrow = s[2]
    col = index mod ncol
    row = (index / ncol) mod nrow
    frame = index / (nrow*ncol)
    badpix=fltarr(s[1],s[2],s[3])
    badpix[*,*,*]=0
    badpix[col,row,frame]=1
    badpix_orig=badpix
    sxaddpar,hdr3,'EXTNAME','BADPIX'

    fits_write,directory+file+'_BADPIX.fits',badpix_orig,extname='BADPIX'
    modfits,directory+file+'_BADPIX.fits',0,h_temp
    modfits,directory+file+'_BADPIX.fits',1,h,extname='BADPIX'
    
    ;need to modify the datacube to mask out the emission regions at the affected wavelengths
    if file_test(directory+'emission_image.fits') eq 1 then begin
      fits_read,directory+'emission_image.fits',emission,h
      limit=5e-8   ;2mig445
      limit=5e-8   ;FCC207
      emission_pixels=where(emission gt 5e-8)
      s = SIZE(emission)
      ncol = s(1)
      colX = emission_pixels MOD ncol
      rowX = emission_pixels / ncol
      em_im=emission
      em_im[*,*]=0
      em_im[colX,rowX]=1
;      em_regions=[229,254,345,360,1629,1676,1767,1797]  ;2mig445
      em_regions=[150,157,310,318,1768,1774,1780,1787,1797,1804,1907,1912,1918,1923,3539,3545]  ;FCC207
      for n=0,n_elements(em_regions)-1,2 do begin
        for m=em_regions[n],em_regions[n+1],1 do badpix[*,*,m]=badpix_orig[*,*,m]+em_im
      endfor
      
      fits_write,directory+file+'_BADPIX_EMISSION.fits',badpix,extname='BADPIX'
      modfits,directory+file+'_BADPIX_EMISSION.fits',0,h_temp
      modfits,directory+file+'_BADPIX_EMISSION.fits',1,h,extname='BADPIX'
      fits_write,directory+file+'tempytemp.fits',em_im



      fits_read,directory+'residual_Ha.fits',em_im,h
      index=where(em_im ge 5e-8)
      s=size(em_im)
      ncol = s(1)
      col = index MOD ncol
      row = index / ncol
      em_im_badpix=fltarr(s[1],s[2])
      em_im_badpix[col,row]=1
      em_im_badpix[*,0:140]=0
      em_im_badpix[*,185:*]=0
      em_im_badpix[0:218,*]=0
      em_im_badpix[241:*,*]=0
      badpix[*,*,*]=0
      
      em_regions=[150,157,310,318,1768,1774,1780,1787,1797,1804,1907,1912,1918,1923,3539,3545]  ;FCC207
      for n=0,n_elements(em_regions)-1,2 do begin
        for m=em_regions[n],em_regions[n+1],1 do badpix[*,*,m]=badpix_orig[*,*,m]+em_im_badpix
      endfor
      
      fits_write,directory+file+'_BADPIX_Ha_residual.fits',badpix,extname='BADPIX'
      modfits,directory+file+'_BADPIX_Ha_residual.fits',0,h_temp
      modfits,directory+file+'_BADPIX_Ha_residual.fits',1,h,extname='BADPIX'
    endif
    
    
  endelse
endif


print,'#############################'
if keyword_set(MANGA) then begin
  print,'## CRVAL3='+string(wavelength[0])+' ##'
  print,'## CDELT3='+string(wavelength[1]-wavelength[0])+' ##'
endif
;if keyword_set(MUSE) then begin
  print,'## CRVAL3='+string(LogLam2[0])+' ##'
  print,'## CDELT3='+string(LogLam2[1]-LogLam2[0])+' ##'
;endif
print,'#############################'

if keyword_set(MANGA) then MaNGA_PSF_datacube,input_file

stop
end