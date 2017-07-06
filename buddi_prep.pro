; quick code to prepare datacubes for BUDDI.
; The code will:
;   - convert flux back into counts in a really dirty way (/counts)
;   - log-rebin the datacube in the wavelength direction if 
;     it's linearly binned (/log_rebin)
; 

;===============================================================


;===============================================================

pro BUDDI_prep,input_file,COUNTS=counts,LOG_REBIN=log_rebin,MANGA=manga,MUSE=muse,JY=jy,BADPIX=badpix
;buddi_prep,'BUDDI_input.txt',/MANGA,/JY,/BADPIX
;buddi_prep,'BUDDI_input2_r2.txt',/LOG_REBIN,/JY,/MUSE
;buddi_prep,'BUDDI_input.txt',/LOG_REBIN,/JY,/MUSE


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

;fits_read,directory+file+'.fits',input,h
fits_read,'/home/bhaeussl/MUSE_data/2MIG_131_DATACUBE_CLEANED.fits',input,h

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
      output[column,row,*]=new_spec
    endfor
  endfor
  sxaddpar,h,'CRVAL3',logLam2[0]
  sxaddpar,h,'CDELT3',logLam2[1]-logLam2[0]
  sxaddpar,h,'CD3_3',logLam2[1]-logLam2[0]
endif else output=input_counts

fits_write,directory+file+'_FLUX.fits',output,extname='FLUX'

temp=mrdfits(directory+file+'.fits',0,h_temp)
if keyword_set(Jy) and keyword_set(MaNGA) then sxaddpar,h_temp,'BUNIT','Jy', 'Specific intensity (per spaxel)'
sxaddpar,h_temp,'CD1_1',sxpar(h,'CD1_1')
sxaddpar,h_temp,'CD2_2',sxpar(h,'CD2_2')

modfits,directory+file+'_FLUX.fits',0,h_temp
modfits,directory+file+'_FLUX.fits',0,h,extname='FLUX'


if keyword_set(MANGA) then begin
  ;extract IVAR datacube and convert to sigma images
  IVAR=MRDFITS(directory+file+'.fits', 2, hdr2)
  s=size(IVAR)
  sigma=fltarr(s[1],s[2],s[3])
  sigma[*,*,*]=1/sqrt(IVAR[*,*,*])
  

  ;fits_read,directory+file+'.fits',input,h2,exten_no=2
  fits_read,'/home/bhaeussl/MUSE_data/2MIG_131_DATACUBE_CLEANED.fits',input,h2,exten_no=2
  if keyword_set(Jy) then begin
    ;for zz=0,z-1,1 do sigma[*,*,zz]=((sigma[*,*,zz]*1e-17)*wave[zz]*wave[zz]/(3.e5))/1e-23
      for zz=0,z-1,1 do sigma[*,*,zz]=(sigma[*,*,zz]*1e-17)*wave[zz]*wave[zz]*3.34e4

  endif
  fits_write,directory+file+'_SIGMA.fits',sigma,extname='SIGMA'
  sxaddpar,h2,'EXTNAME','SIGMA'

  temp=mrdfits(directory+file+'.fits',0,h_temp)
  if keyword_set(Jy) then sxaddpar,h_temp,'BUNIT','Jy', 'Specific intensity (per spaxel)'
  modfits,directory+file+'_SIGMA.fits',0,h_temp
  modfits,directory+file+'_SIGMA.fits',1,h2,extname='SIGMA'
endif


if keyword_set(MUSE) then begin
  ;extract IVAR datacube and convert to sigma images
  ;for MUSE data, the STAT extension shows sigma^2
  IVAR=MRDFITS(directory+file+'.fits', 2, hdr2)
  s=size(IVAR)
  sigma=fltarr(s[1],s[2],s[3])
  sigma[*,*,*]=sqrt(IVAR[*,*,*])
  
  
;  fits_read,directory+file+'.fits',input,h2,exten_no=2
  fits_read,'/home/bhaeussl/MUSE_data/2MIG_131_DATACUBE_CLEANED.fits',input,h2,exten_no=2
  if keyword_set(Jy) then begin
    ;for zz=0,z-1,1 do sigma[*,*,zz]=((sigma[*,*,zz]*1e-17)*wave[zz]*wave[zz]/(3.e5))/1e-23
    for zz=0,z-1,1 do sigma[*,*,zz]=(sigma[*,*,zz]*1e-20)*wavelength[zz]*wavelength[zz]*3.34e4
    
  endif
  fits_write,directory+file+'_SIGMA.fits',sigma,extname='SIGMA'
  sxaddpar,h2,'EXTNAME','SIGMA'
  
  temp=mrdfits(directory+file+'.fits',0,h_temp)
  if keyword_set(Jy) then sxaddpar,h_temp,'BUNIT','Jy', 'Specific intensity (per spaxel)'
  modfits,directory+file+'_SIGMA.fits',0,h_temp
  modfits,directory+file+'_SIGMA.fits',1,h2,extname='SIGMA'
endif




;badpix=fltarr(x,y,z)
if keyword_set(BADPIX) then begin
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
    sxaddpar,hdr3,'EXTNAME','BADPIX'
    fits_write,directory+file+'_BADPIX.fits',badpix,extname='BADPIX'
    modfits,directory+file+'_BADPIX.fits',0,h_temp
    modfits,directory+file+'_BADPIX.fits',0,h,extname='BADPIX'
  endif
endif


print,'#############################'
if keyword_set(MANGA) then begin
  print,'## CRVAL3='+string(wavelength[0])+' ##'
  print,'## CDELT3='+string(wavelength[1]-wavelength[0])+' ##'
endif
if keyword_set(MUSE) then begin
  print,'## CRVAL3='+string(LogLam2[0])+' ##'
  print,'## CDELT3='+string(LogLam2[1]-LogLam2[0])+' ##'
endif
print,'#############################'

if keyword_set(MANGA) then MaNGA_PSF_datacube,input_file
stop
end