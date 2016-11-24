; 
; This code will print out the GalfitM feedme files for multi- 
; band fits for IFU data. Note that GALFITM is being used fo all 
; fits to maintain consistency in the results and in the codes.
; 
;  BINNED keyword- to be used with the binned images
;  SLICES keyword- to be used with individual image slices
; 
; scale
;
pro galfitm_multiband_psf,output,x,y,scale,$
  estimates,disk_mag_polynomial_in,galfitm,BINNED=binned,$
  SLICES=slices,FILE=file,HEADER=header

ugriz=[4770,6231,7625,9134]
PSF_files=['Gpsf.fits','Rpsf.fits','Ipsf.fits','Zpsf.fits']

disk_type='psf'

x1=4
nband=4;info[2]
imgblock='imgblock_psf'




openw,60,output+'galfitm.feedme'

printf,60,'==============================================================================='
printf,60,'# IMAGE and GALFIT CONTROL PARAMETERS';+$
printf, 60, 'A) Gpsf.fits,Rpsf.fits,Ipsf.fits,Zpsf.fits           # Input data image (FITS file)'
printf, 60, 'A1) 1,2,3,4             # Band labels (can be omitted if fitting a single band)'
printf, 60, 'A2) 4770,6231,7625,9134             # Band wavelengths'
printf, 60, 'B) imgblock_psf.fits       # Output data image block
printf, 60, 'C) none                # Sigma image name (made from data if blank or "none") 
printf, 60, 'D) none,none,none,none           # Input PSF image and (optional) diffusion kernel
printf, 60, 'E) 1                   # PSF fine sampling factor relative to data 
printf, 60, 'F) none,none,none,none                # Bad pixel mask (FITS image or ASCII coord list)
printf, 60, 'G) none                # File with parameter constraints (ASCII file)'  
printf, 60, 'H) 1    72   1  72    # Image region to fit (xmin xmax ymin ymax)'
printf, 60, 'I) 72 72      # Size of the convolution box (x y)'
printf, 60, 'J) 5,5,5,5              # Magnitude photometric zeropoint '
printf, 60, 'K) 0.5  0.5        # Plate scale (dx dy)    [arcsec per pixel]'
printf, 60, 'O) regular             # Display type (regular, curses, both)'
printf, 60, 'P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps'
     
printf, 60, ' '
printf, 60, ' '
      
printf, 60, '  # INITIAL FITTING PARAMETERS'
printf, 60, '#'
printf, 60, '#   For object type, the allowed functions are: '
printf, 60, '#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat, '
printf, 60, '#       ferrer, powsersic, sky, and isophote. '
printf, 60, '#  '
printf, 60, '#   Hidden parameters will only appear when theyre specified:'
printf, 60, '#       C0 (diskyness/boxyness), '
printf, 60, '#       Fn (n=integer, Azimuthal Fourier Modes),'
printf, 60, '#       R0-R10 (PA rotation, for creating spiral structures).'
printf, 60, '# '
printf, 60, '# -----------------------------------------------------------------------------'
printf, 60, '#   par)    par value(s)    fit toggle(s)    # parameter description '
printf, 60, '# -----------------------------------------------------------------------------'

printf, 60, ' '
printf, 60, ' '
printf, 60, ' '
printf, 60, ' '
    

       
printf, 60, '# Object number: 1'    ;disc
printf, 60, ' 0) psf                 #  object type'
printf, 60, ' 1) '+x+','+x+','+x+','+x+'   1 band  #  position x, y'
printf, 60, ' 2) '+y+','+y+','+y+','+y+'   1 band  #  position x, y'
printf, 60, ' 3) '+estimates[1]+','+estimates[1]+','+estimates[1]+','+estimates[1]+'        3 band  #  Integrated magnitude' 


close,60
end

