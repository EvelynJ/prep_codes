; quick code to create PSF datacubes for MaNGA datacubes
; in preparation of fitting the datacube with BUDDI


FUNCTION maggies2cts, maggies, expt, zp
;zp: zeropoint is a positive number
   return, maggies * expt / 10^(-0.4*zp)
END

pro MaNGA_PSF_datacube,input_file

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

;Best way to do this for SDSS MaNGA
;data is to read in the FWHM for each waveband, and then interpolate 
;between their central wavelengths.
;
;read in PSF extensions for each filter

fits_read,directory+file+'.fits',input_g,header_g,extname='GPSF'
fits_read,directory+file+'.fits',input_r,header_r,extname='RPSF'
fits_read,directory+file+'.fits',input_i,header_i,extname='IPSF'
fits_read,directory+file+'.fits',input_z,header_z,extname='ZPSF'

h = headfits(directory+file+'.fits')
h2 = headfits(directory+file+'.fits',exten=1)
sxaddpar,h,'EXPTIME',1
sxaddpar,h,'GAIN',1
x_side=sxpar(header_g,'NAXIS1')
y_side=sxpar(header_g,'NAXIS2')
z_side=sxpar(h2,'NAXIS3')
ugriz_wave=[3543,4770,6231,7635,9134]



;for each image slice, use the wavelength and the filter transmission curves 
;to determine the fraction of light from each filter and coadd them in the 
;correct proportions.
;
readcol,'/home/ejohnsto/Manga/g.dat',format='F,X,F,X,X',wave_g,transmission_g,/SILENT
readcol,'/home/ejohnsto/Manga/r.dat',format='F,X,F,X,X',wave_r,transmission_r,/SILENT
readcol,'/home/ejohnsto/Manga/i.dat',format='F,X,F,X,X',wave_i,transmission_i,/SILENT
readcol,'/home/ejohnsto/Manga/z.dat',format='F,X,F,X,X',wave_z,transmission_z,/SILENT

wavelength0_lin=sxpar(h2,'CRVAL3')
wavelength0=alog10(wavelength0_lin)
step_lin=sxpar(h2,'CD3_3')
step=alog10(wavelength0_lin+step_lin)-wavelength0
wavelength_arr=fltarr(sxpar(h2,'NAXIS3'))
for m=0,sxpar(h2,'NAXIS3')-1,1 do wavelength_arr[m]=wavelength0+(m*step)


psf_out=fltarr(x_side,y_side,z_side)
psf_temp=fltarr(x_side,y_side)

for j=0,z_side-1,1 do begin
  wave=10^wavelength_arr[j]
  total_light=0
  if wave ge wave_g[0] and wave le wave_g[-1] then $
    total_g=linear_interpolate(wave,wave_g,transmission_g)$
    else total_g=0
  if wave ge wave_r[0] and wave le wave_r[-1] then $
    total_r=linear_interpolate(wave,wave_r,transmission_r)$
    else total_r=0
  if wave ge wave_i[0] and wave le wave_i[-1] then $
    total_i=linear_interpolate(wave,wave_i,transmission_i)$
    else total_i=0
  if wave ge wave_z[0] and wave le wave_z[-1] then $
    total_z=linear_interpolate(wave,wave_z,transmission_z)$
    else total_z=0
  total_light=total_g+total_r+total_i+total_z

  for x=0,x_side-1,1 do begin
    for y=0,y_side-1,1 do psf_temp[x,y]=((total_g/total_light)*input_g[x,y])+$
      ((total_r/total_light)*input_r[x,y])+((total_i/total_light)*input_i[x,y])+$
      ((total_z/total_light)*input_z[x,y])
  endfor
  psf_out[*,*,j]=psf_temp
  
endfor


;spawn,'mkdir '+directory+'PSF_datacube_creation/'
;
;mkhdr,h,input_g
;sxaddpar,h,'CRPIX1',sxpar(header_g,'CRPIX1')
;sxaddpar,h,'CRPIX2',sxpar(header_g,'CRPIX2')
;sxaddpar,h,'CRVAL1',sxpar(header_g,'CRVAL1')
;sxaddpar,h,'CRVAL2',sxpar(header_g,'CRVAL2')
;sxaddpar,h,'CD1_1',sxpar(header_g,'CD1_1')
;sxaddpar,h,'CD2_2',sxpar(header_g,'CD2_2')
;sxaddpar,h,'CTYPE1',sxpar(header_g,'CTYPE1')
;sxaddpar,h,'CTYPE2',sxpar(header_g,'CTYPE2')
;sxaddpar,h,'CUNIT1',sxpar(header_g,'CUNIT1')
;sxaddpar,h,'CUNIT2',sxpar(header_g,'CUNIT2')
;sxaddpar,h,'IFURA',sxpar(header_g,'IFURA')
;sxaddpar,h,'IFUDEC',sxpar(header_g,'IFUDEC')
;sxaddpar,h,'OBJRA',sxpar(header_g,'OBJRA')
;sxaddpar,h,'OBJDEC',sxpar(header_g,'OBJDEC')
;sxaddpar,h,'BUNIT',sxpar(header_g,'BUNIT')
;sxaddpar,h,'GFWHM',sxpar(header_g,'GFWHM')
;sxaddpar,h,'RFWHM',sxpar(header_r,'RFWHM')
;sxaddpar,h,'IFWHM',sxpar(header_i,'IFWHM')
;sxaddpar,h,'ZFWHM',sxpar(header_z,'ZFWHM')
;sxaddpar,h,'EXPTIME',900
;
;for x=0,sxpar(header_g,'NAXIS1')-1 do begin
;  for y=0,sxpar(header_g,'NAXIS2')-1 do begin
;    input_g[x,y]=maggies2cts(input_g[x,y], 900., 5.)
;    input_r[x,y]=maggies2cts(input_r[x,y], 900., 5.)
;    input_i[x,y]=maggies2cts(input_i[x,y], 900., 5.)
;    input_z[x,y]=maggies2cts(input_z[x,y], 900., 5.)
;  endfor
;endfor
;tot_g=total(input_g)
;tot_r=total(input_r)
;tot_i=total(input_i)
;tot_z=total(input_z)
;;for x=0,sxpar(header_g,'NAXIS1')-1 do begin
;;  for y=0,sxpar(header_g,'NAXIS2')-1 do begin
;;    input_g[x,y]/=tot_g
;;    input_r[x,y]/=tot_r
;;    input_i[x,y]/=tot_i
;;    input_z[x,y]/=tot_z
;;  endfor
;;endfor
;
;fits_write,directory+'PSF_datacube_creation/Gpsf.fits',input_g,h
;fits_write,directory+'PSF_datacube_creation/Rpsf.fits',input_r,h
;fits_write,directory+'PSF_datacube_creation/Ipsf.fits',input_i,h
;fits_write,directory+'PSF_datacube_creation/Zpsf.fits',input_z,h
;
;x_centre='37'
;y_centre='37'
;scale=['0.5','0.5']
;estimates=['0','4']
;disk_mag_polynomial='3'
;galfitm='/Users/ejohnsto/Galfitm_files/galfitm-1.2.0-osx'
;galfitm_multiband_psf,directory+'PSF_datacube_creation/',x_centre,$
;          y_centre,scale,estimates, $
;          disk_mag_polynomial,$
;          galfitm,/file
;
;cd,directory+'PSF_datacube_creation/'
;
;spawn,galfitm+' galfitm.feedme'
;
;cd,'../'
;stop

mkhdr,h,psf_out
sxaddpar,h,'CRPIX1',sxpar(header_g,'CRPIX1')
sxaddpar,h,'CRPIX2',sxpar(header_g,'CRPIX2')
sxaddpar,h,'CRVAL1',sxpar(header_g,'CRVAL1')
sxaddpar,h,'CRVAL2',sxpar(header_g,'CRVAL2')
sxaddpar,h,'CD1_1',sxpar(header_g,'CD1_1')
sxaddpar,h,'CD2_2',sxpar(header_g,'CD2_2')
sxaddpar,h,'CTYPE1',sxpar(header_g,'CTYPE1')
sxaddpar,h,'CTYPE2',sxpar(header_g,'CTYPE2')
sxaddpar,h,'CUNIT1',sxpar(header_g,'CUNIT1')
sxaddpar,h,'CUNIT2',sxpar(header_g,'CUNIT2')
sxaddpar,h,'IFURA',sxpar(header_g,'IFURA')
sxaddpar,h,'IFUDEC',sxpar(header_g,'IFUDEC')
sxaddpar,h,'OBJRA',sxpar(header_g,'OBJRA')
sxaddpar,h,'OBJDEC',sxpar(header_g,'OBJDEC')
sxaddpar,h,'BUNIT',sxpar(header_g,'BUNIT')
sxaddpar,h,'GFWHM',sxpar(header_g,'GFWHM')
sxaddpar,h,'RFWHM',sxpar(header_r,'RFWHM')
sxaddpar,h,'IFWHM',sxpar(header_i,'IFWHM')
sxaddpar,h,'ZFWHM',sxpar(header_z,'ZFWHM')
sxaddpar,h,'EXPTIME',900



sxaddpar,h,'CRVAL3',wavelength_arr[0]
sxaddpar,h,'CDELT3',wavelength_arr[1]-wavelength_arr[0]
sxaddpar,h,'CD3_3',wavelength_arr[1]-wavelength_arr[0]

for n=0,sxpar(h,'NAXIS3')-1 do begin
  if total(total(psf_out[*,*,n])) gt 0 then psf_out[*,*,n]=psf_out[*,*,n]/total(psf_out[*,*,n])
endfor
fits_write,directory+'psf.fits', psf_out, h


end