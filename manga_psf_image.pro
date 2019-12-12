; quick code to create PSF datacubes for MaNGA datacubes
; in preparation of fitting the datacube with BUDDI



pro MaNGA_PSF_image

readcol,'star_centres.txt',format='f,f',xstar,ystar,/silent


fits_read,'science_reduced_img_b.fits',input_b,header_b
fits_read,'science_reduced_img_i.fits',input_i,header_i


x_side=sxpar(header_b,'NAXIS1')
y_side=sxpar(header_b,'NAXIS2')


psf_out=fltarr(x_side,y_side)
psf_temp=fltarr(x_side,y_side)
side=50

PSF_EXTRACT,xstar,ystar,0,0,input_b,2*side+1,psf_new,Psf_fwhm,Background,/SKY_MEDIAN
fits_write,'psf_b.fits',psf_new,header_b
PSF_EXTRACT,xstar,ystar,0,0,input_i,2*side+1,psf_new,Psf_fwhm,Background,/SKY_MEDIAN
fits_write,'psf_i.fits',psf_new,header_i

end