;create a badpiel mask from a region file

;cd Data_Reduction
;idl
;.compile badpixelmask_MUSE
;badpixelmask_MUSE

;save regions as ds9 format, x and y (or image)


pro badpixelmask_MUSE

dir='/raid/ejohnston/MUSE_S0s/CCC158/IFU_decomp/decomposed_data/images/'
;dir='/data2/ejohnston/Fornax_new_sample/FCC306/'
;dir='/raid/ejohnston/double_nucleated/0338/'
;reg='BUDDI/ds9.reg'
reg='ds9.reg'
;reg='products/Galfit_fits/i/ds9_circles.reg'
;im='products/IMAGE_FOV_0007'
im='orig'
;im='IFU_decomp_2comp/median_image/image'
image = mrdfits(dir+im+'.fits',1)
hdr = headfits(dir+im+'.fits',exten=1)
;mask=fltarr(sxpar(hdr,'NAXIS1'),sxpar(hdr,'NAXIS2'))
image[*,*]=0
newmask=MASKREGION(image, dir+reg,BADVALUE=1)
;newmask[292:*,*]=1   ;NGC1396
;newmask[0:325,*]=1    ;CCC43


;fits_read,'/data2/ejohnston/MUSE_Es/J020536_reduced/BUDDI/IFU_decomp/median_image/residual.fits',resid
;emission=where(resid[*,191:*] gt 2e-8)
;s = SIZE(resid)
;ncol = s(1)
;col = emission MOD ncol
;row = emission / ncol
;newmask[col,row+191]=1


;fits_write,dir+'BUDDI/badpix_2D.fits',newmask,hdr
fits_write,dir+'badpix_2D.fits',newmask,hdr
;fits_write,dir+'badpix_Ha_residual.fits',newmask,hdr
;fits_write,dir+'products/Galfit_fits/i/badpix_2D_new.fits',newmask,hdr


masked=where(newmask eq 1,count)
s = SIZE(image)
ncol = s(1)
col = masked MOD ncol
row = masked / ncol
;openw,01,dir+'BUDDI/badpix_new.txt'
openw,01,dir+'badpix_new.txt'
for n=0,count-1,1 do printf,01,col[n],row[n]


close,01

end


