;create mock datacubes to test BUDDI on nucleated dwarf galaxies.

pro mock_MUSE_cubes

;read in stellar populations for the models.
; best to use some variations to make the points stand out better in the paper plots
Miles='/raid/ejohnston/IDLWorkspace84/ppxf/miles_models/'
Miles2='/raid/ejohnston/IDLWorkspace84/PPXF_v4_79/miles_models_full/'
fits_read,Miles+'Mun1.30Zm1.31T03.1623_iPp0.00_baseFe_linear_FWHM_2.51.fits',age_y1_met_l,h_spec
fits_read,Miles+'Mun1.30Zm0.71T03.9811_iPp0.00_baseFe_linear_FWHM_2.51.fits',age_y2_met_s
fits_read,Miles+'Mun1.30Zp0.00T02.5119_iPp0.00_baseFe_linear_FWHM_2.51.fits',age_y3_met_h;
fits_read,Miles+'Mun1.30Zm1.71T07.9433_iPp0.00_baseFe_linear_FWHM_2.51.fits',age_s1_met_l
;fits_read,Miles+'Mun1.30Zm0.71T07.9433_iPp0.00_baseFe_linear_FWHM_2.51.fits',age_s2_met_s
fits_read,Miles2+'Mun1.30Zm0.40T07.0795_iPp0.00_baseFe.fits',age_s2_met_s
fits_read,Miles+'Mun1.30Zp0.00T07.9433_iPp0.00_baseFe_linear_FWHM_2.51.fits',age_s3_met_h
;fits_read,Miles+'Mun1.30Zm1.71T15.8489_iPp0.00_baseFe_linear_FWHM_2.51.fits',age_o1_met_l;
fits_read,Miles2+'Mun1.30Zm0.71T11.2202_iPp0.00_baseFe.fits',age_o1_met_l;
fits_read,Miles+'Mun1.30Zm0.71T12.5893_iPp0.00_baseFe_linear_FWHM_2.51.fits',age_o2_met_s
fits_read,Miles+'Mun1.30Zp0.22T10.0000_iPp0.00_baseFe_linear_FWHM_2.51.fits',age_o3_met_h

fits_read,Miles+'Mun1.30Zm1.71T12.5893_iPp0.00_baseFe_linear_FWHM_2.51.fits',age_o1_met_l2;
fits_read,Miles+'Mun1.30Zm0.40T06.3096_iPp0.00_baseFe_linear_FWHM_2.51.fits',age_y2_met_s2


;log-rebin the model spectra
lamRange2 = sxpar(h_spec,'CRVAL1') + [0d,sxpar(h_spec,'CDELT1')*(sxpar(h_spec,'NAXIS1')-1d)]
log10_rebin, lamRange2, age_y1_met_l, age_y1_met_l_log, logLam2;, VELSCALE=velScale
log10_rebin, lamRange2, age_y2_met_s, age_y2_met_s_log, logLam2;, VELSCALE=velScale
log10_rebin, lamRange2, age_y3_met_h, age_y3_met_h_log, logLam2;, VELSCALE=velScale
log10_rebin, lamRange2, age_s1_met_l, age_s1_met_l_log, logLam2;, VELSCALE=velScale
log10_rebin, lamRange2, age_s2_met_s, age_s2_met_s_log, logLam2;, VELSCALE=velScale
log10_rebin, lamRange2, age_s3_met_h, age_s3_met_h_log, logLam2;, VELSCALE=velScale
log10_rebin, lamRange2, age_o1_met_l, age_o1_met_l_log, logLam2;, VELSCALE=velScale
log10_rebin, lamRange2, age_o2_met_s, age_o2_met_s_log, logLam2;, VELSCALE=velScale
log10_rebin, lamRange2, age_o3_met_h, age_o3_met_h_log, logLam2;, VELSCALE=velScale

log10_rebin, lamRange2, age_o1_met_l2, age_o1_met_l2_log, logLam2;, VELSCALE=velScale
log10_rebin, lamRange2, age_y2_met_s2, age_y2_met_s2_log, logLam2;, VELSCALE=velScale
print,logLam2[0],logLam2[-1]

wave0=logLam2[0]
step=logLam2[1]-logLam2[0]
wave=fltarr(n_elements(logLam2))
for n=0,n_elements(logLam2)-1,1 do wave[n]=(wave0+n*step)
wave_lin=10^wave
sample=where(wave_lin ge 4800)

age_y1_met_l_log/=median(age_y1_met_l_log[sample])
age_y2_met_s_log/=median(age_y2_met_s_log[sample])
age_y3_met_h_log/=median(age_y3_met_h_log[sample])
age_s1_met_l_log/=median(age_s1_met_l_log[sample])
age_s2_met_s_log/=median(age_s2_met_s_log[sample])
age_s3_met_h_log/=median(age_s3_met_h_log[sample])
age_o1_met_l_log/=median(age_o1_met_l_log[sample])
age_o2_met_s_log/=median(age_o2_met_s_log[sample])
age_o3_met_h_log/=median(age_o3_met_h_log[sample])

age_o1_met_l2_log/=median(age_o1_met_l2_log[sample])
age_y2_met_s2_log/=median(age_y2_met_s2_log[sample])

;read in galfit models for each component, and a residuals and sky residuals datacube for FCC245?
fits_read,'/raid/ejohnston/FCC245/BUDDI/IFU_decomp/median_image/temp/subcomps.fits',c1_image,h_c1,extname='COMPONENT_1_sersic'
fits_read,'/raid/ejohnston/FCC245/BUDDI/IFU_decomp/median_image/temp/subcomps.fits',c2_image,h_c2,extname='COMPONENT_2_psf'

fits_read,'/raid/ejohnston/FCC245/BUDDI/IFU_decomp/decomposed_data/residuals_cube.fits',resid,h_resid
fits_read,'/raid/ejohnston/FCC245/BUDDI/IFU_decomp/decomposed_data/residual_sky_cube.fits',resid_sky,h_sky

s=size(resid)
c1_cube=fltarr(s[1],s[2],2500);s[3])
c2_cube=fltarr(s[1],s[2],2500);s[3])
;for j=0,s[3]-1,1 do begin


for j=0,2499,1 do begin
  c1_cube[*,*,j]=c1_image*age_o1_met_l_log[sample[j]]
  c2_cube[*,*,j]=c2_image*age_y3_met_h_log[sample[j]]
endfor
Mock_Host_Mlower_Aolder=resid+resid_sky+c1_cube+c2_cube


for j=0,2499,1 do begin
  c1_cube[*,*,j]=c1_image*age_y1_met_l_log[sample[j]]
  c2_cube[*,*,j]=c2_image*age_o3_met_h_log[sample[j]]
endfor
Mock_Host_Mlower_Ayounger=resid+resid_sky+c1_cube+c2_cube


for j=0,2499,1 do begin
  c1_cube[*,*,j]=c1_image*age_s1_met_l_log[sample[j]]
  c2_cube[*,*,j]=c2_image*age_s3_met_h_log[sample[j]]
endfor
Mock_Host_Mlower_Asame=resid+resid_sky+c1_cube+c2_cube


for j=0,2499,1 do begin
  c1_cube[*,*,j]=c1_image*age_o3_met_h_log[sample[j]]
  c2_cube[*,*,j]=c2_image*age_y1_met_l_log[sample[j]]
endfor
Mock_Host_Mhigher_Aolder=resid+resid_sky+c1_cube+c2_cube


for j=0,2499,1 do begin
  c1_cube[*,*,j]=c1_image*age_y3_met_h_log[sample[j]]
  c2_cube[*,*,j]=c2_image*age_o1_met_l_log[sample[j]]
endfor
Mock_Host_Mhigher_Ayounger=resid+resid_sky+c1_cube+c2_cube


for j=0,2499,1 do begin
  c1_cube[*,*,j]=c1_image*age_s3_met_h_log[sample[j]]
  c2_cube[*,*,j]=c2_image*age_s1_met_l_log[sample[j]]
endfor
Mock_Host_Mhigher_Asame=resid+resid_sky+c1_cube+c2_cube


for j=0,2499,1 do begin
  c1_cube[*,*,j]=c1_image*age_o2_met_s_log[sample[j]]
  c2_cube[*,*,j]=c2_image*age_y2_met_s_log[sample[j]]
endfor
Mock_Host_Msame_Aolder=resid+resid_sky+c1_cube+c2_cube


for j=0,2499,1 do begin
  c1_cube[*,*,j]=c1_image*age_y2_met_s_log[sample[j]]
  c2_cube[*,*,j]=c2_image*age_o2_met_s_log[sample[j]]
endfor
Mock_Host_Msame_Ayounger=resid+resid_sky+c1_cube+c2_cube


for j=0,2499,1 do begin
  c1_cube[*,*,j]=c1_image*age_s2_met_s_log[sample[j]]
  c2_cube[*,*,j]=c2_image*age_s2_met_s_log[sample[j]]
endfor
Mock_Host_Msame_Asame=resid+resid_sky+c1_cube+c2_cube


for j=0,2499,1 do begin
  c1_cube[*,*,j]=c1_image*age_o1_met_l2_log[sample[j]]
  c2_cube[*,*,j]=c2_image*age_y2_met_s2_log[sample[j]]
endfor
Mock_Host_extra=resid+resid_sky+c1_cube+c2_cube



;convert the Galfit models to datacubes by adding the spectra

sxaddpar,h_resid,'CRVAL1',wave[sample[0]]
sxaddpar,h_resid,'CDELT1',step
sxaddpar,h_resid,'CD1_1',step

;save datacubes to their own diretctories
fits_write,'/raid/ejohnston/mock_cubes/1_Mock_Host_Mlower_Aolder.fits',Mock_Host_Mlower_Aolder,h_resid
;fits_write,'/raid/ejohnston/mock_cubes/2_Mock_Host_Mlower_Ayounger.fits',Mock_Host_Mlower_Ayounger,h_resid
;fits_write,'/raid/ejohnston/mock_cubes/3_Mock_Host_Mlower_Asame.fits',Mock_Host_Mlower_Asame,h_resid
;fits_write,'/raid/ejohnston/mock_cubes/4_Mock_Host_Mhigher_Aolder.fits',Mock_Host_Mhigher_Aolder,h_resid
fits_write,'/raid/ejohnston/mock_cubes/5_Mock_Host_Mhigher_Ayounger.fits',Mock_Host_Mhigher_Ayounger,h_resid
;fits_write,'/raid/ejohnston/mock_cubes/6_Mock_Host_Mhigher_Asame.fits',Mock_Host_Mhigher_Asame,h_resid
;fits_write,'/raid/ejohnston/mock_cubes/7_Mock_Host_Msame_Aolder.fits',Mock_Host_Msame_Aolder,h_resid
;fits_write,'/raid/ejohnston/mock_cubes/8_Mock_Host_Msame_Ayounger.fits',Mock_Host_Msame_Ayounger,h_resid
fits_write,'/raid/ejohnston/mock_cubes/9_Mock_Host_Msame_Asame.fits',Mock_Host_Msame_Asame,h_resid

;fits_write,'/raid/ejohnston/mock_cubes/0_Mock_Host_extra.fits',Mock_Host_extra,h_resid

stop



end