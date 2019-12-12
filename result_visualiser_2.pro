


pro result_visualiser_2,setup,info,start_wavelength,end_wavelength,wavelength_in,$
  original_datacube,bestfit_datacube,residual_datacube,disk_datacube,$
  residual_sky_datacube,bulge_datacube,comp3_datacube,MANGA=manga,CALIFA=califa
  
  root=setup.root
  decomp=setup.decomp
  galaxy_ref=setup.galaxy_ref
  slices_dir=setup.slices_dir
  binned_dir=setup.binned_dir
  median_dir=setup.median_dir
  decomp_dir=setup.decomp_dir
  x_centre=fix(setup.x_centre-1)             ;x position of centre of galaxy, -1 to convert to position in array
  y_centre=fix(setup.y_centre-1)             ;y position of centre of galaxy, -1 to convert to position in array
  Redshift=setup.Redshift
  n_comp=setup.n_comp
  comp3_type=setup.comp3_type
  comp4_type=setup.comp4_type
  no_slices=setup.no_slices
  
  
  
  
first_image=info[0]
final_image=info[1]
no_bins=info[2]
images_per_bin=info[3]
start_wavelength=info[4]
end_wavelength=info[5]
total_images=final_image-first_image+1
no_images_final=total_images mod no_slices

wavelength_in=10^(wavelength_in)

; read in datacubes
fits_read,root+decomp+'decomposed_data/original_cube.fits',original_datacube,h_orig
fits_read,root+decomp+'decomposed_data/bestfit_cube.fits',bestfit_datacube,h_bestfit
fits_read,root+decomp+'decomposed_data/residuals_cube.fits',residual_datacube,h_resid

fits_read,root+decomp+'decomposed_data/component1_cube.fits',disk_datacube,h_disk
fits_read,root+decomp+'decomposed_data/residual_sky_cube.fits',residual_sky_datacube,h_sky
if n_comp ge 1100 then fits_read,root+decomp+'decomposed_data/component2_cube.fits',bulge_datacube,h_bulge
  if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then fits_read,root+decomp+'decomposed_data/comp3.fits',comp3_datacube,h_comp3
  if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011 then fits_read,root+decomp+'decomposed_data/comp4.fits',comp4_datacube,h_comp4

npix=sxpar(h_disk,'NAXIS3')
bulge_1D=fltarr(npix)
disk_1D=fltarr(npix)
comp3_1D=fltarr(npix)
comp4_1D=fltarr(npix)
orig_1D=fltarr(npix)
bestfit_1D=fltarr(npix)
resid_1D=fltarr(npix)
sky=fltarr(npix)
resid_sky_1D=fltarr(npix)

bulge_Re=fltarr(npix)
disk_Re=fltarr(npix)
wavelength=fltarr(npix)
for n=0,npix-1,1 do wavelength[n]=10^(sxpar(h_disk,'CRVAL3')+n*sxpar(h_disk,'CD3_3'))

;read in S/N values to identify good pixels for measuring flux
;readcol,root+galaxy_ref+'_S_N_array.txt',format='F,F,F,F',Xpix,Ypix,signal,noise,comment='#'
result=file_test(root+galaxy_ref+'_voronoi_2d_binning_output.txt')

if result eq 1 then begin
  readcol,root+galaxy_ref+'_voronoi_2d_binning_output.txt',format='F,F,F',Xpix,Ypix,BINpix,comment='#'
  
  
  ;calculate integrated bulge and disc flux
  ;first identify pixels affected by foreground/background stars so that their flux isn't included
  
  
  
  if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011 then begin
    x_temp=-long(x_centre-(comp4_x-1))
    y_temp=-long(y_centre-(comp4_y-1))
    temperoo_x=where(Xpix eq x_temp)
    temperoo_y=where(Ypix eq y_temp)
    match, temperoo_x,temperoo_y,x_element,y_element
    bad_bin=BINpix[temperoo_x[x_element]]
  endif else begin
    bad_bin=-99
    BINpix[*]=-999
  endelse
  
  FOR n=0,n_elements(Xpix)-1,1 do begin
      if BINpix[n] ne bad_bin then begin
        if n_comp ge 1100  then bulge_1D[*]+=bulge_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then comp3_1D[*]+=comp3_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011  then comp4_1D[*]+=comp4_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        disk_1D+=disk_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        orig_1D+=original_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        bestfit_1D+=bestfit_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        resid_1D+=residual_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        resid_sky_1D+=residual_sky_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        ;print,original_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        ;print,bulge_datacube[Xpix[n],Ypix[n],*]
      endif
  ENDFOR
endif else begin
  ;***not working yet!!!***
  for aaa=0,npix-1,1 do begin
    disk_1D[aaa]=total(disk_datacube[*,*,aaa])
    orig_1D[aaa]=total(original_datacube[*,*,aaa])
    bestfit_1D[aaa]=total(bestfit_datacube[*,*,aaa])
    resid_1D[aaa]=total(residual_datacube[*,*,aaa])
    resid_sky_1D[aaa]=total(residual_sky_datacube[*,*,aaa])
    if n_comp ge 1100  then bulge_1D[aaa]=total(bulge_datacube[*,*,aaa])
  endfor
endelse
;convert magnitude units into flux units where necessary, Magzp is 15 in the feedme files
;bulge_1D=10^((bulge_1D-15)/(-2.5))
;disk_1D=10^((disk_1D-15)/(-2.5))


;create wavelength array
;wavelength=fltarr(npix)
;wavelength0=sxpar(h_disk,'CRVAL3')
;step=sxpar(h_disk,'CD3_3')
;if keyword_set(manga) then begin
;  for m=0,npix-1,1 do wavelength[m]=10^(wavelength0+(step*m))
;  wavelength_log=alog10(wavelength)
;endif else if keyword_set(califa) then begin
;  for m=0,npix-1,1 do wavelength[m]=exp(wavelength0+(step*m))
;  wavelength_log=alog(wavelength)
;endif

temp = file_search(root+decomp+slices_dir+'galfitm_*.feedme',COUNT=nfiles)

for n=0,nfiles-1,1 do begin
  result=file_test(root+decomp+slices_dir+'subcomps_'+string(n,format='(I4.4)')+'.fits')
  if n ne nfiles-1 then nbands=no_slices else nbands=no_images_final   
  a=n*no_slices
  b=a+nbands-1
  if result eq 0 then begin
;    bulge_1D[a:b]=-99;bulge_1D[0]
;    disk_1D[a:b]=-99;disk_1D[0]
    bulge_Re[a:b]=-99
    disk_Re[a:b]=-99
  endif else begin
;    if n ne nfiles-1 then nbands=no_slices else nbands=no_images_final 
    
    

  if n_comp eq 1000 then res=read_sersic_results_2comp(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1100 then res=read_sersic_results_2comp(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1101 and comp4_type eq 'psf' then res=read_sersic_results_3psf(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1101 and comp4_type eq 'sersic' then res=read_sersic_results_3sersic(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1001 and comp4_type eq 'psf' then res=read_sersic_results_3psf(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1001 and comp4_type eq 'sersic' then res=read_sersic_results_3sersic(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  
  else if n_comp eq 1010  and comp3_type eq 'psf' then res=read_sersic_results_2comp_p(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1010  and comp3_type eq 'sersic' then res=read_sersic_results_2comp_s(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1110  and comp3_type eq 'psf' then res=read_sersic_results_2comp_p(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1110  and comp3_type eq 'sersic' then res=read_sersic_results_2comp_s(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'psf' and comp3_type eq 'psf' then res=read_sersic_results_3psf_p(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'psf' and comp3_type eq 'sersic' then res=read_sersic_results_3psf_s(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'sersic' and comp3_type eq 'psf' then res=read_sersic_results_3sersic_p(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'sersic' and comp3_type eq 'sersic' then res=read_sersic_results_3sersic_s(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1011 and comp4_type eq 'psf' and comp3_type eq 'psf' then res=read_sersic_results_3psf_p(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1011 and comp4_type eq 'psf' and comp3_type eq 'sersic' then res=read_sersic_results_3psf_s(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1011 and comp4_type eq 'sersic' and comp3_type eq 'psf' then res=read_sersic_results_3sersic_p(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1011 and comp4_type eq 'sersic' and comp3_type eq 'sersic' then res=read_sersic_results_3sersic_s(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) 


    exptime=sxpar(h_disk,'EXPTIME')
;    print,a,b
;    disk_1D[a:b]=10^((res.mag_galfit_band_d[0:nbands-1]-15)/(-2.5))*exptime
    disk_Re[a:b]=res.Re_galfit_band_d[0:nbands-1]
    sky[a:b]=res.sky_galfit_band[0:nbands-1]*n_elements(Xpix)
    
    if n_comp ge 1100 then begin
;      bulge_1D[a:b]=10^((res.mag_galfit_band_b[0:nbands-1]-15)/(-2.5))*exptime
      bulge_Re[a:b]=res.Re_galfit_band_b[0:nbands-1]
    endif
;    res=read_sersic_results(root+decomp+slices_dir+'models/subcomps_'+string(n,format='(I4.4)')+'.fits', nband, bd=1)
    
;    bulge_1D[n-first_image]=10^((res.mag_galfit_band_b-15)/(-2.5))*exptime
;    disk_1D[n-first_image]=10^((res.mag_galfit_band_d-15)/(-2.5))*exptime
    ;print,bulge_1D[n-944],res.mag_galfit_b
  endelse
endfor    

if n_comp eq 1000 or n_comp eq 1001 then disk_1d=disk_1d+sky $
else if n_comp eq 1100 or n_comp eq 1101 then begin
  disk_1D_orig=disk_1D
  bulge_1D_orig=bulge_1D
  disk_1D=disk_1d+(0.5*sky)
  bulge_1D=bulge_1D+(0.5*sky)
endif else if n_comp eq 1111  or n_comp eq 1110 then begin
  disk_1D_orig=disk_1D
  bulge_1D_orig=bulge_1D
  comp3_1D_orig=comp3_1D
  disk_1D=disk_1d+(0.333*sky)
  bulge_1D=bulge_1D+(0.333*sky)
  comp3_1D=comp3_1D+(0.333*sky)  
endif else begin
  disk_1D_orig=disk_1D
  comp3_1D_orig=comp3_1D
  disk_1D=disk_1d+(0.5*sky)
  comp3_1D=comp3_1D+(0.5*sky)  
endelse

;if galaxy_ref eq 'manga-7443-9102' then begin
;  for j=30,n_elements(disk_1d)-1,30 do begin
;    disk_1d[j]=0.5*(disk_1d[j-1]+disk_1d[j+1])
;    bulge_1d[j]=0.5*(bulge_1d[j-1]+bulge_1d[j+1])
;    comp3_1D[j]=0.5*(comp3_1D[j-1]+comp3_1D[j+1])
;  endfor
;endif

;obtain new bulge spectrum within 3Re
res2=read_sersic_results_2comp(root+decomp+'binned_images/imgblock.fits', nband, bd=1)
readcol,root+decomp+slices_dir+'info.txt',format='X,F',info
PA=res2.PA_GALFIT_BAND_D[0:nbands-1]
no_slices=info[2]
h_temp=headfits(root+decomp+'binned_images/image_'+string(0,format='(I4.4)')+'.fits')
wave1=sxpar(h_temp,'WAVELENG')
h_temp=headfits(root+decomp+'binned_images/image_'+string(no_slices-1,format='(I4.4)')+'.fits')
wave2=sxpar(h_temp,'WAVELENG')
wavelength_Re=5500
Re_b=chebeval(wavelength_Re,res2.RE_GALFIT_CHEB_B,INTERVAL=[wave1,wave2])
bulge_1D_small=fltarr(n_elements(bulge_1D))
;stop
measure_circular_radius,indgen(n_elements(Xpix)),Xpix,Ypix,0,0,PA[0], radii
FOR n=0,n_elements(Xpix)-1,1 do begin
    if n_comp ge 1100 and radii[n] le 3*Re_b then bulge_1D_small[*]+=bulge_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
endfor




set_plot,'ps'
;device,file='/Users/ejohnsto/Dropbox/papers/Paper4/decomposed_spectra_1D_blue_'+galaxy_ref+'.eps',xsize=19.5,ysize=11,/portrait;,/landscape
;device,file=root+decomp+'decomposed_data/Spectra_integrated_blue.eps',xsize=19.5,ysize=8,/portrait;,/landscape;,xsize=11,ysize=8,/inches,/color;,/landscape
;!P.thick=2
;!p.charthick=2
;!p.charsize=1.0
;!p.multi=0;[0,1,4]
;;start_wavelength=4600
;end_wavelength=7000
;
;
;
;plot,wavelength,disk_1D,/NODATA,yrange=[-0.1,1.9],$
;    xrange=[start_wavelength-100,end_wavelength+100],$
;    /xstyle,/ystyle,xthick=3,ythick=3,$;ytickinterval=30,$
;   ; ytickname=['Residuals','Galaxy + !CBest Fit','Disc','Bulge'],$
;    xtitle='Wavelength ('+cgSymbol("angstrom")+')',ytitle='Relative Flux';,title=galaxy_ref
;
;resid_smooth=fltarr(n_elements(resid_1D))
;for m=40,n_elements(resid_1D)-41,1 do resid_smooth[m]=median(resid_1D[m-40:m+40])
;resid_smooth[0:39]=resid_smooth[40]
;resid_smooth[-40:-1]=resid_smooth[-41]
;
;if n_comp ge 1100 then oplot,wavelength,(bulge_1D_orig/median(orig_1D)),color=cgcolor('red');/10000;+90
;oplot,wavelength,(disk_1D_orig/median(orig_1D)),color=cgcolor('blue');/10000;+60
;oplot,wavelength,(orig_1D/median(orig_1D));/10000;+30
;oplot,wavelength,((disk_1D_orig+bulge_1D_orig+resid_sky_1D)/median(orig_1D)),color=cgcolor('purple');/10000;+30,color=cgcolor('red')
;;oplot,wavelength,((bulge_1D+disk_1D)-median(bulge_1D+disk_1D))/10+10,color=cgcolor('red')
;oplot,wavelength,(resid_sky_1D/median(orig_1D)),color=cgcolor('dark grey');/10000,color=cgcolor('green')
;oplot,wavelength,(resid_1D/median(orig_1D)),color=cgcolor('olive');/10000,color=cgcolor('green')
;
;
;if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then oplot,wavelength,(comp3_1D/median(orig_1D)),color=cgcolor('skyblue')
;
;if n_comp eq 1000 or n_comp eq 1001 then legend,['Integrated spectrum from datacube','Bulge + Disc','Disc','Residual sky','Residuals'],linestyle=[0,0,0,0,0],$
;  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('blue'),cgcolor('dark grey'),cgcolor('olive')],charsize=0.9,box=0,/left,/top
;
;if n_comp eq 1100 or n_comp eq 1101 then legend,['Integrated spectrum from datacube','Galaxy + Nucleus','Galaxy','Nucleus','Residual sky','Residuals'],linestyle=[0,0,0,0,0,0],$
;  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('red'),cgcolor('blue'),cgcolor('dark grey'),cgcolor('olive')],charsize=0.8,box=0,/left,/top
;
;if n_comp eq 1010 or n_comp eq 1011 then legend,['Integrated spectrum from datacube','Centre + Disc','Centre','Disc','Residual sky','Residuals'],linestyle=[0,0,0,0,0,0],$
;  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('skyblue'),cgcolor('blue'),cgcolor('dark grey'),cgcolor('olive')],charsize=0.9,box=0,/left,/top
;
;if n_comp eq 1110 or n_comp eq 1111 then legend,['Integrated spectrum from datacube','Centre + Bulge + Disc','Bulge','Centre','Disc','Residual sky','Residuals'],linestyle=[0,0,0,0,0,0,0],$
;  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('red'),cgcolor('skyblue'),cgcolor('blue'),cgcolor('dark grey'),cgcolor('olive')],charsize=0.9,box=0,/left,/top
;
;;a=(Redshift+1)*(Redshift+1)
;;Velocity=3e5*((a-1)/(a+1))
;CaII_new=3934+(3934*Redshift)
;Ha_new=6563+(6563*Redshift)
;Hb_new=4861+(4861*Redshift)
;Hg_new=4341+(4341*Redshift)
;Hd_new=4102+(4102*Redshift)
;Mgb_new=5177+(5177*Redshift)
;Fe5335_new=5335+(5335*Redshift)
;Fe5270_new=5270+(5270*Redshift)
;Na_new= 5895.92+(5895.92*Redshift)
;CaII2_new= 8542.09+(8542.09*Redshift)
;Pae_new= 9546+(9546*Redshift)
;Paz_new= 9229+(9229*Redshift)
;
;
;;oplot,[Ha_new,Ha_new],[-5,5],linestyle=1
;;oplot,[Hb_new,Hb_new],[-5,5],linestyle=1
;;oplot,[Hg_new,Hg_new],[-5,5],linestyle=1
;;oplot,[Hd_new,Hd_new],[-5,5],linestyle=1
;;oplot,[Mgb_new,Mgb_new],[-5,5],linestyle=1
;;oplot,[Fe5335_new,Fe5335_new],[-5,5],linestyle=1
;;oplot,[Fe5270_new,Fe5270_new],[-5,5],linestyle=1
;
a=where(wavelength ge Ha_new)
b=where(wavelength ge Hb_new)
c=where(wavelength ge Hg_new)
d=where(wavelength ge Hd_new)
e=where(wavelength ge Mgb_new)
f=where(wavelength ge Fe5335_new)
g=where(wavelength ge Fe5270_new)
h=where(wavelength ge Na_new)
i=where(wavelength ge CaII_new)
j=where(wavelength ge CaII2_new)
k=where(wavelength ge Pae_new)
l=where(wavelength ge Paz_new)

if a[0]-30 gt 0 and a[0]+30 lt n_elements(orig_1D)-1 then $
  aa=mean((orig_1D[a[0]-30:a[0]+30]/median(orig_1D)))+0.2 $
  else aa=-1000
if b[0]-30 gt 0 and b[0]+30 lt n_elements(orig_1D)-1 then $
  bb=mean((orig_1D[b[0]-30:b[0]+30]/median(orig_1D)))+0.2 $
  else bb=-1000
if c[0]-30 gt 0 and c[0]+30 lt n_elements(orig_1D)-1 then $
  cc=mean((orig_1D[c[0]-30:c[0]+30]/median(orig_1D)))+0.2 $
  else cc=-1000
if d[0]-30 gt 0 and d[0]+30 lt n_elements(orig_1D)-1 then $
  dd=mean((orig_1D[d[0]-30:d[0]+30]/median(orig_1D)))+0.2 $
  else dd=-1000
if e[0]-30 gt 0 and e[0]+30 lt n_elements(orig_1D)-1 then $
  ee=mean((orig_1D[e[0]-30:e[0]+30]/median(orig_1D)))+0.2 $
  else ee=-1000
if f[0]-30 gt 0 and f[0]+30 lt n_elements(orig_1D)-1 then $
  ff=mean((orig_1D[f[0]-30:f[0]+30]/median(orig_1D)))+0.2 $
  else ff=-1000
if g[0]-30 gt 0 and g[0]+30 lt n_elements(orig_1D)-1 then $
  gg=mean((orig_1D[g[0]-30:g[0]+30]/median(orig_1D)))+0.2 $
  else gg=-1000
if h[0]-30 gt 0 and h[0]+30 lt n_elements(orig_1D)-1 then $
  hh=mean((orig_1D[h[0]-30:h[0]+30]/median(orig_1D)))+0.2 $
  else hh=-1000
if i[0]-30 gt 0 and i[0]+30 lt n_elements(orig_1D)-1 then $
  ii=mean((orig_1D[i[0]-30:i[0]+30]/median(orig_1D)))+0.2 $
  else ii=-1000
if j[0]-30 gt 0 and j[0]+30 lt n_elements(orig_1D)-1 then $
  jj=mean((orig_1D[j[0]-30:j[0]+30]/median(orig_1D)))+0.2 $
  else jj=-1000
if k[0]-30 gt 0 and k[0]+30 lt n_elements(orig_1D)-1 then $
  kk=mean((orig_1D[k[0]-30:k[0]+30]/median(orig_1D)))+0.2 $
  else kk=-1000
if l[0]-30 gt 0 and l[0]+30 lt n_elements(orig_1D)-1 then $
  ll=mean((orig_1D[l[0]-30:l[0]+30]/median(orig_1D)))+0.2 $
  else ll=-1000
;
;
;xyouts,Ha_new-30,aa,'H'+greek('alpha'),charsize=0.8
;xyouts,Hb_new-30,bb,'H'+greek('beta'),charsize=0.8
;xyouts,Hg_new-30,cc,'H'+greek('gamma'),charsize=0.8
;xyouts,Hd_new-30,dd,'H'+greek('delta'),charsize=0.8
;xyouts,Mgb_new-30,ee,'Mg',charsize=0.8
;xyouts,Fe5335_new-30,ff,'Fe',charsize=0.8
;xyouts,Fe5270_new-30,gg,'Fe',charsize=0.8
;xyouts,Na_new-30,hh,'Na D',charsize=0.8
;xyouts,CaII_new-30,0.95,'Ca II',charsize=0.8
;;xyouts,CaII_new-30,ii,'Ca II',charsize=0.9
;;xyouts,CaII2_new-30,jj,'Ca II',charsize=0.9
;
;;!p.multi=0
;device,/close
;stop
;device,file='/Users/ejohnsto/Dropbox/papers/Paper4/decomposed_spectra_1D_red_'+galaxy_ref+'.eps',xsize=19.5,ysize=11,/portrait;,/landscape
;;device,file=root+decomp+'decomposed_data/Spectra_integrated_3.eps',/landscape;,xsize=11,ysize=8,/inches,/color;,/landscape
;!P.thick=3
;!p.charthick=3
;!p.charsize=1.0
;!p.multi=0;[0,1,4]
;start_wavelength=6700
;end_wavelength=10000
;
;sample=where(wavelength le end_wavelength)
;
;plot,wavelength[sample],disk_1D[sample],/NODATA,yrange=[-0.1,1.9],$
;    xrange=[start_wavelength-100,end_wavelength+100],$
;    /xstyle,/ystyle,xthick=3,ythick=3,$;ytickinterval=30,$
;   ; ytickname=['Residuals','Galaxy + !CBest Fit','Disc','Bulge'],$
;    xtitle='Wavelength ('+cgSymbol("angstrom")+')',ytitle='Relative Flux';,title=galaxy_ref
;
;
;if n_comp ge 1100 then oplot,wavelength[sample],(bulge_1D_orig/median(orig_1D)),color=cgcolor('red');/10000;+90
;oplot,wavelength[sample],(disk_1D_orig/median(orig_1D)),color=cgcolor('blue');/10000;+60
;oplot,wavelength[sample],(orig_1D/median(orig_1D));/10000;+30
;oplot,wavelength[sample],((bulge_1D_orig+disk_1D_orig+resid_sky_1D)/median(orig_1D)),color=cgcolor('purple');/10000;+30,color=cgcolor('red')
;;oplot,wavelength,((bulge_1D+disk_1D)-median(bulge_1D+disk_1D))/10+10,color=cgcolor('red')
;oplot,wavelength[sample],(resid_sky_1D/median(orig_1D)),color=cgcolor('dark grey');/10000,color=cgcolor('green')
;oplot,wavelength[sample],(resid_1D/median(orig_1D)),color=cgcolor('olive');/10000,color=cgcolor('green')
;
;
;if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then oplot,wavelength[sample],(comp3_1D/median(orig_1D)),color=cgcolor('skyblue')
;
;if n_comp eq 1000 or n_comp eq 1001 then legend,['Integrated spectrum from datacube','Bulge + Disc','Disc','Residual sky','Residuals'],linestyle=[0,0,0,0,0],$
;  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('blue'),cgcolor('dark grey'),cgcolor('olive')],charsize=0.9,box=0,/left,/top
;
;if n_comp eq 1100 or n_comp eq 1101 then legend,['Integrated spectrum from datacube','Bulge + Disc','Bulge','Disc','Residual sky','Residuals'],linestyle=[0,0,0,0,0,0],$
;  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('red'),cgcolor('blue'),cgcolor('dark grey'),cgcolor('olive')],charsize=0.9,box=0,/left,/top
;
;if n_comp eq 1010 or n_comp eq 1011 then legend,['Integrated spectrum from datacube','Centre + Disc','Centre','Disc','Residual sky','Residuals'],linestyle=[0,0,0,0,0,0],$
;  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('skyblue'),cgcolor('blue'),cgcolor('dark grey'),cgcolor('olive')],charsize=0.9,box=0,/left,/top
;
;if n_comp eq 1110 or n_comp eq 1111 then legend,['Integrated spectrum from datacube','Centre + Bulge + Disc','Bulge','Centre','Disc','Residual sky','Residuals'],linestyle=[0,0,0,0,0,0,0],$
;  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('red'),cgcolor('skyblue'),cgcolor('blue'),cgcolor('dark grey'),cgcolor('olive')],charsize=0.9,box=0,/left,/top
;
;
;xyouts,Ha_new-30,aa,'H'+greek('alpha'),charsize=0.9
;;xyouts,Hb_new-30,bb,'H'+greek('beta'),charsize=0.9
;;xyouts,Hg_new-30,cc,'H'+greek('gamma'),charsize=0.9
;;xyouts,Hd_new-30,dd,'H'+greek('delta'),charsize=0.9
;;xyouts,Mgb_new-30,ee,'Mg',charsize=0.9
;;xyouts,Fe5335_new-30,ff,'Fe',charsize=0.9
;;xyouts,Fe5270_new-30,gg,'Fe',charsize=0.9
;;xyouts,Na_new-30,hh,'Na D',charsize=0.9
;;xyouts,CaII_new-30,ii,'Ca II',charsize=0.9
;xyouts,CaII2_new-30,jj,'Ca II',charsize=0.9
;if k[0] ne -1 then xyouts,Pae_new-30,kk,'Pa'+greek('epsilon'),charsize=0.9
;if l[0] ne -1 then xyouts,Paz_new-30,ll,'Pa'+greek('zeta'),charsize=0.9
;
;
;device,/close



device,file=root+decomp+decomp_dir+'decomposed_spectra_1D.eps',xsize=26,ysize=12,/portrait;,/landscape

;device,file='/Users/ejohnsto/Dropbox/papers/Paper4/decomposed_spectra_1D_full_'+galaxy_ref+'.eps',xsize=50,ysize=10,/portrait;,/landscape
;device,file=root+decomp+'decomposed_data/Spectra_integrated_3.eps',/landscape;,xsize=11,ysize=8,/inches,/color;,/landscape
!P.thick=4
!p.charthick=4
!p.charsize=1.0
!p.multi=0;[0,1,4]
start_wavelength=4700
end_wavelength=7220

sample=where(wavelength le end_wavelength)

plot,wavelength[sample],disk_1D[sample],/NODATA,yrange=[-0.1,1.9],$
    xrange=[start_wavelength,end_wavelength],$
    /xstyle,/ystyle,xthick=4,ythick=4,$;ytickinterval=30,$
   ; ytickname=['Residuals','Galaxy + !CBest Fit','Disc','Bulge'],$
    xtitle='Wavelength ('+cgSymbol("angstrom")+')',ytitle='Relative Flux';,title=galaxy_ref


if n_comp ge 1100 then oplot,wavelength[sample],(bulge_1D_orig/median(orig_1D)),color=cgcolor('red');/10000;+90
oplot,wavelength[sample],(disk_1D_orig/median(orig_1D)),color=cgcolor('blue');/10000;+60
oplot,wavelength[sample],(orig_1D/median(orig_1D));/10000;+30
oplot,wavelength[sample],((bulge_1D_orig+disk_1D_orig+resid_smooth)/median(orig_1D)),color=cgcolor('purple');/10000;+30,color=cgcolor('red')
;oplot,wavelength,((bulge_1D+disk_1D)-median(bulge_1D+disk_1D))/10+10,color=cgcolor('red')
oplot,wavelength[sample],(resid_sky_1D/median(orig_1D)),color=cgcolor('dark grey');/10000,color=cgcolor('green')
oplot,wavelength[sample],(resid_1D/median(orig_1D)),color=cgcolor('olive');/10000,color=cgcolor('green')


if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then oplot,wavelength[sample],(comp3_1D/median(orig_1D)),color=cgcolor('skyblue')

if n_comp eq 1000 or n_comp eq 1001 then legend,['Integrated spectrum from datacube','Bulge + Disc','Disc','Residual sky','Residuals'],linestyle=[0,0,0,0,0],$
  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('blue'),cgcolor('dark grey'),cgcolor('olive')],charsize=0.9,box=0,/left,/top

if n_comp eq 1100 or n_comp eq 1101 then legend,['Integrated spectrum from datacube','Bulge + Disc','Bulge','Disc','Residual sky','Residuals'],linestyle=[0,0,0,0,0,0],$
  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('red'),cgcolor('blue'),cgcolor('dark grey'),cgcolor('olive')],charsize=0.9,box=0,/left,/top

if n_comp eq 1010 or n_comp eq 1011 then legend,['Integrated spectrum from datacube','Centre + Disc','Centre','Disc','Residual sky','Residuals'],linestyle=[0,0,0,0,0,0],$
  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('skyblue'),cgcolor('blue'),cgcolor('dark grey'),cgcolor('olive')],charsize=0.9,box=0,/left,/top

if n_comp eq 1110 or n_comp eq 1111 then legend,['Integrated spectrum from datacube','Centre + Bulge + Disc','Bulge','Centre','Disc','Residual sky','Residuals'],linestyle=[0,0,0,0,0,0,0],$
  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('red'),cgcolor('skyblue'),cgcolor('blue'),cgcolor('dark grey'),cgcolor('olive')],charsize=0.9,box=0,/left,/top


xyouts,Ha_new-30,aa,'H'+greek('alpha'),charsize=0.9
xyouts,Hb_new-30,bb,'H'+greek('beta'),charsize=0.9
xyouts,Hg_new-30,cc,'H'+greek('gamma'),charsize=0.9
xyouts,Hd_new-30,dd,'H'+greek('delta'),charsize=0.9
xyouts,Mgb_new-30,ee,'Mg',charsize=0.9
xyouts,Fe5335_new-30,ff,'Fe',charsize=0.9
xyouts,Fe5270_new-30,gg,'Fe',charsize=0.9
xyouts,Na_new-30,hh,'Na D',charsize=0.9
xyouts,CaII_new-30,ii,'Ca II',charsize=0.9
xyouts,CaII2_new-30,jj,'Ca II',charsize=0.9
if k[0] ne -1 then xyouts,Pae_new-30,kk,'Pa'+greek('epsilon'),charsize=0.9
if l[0] ne -1 then xyouts,Paz_new-30,ll,'Pa'+greek('zeta'),charsize=0.9


device,/close

stop
end

