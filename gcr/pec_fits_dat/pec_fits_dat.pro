pro PEC_FITS_DAT

;--------------------------------------------------------------------------
; Program to calculate ionization balance, excited populations and
; spectrum.  The user specifies
; (1) The element name
; (2) The atomic mass number for the element
; (3) The atomic datafile
; (4) The wavelength range for the spectral calculation  
;--------------------------------------------------------------------------

;--------------------------------------------------------------------------
; Choose element and datafile [Can be changed by user]
;--------------------------------------------------------------------------

element     = 'Ar'
ion_stage   = '0+'
atomic_mass = 39.948d0
datafile    = 'adf04'

wave_min    = 6500.d0  ;; <--------------change these to select a
wave_max    = 9500.d0  ;; <--------------new wavelength range
ntemps      = 35
ndens       = 24
n_meta      = [1,2,4]


;--------------------------------------------------------------------------
; Define temperature (eV) and density arrays (cm^-3) [Can be changed by user]
;--------------------------------------------------------------------------

tmin         = 0.5    ;; <--------------change these to select a
tmax         = 10.5    ;; <--------------new temperature range

te_ev        = adas_vector(low=tmin, high=tmax, num=ntemps, /linear)

dens_min     = 1.e8   ;; <--------------change these to select a
dens_max     = 1.e11  ;; <--------------new density range

dens         = adas_vector(low=dens_min, high=dens_max, num=ndens)
;dens         = findgen(ndens, increment = (dens_max - dens_min) / ndens start = dens_min)
;dens= [1e9, 2e9, 3e9, 4e9, 5e9, 6e9, 7e9, 9e9, 1e10, 2e10, 3e10, 4e10, 5e10, 6e10, 7e10, 8e10, 9e10, 1e11]
;Use the resolution of the spectrometer
lambda_d=2.0d0
fwhm = 2.d0*(alog(2.d0))^0.5d0*lambda_d

;--------------------------------------------------------------------------
; Convert amu to kg and work out mc^2 in eV and choose 
; min A-value for spectral lines.
;--------------------------------------------------------------------------

atomic_mass = 1.66053886d-27*atomic_mass
mcsq        = atomic_mass*3.d8*3.d8/1.6d-19
amin        = 2.e-30

;--------------------------------------------------------------------------
; Calculate excited populations for Ar2+
;--------------------------------------------------------------------------

cd,'./',current=pwd
run_adas208,    adf04 = datafile        , adf18     = ' '     , $
     		te    = te_ev   	, dens      = dens    , $
		meta  = n_meta		, zeff      = 1       , $
     		ion   = 0	   	, cx        = 0       , $
		rec   = 0               , pec       = 1       , $
                sxb   = 0               , gcr       = 0	      , $
		pop   = pop	        , wmin      = wave_min, $
                wmax  = wave_max        , amin     = amin     , $
		log    = 'paper.txt'    , pass_dir = pwd
stop
;--------------------------------------------------------------------------
;--------------------------------------------------------------------------
; GENERATE A SPECTRUM 
;--------------------------------------------------------------------------
;--------------------------------------------------------------------------

;--------------------------------------------------------------------------
;Read the PEC file and get all of the spectral lines
;--------------------------------------------------------------------------
read_adf15,file='pec.pass',te=te_ev,dens=dens,$
   	fulldata=pec,/all

;--------------------------------------------------------------------------
;  Choose wavelength mesh to map out line profiles
;--------------------------------------------------------------------------
dwavel = min(fwhm)/20.d0
n_wavel=ulong((wave_max-wave_min)/dwavel)
lambda=dblarr(n_wavel)
flux=dblarr(ntemps,ndens,n_wavel)
flux_tst=dblarr(n_wavel)

;--------------------------------------------------------------------------
;  Doppler broaden each line profile and store results in the array flux
;--------------------------------------------------------------------------
size_meta = n_elements(n_meta)
size_wave = n_elements(pec.wave)/size_meta
wavelengths = pec.wave[0:size_wave-1]


;lambda[0] = wave_min
;for i=1L,n_wavel-1 do begin
;    lambda[i] = wave_min+dwavel*double(i)
;   for j=0,n_elements(pec.wave)-1 do begin
;       if ((lambda[i] gt (pec.wave[j]-fwhm*5.d0)) and (lambda[i] lt (pec.wave[j]+fwhm*5.d0)))  then begin
;	       for ii=0,ntemps-1 do begin
;                   for jj=0,ndens-1 do begin
;          flux[ii,jj,i] = flux[ii,jj,i] + pec.pec[ii,jj,j]*(2.d0*asin(1.d0))^(-0.5d0)/lambda_d*exp(-abs((lambda[i]-pec.wave[j])/lambda_d)^2.d0)
;	           endfor
;	       endfor
;       endif
;    endfor
;endfor

;--------------------------------------------------------------------------
; Create a reduced array, flux_reduced, that only contains non-zero
; flux points.
;--------------------------------------------------------------------------
index=where(flux[0,0,*] ne 0.0,n_w)
flux_reduced=flux[*,*,index]
lambda=lambda/10
lambda_red=lambda[index]
pec_wave=pec.wave[0:size_wave-1]/10
;--------------------------------------------------------------------------
; WRITE THE DATA TO A FITS FILE
; The /APPEND keyword adds the next set of data as a new HDU, and retains 
; info (shape, etc).
;--------------------------------------------------------------------------
;--------------------------------------------------------------------------
writefits, 'pec_dat.fits', te_ev
writefits, 'pec_dat.fits', n_meta, /APPEND
writefits, 'pec_dat.fits', dens, /APPEND
writefits, 'pec_dat.fits', pec_wave, /APPEND
writefits, 'pec_dat.fits', pec.pec, /APPEND

print, '' & print, 'Type .c to continue...' & print,''
stop

print,'' & print,'YOU HAVE SELECTED Te = ',te_ev,' eV (',ntemps,'temps)'
print,''
print, 'and Ne = ',dens,' cm^-3 (',ndens,' densities)' & print,''
print,'' & print,'NUMBER OF WAVELENGTH POINTS = ',n_wavel & print,''
print,'' & print,'END OF PROGRAM' & print,''

END


;Optional Code

;openw,fd_data,'lambda.dat',/get_lun
;printf,fd_data,lambda, format='(f10.3)'
;close,fd_data
;
;openw,fd_data,'lambda_red.dat',/get_lun
;printf,fd_data,lambda_red, format='(f10.3)'
;close,fd_data
;
;openw,fd_data,'pec_wave.dat',/get_lun
;printf,fd_data,pec.wave/10.0, format='(f10.3)'
;close,fd_data
;
;openw,fd_data,'pec_dens.dat',/get_lun
;printf,fd_data,dens, format='(f18.3)'
;close,fd_data
;
;openw,fd_data,'pec_temps.dat',/get_lun
;printf,fd_data,te_ev, format='(f10.3)'
;close,fd_data
;
;writefits,'flux_dlam_20.fits',flux
;writefits,'flux_red_dlam_20.fits',flux_reduced
;writefits,'pec.fits',pec.pec


;surface,pec.pec[*,*,33]/pec.pec[*,*,23],te_ev,dens,/ylog,charsize=4.,az=10,$
;ytitle='Density (cm^-3)',xtitle='Te (eV)',ztitle='Line ratio'
;plot,te_ev,pec.pec[*,0,33]/pec.pec[*,0,23],xr=[2,5],yr=[4,7]


;--------------------------------------------------------------------------
; The following loop generates data for the line ratio in third index of 
; pec.pec.... i.e  pec.pec[i,j,29]/pec.pec[i,j,19] is the ratio of the
; 29th to the 19th line.
;--------------------------------------------------------------------------
;openw,1,'pec_ratio.dat'
;for i=0,n_elements(te_ev)-1 do begin
;     for j=0,n_elements(dens)-1 do begin
;          printf,1,dens[j],te_ev[i],pec.pec[i,j,33]/pec.pec[i,j,23]
;    endfor
;    printf,1,' '
;endfor    	   
;stop

; Calculate FWHM and Doppler widths for each spectral line
;fwhm=dblarr(n_elements(pec.wave))
;lambda_d=dblarr(n_elements(pec.wave))
;for j=0,n_elements(pec.wave)-1 do begin
;    lambda_d[j] =(2.d0*te_ev/mcsq)^0.5d0*pec.wave[j]    
;    fwhm[j] = 2.d0*(alog(2.d0))^0.5d0*lambda_d[j]
;endfor


; Make the plot.  
;ymax=max(pec.pec[te_indx,dens_indx,*])*10.
;ymin=min(pec.pec[te_indx,dens_indx,*])/1.e4
;;title='Spectrum for '+element+' at Te = '+strtrim(string(te_ev),2)+' (eV) and Ne = '+strtrim(string(dens),2)+' cm!U-3!N'
;plot,[pec.wave[0],pec.wave[0]],[ymin,pec.pec[te_indx,dens_indx,0]],xtitle='Wavelength (A)',$
;	ytitle='PEC (Ph/s/cm!U3!N)',$
;	charsize=1.5,yr=[ymin,ymax],thick=3,$
;	title=title,xstyle=1,xr=[wave_min,wave_max],/ylog,ystyle=1
;for i=1,n_elements(pec.wave)-1 do oplot,[pec.wave[i],pec.wave[i]],$
;	[ymin,pec.pec[te_indx,dens_indx,i]],thick=3
;oplot,lambda,flux
;
;plot,lambda,flux,xtitle='Wavelength (A)',$
;	ytitle='PEC (Ph/s/cm!U3!N)',charsize=1.5
;oplot,[pec.wave[0],pec.wave[0]],[ymin,pec.pec[te_indx,dens_indx,0]]
;for i=1,n_elements(pec.wave)-1 do oplot,[pec.wave[i],pec.wave[i]],$
;	[ymin,pec.pec[te_indx,dens_indx,i]],thick=3

;set_plot,'ps'
;filenm='plot_spectrum.ps'
;device,filename=filenm,/landscape
;plot,[pec.wave[0],pec.wave[0]],[ymin,pec.pec[te_indx,dens_indx,0]],xtitle='Wavelength (A)',$
;	ytitle='PEC (Ph/s/cm!U3!N)',$
;	charsize=1.5,yr=[ymin,ymax],thick=3,$
;	title=title,xstyle=1,xr=[wave_min,wave_max],/ylog,ystyle=1
;for i=1,n_elements(pec.wave)-1 do oplot,[pec.wave[i],pec.wave[i]],$
;	[ymin,pec.pec[te_indx,dens_indx,i]],thick=3
;oplot,lambda,flux
;device,/close
;set_plot,'x'

;openw,2,'spec_adas.dat'
;for j=0,ndens -1 do begin
;	for i=0,n_elements(lambda)-1 do printf,2,lambda[i],(flux[k,j,i],k=0,ntemps),format='(31(1e12.4,2x))'
;		printf,2,' '
;close,2


;--------------------------------------------------------------------------
; Close windows and exit program 

;--------------------------------------------------------------------------
;wdelete,0
