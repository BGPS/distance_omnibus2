

restore,'/d1/BGPS/distance-omnibus/local/BGPS_V2_grs_spectra_r46.sav',/ver
restore,'/d1/BGPS/distance-omnibus/local/BGPS_V2_velocities.sav',/ver
s = read_mrt('/d1/BGPS/distance-omnibus/local/mrt/bolocatV210.mrt')

help,grs,/str

fn = FILE_SEARCH('/d1/BGPS/edmonton2013/v2.0*labelmask.fits',count=n)

;; Find the CNUM range for these fields


cn = lonarr(2,n)

FOR i=0,n-1 DO BEGIN
   
   lab = readfits(fn[i],lhd)
   lab = lab[where(lab NE 0)]
   
   cn[*,i] = minmax(lab)
   
ENDFOR

cnums = minmax(cn)
print,cnums


ind = where(grs.cnum GE cnums[0] AND grs.cnum LE cnums[1], n)
print,n

grs = grs[ind]
v   = v[ind]
s   = s[ind]

print,minmax(v.l)

i2 = where(v.vlsr GE -200, n2)
print,float(n2)/float(n)

v_std = findgen(1500)*0.2 + (-100.)

SKIPPLOT = 1b


hcop   = v.mol[0].tmb * (~v.mol[0].multiv)
thco   = dblarr(n)
covlsr = dblarr(n)
colw   = dblarr(n)
thcoon = dblarr(n)

FOR i=0,n-1 DO BEGIN
   
   sp   = grs[i].spectrum
   spon = grs[i].onspec
   
   IF ~skipplot THEN BEGIN
      
      cgPlot,v_std,sp,yr=[0,3],xr=[-30,130],/xst,$
             title='G'+string(grs[i].l,grs[i].b,format="(F0.3,F+0.3)")
      vline,/h,0.3,color='blue'
      ;; vline,/h,0.4,color='red'
      
      vline,-5,linestyle=1
      vline,135,linestyle=1
      
      vline,v[i].vlsr,color='green',thick=2,linestyle=3
      
   ENDIF
   
   ;; Find contiguous regions above 0.30K
   
   sigind = where(sp GE 0.30, ns)
   IF ns GT 1 THEN BEGIN        ; Only consider detections
      print,[ns,minmax(v_std[sigind]),max(v_std[sigind])-min(v_std[sigind])]
      
      step = sigind[1:*]-sigind[0:*]
      IF max(step) EQ 1 THEN BEGIN ; Single peak
         
         ;; print,'PLOT!!!'
         estimates = [max(sp),v_std[median(sigind)],2]
         yf = mpfitpeak(v_std,sp,A,NTERMS=3,$
                        estimates=estimates)
         covlsr[i] = A[1]
         colw[i]   = A[2] * 2.355
         thco[i]   = max(sp,jjj)
         thcoon[i] = spon[jjj]
         
         ;; print,estimates
         ;; print,A
         IF ~skipplot THEN BEGIN
            cgOplot,v_std,yf,color='orange',thick=2
            vline,A[1],color='orange',thick=2,linestyle=4
         ENDIF
      ENDIF
      
   ENDIF
   
   IF ~skipplot THEN wait,0.5
ENDFOR

;;===================================================================
;; CO TA vs. BGPS Flux Density
;; window,0,xsize=800,ysize=800
myps,'../plots/bgps_thco_fd.eps'
cgPlot,thco,s.flux,symsize=0.5,psym=16,xtit='!u13!nCO On-Off T!dA!u*!n  [K]',$
       ytit='BGPS Flux Density  [Jy]',/ylog,ytickformat='exponent10',$
       charsize=1.0,xr=[0.1,4.5],/xst ;,xr=[0,2],yr=[0,2]
;; one2one,color='red',linestyle=3
myps,/done

;;===================================================================
;; CO TA ON vs. BGPS Flux Density
;; window,0,xsize=800,ysize=800
myps,'../plots/bgps_thcoon_fd.eps'
cgPlot,thcoon,s.flux,symsize=0.5,psym=16,xtit='!u13!nCO T!dA!u*!n  [K]',$
       ytit='BGPS Flux Density  [Jy]',/ylog,ytickformat='exponent10',$
       charsize=1.0,xr=[0.1,5.5],/xst ;,xr=[0,2],yr=[0,2]
;; one2one,color='red',linestyle=3
myps,/done

;;===================================================================
;; HCOP vs. BGPS Flux Density
;; window,0,xsize=800,ysize=800
myps,'../plots/bgps_hcop_fd.eps'
cgPlot,hcop,s.flux,symsize=0.5,psym=16,xtit='HCO!u+!n T!dA!u*!n  [K]',$
       ytit='BGPS Flux Density  [Jy]',/ylog,ytickformat='exponent10',$
       charsize=1.0,xr=[0.1,5.5],/xst ;,xr=[0,2],yr=[0,2]
;; one2one,color='red',linestyle=3
myps,/done

;;===================================================================
;; Antenna Temperature
;; window,0,xsize=800,ysize=800
myps,'../plots/hcop_thco_ta.eps'
cgPlot,thco,hcop,symsize=0.5,psym=16,xtit='!u13!nCO On-Off T!dA!u*!n  [K]',$
       ytit='HCO!u+!n(3-2) T!dA!u*!n  [K]',yr=[-0.5,5.5],/yst,$
       charsize=1.0,xr=[-0.5,4.5] ;,xr=[0,2],yr=[0,2]
one2one,color='red',linestyle=3
myps,/done

;; window,1,xsize=800,ysize=800
myps,'../plots/hcop_thco_cdf.eps'
jj = where(thco NE 0)
thcojj = thco[jj]
hcopjj = hcop[jj]
i0 = where(hcopjj EQ 0, n0)   
i1 = where(hcopjj NE 0, n1)
thi0 = thcojj[i0]     
thi1 = thcojj[i1]
cdf1 = ccdf(thi1,/cdf)
cdf0 = ccdf(thi0,/cdf)
cgPlot,thi1,cdf1,xtit='!u13!nCO On-Off T!dA!u*!N  [K]',$
       ytit='Cumulative Distribution',charsize=1.0
cgOplot,thi0,cdf0,color='blue'                        
al_legend,color=['black','blue'],linestyle=0,['Y','N']+' HCO!u+!n',$
          /bottom,/right,box=0,linsize=0.5
myps,/done


;;===================================================================
;; Velocities
myps,'../plots/hcop_thco_vlsr.eps',xsize=10
multiplot_xm,[2,1],xgap=0.04,mpcharsize=1.0
;; cgPlot,covlsr,v.mol[0].vlsr,psym=16,symsize=0.4,xr=[-50,150],$
;;        xtit='!u13!nCO V!dLSR!n  [km s!u-1!n]',$
;;        ytit='HCO!u+!n V!dLSR!n  [km s!u-1!n]',charsize=1.0

ikp = where(covlsr NE 0 AND v.mol[0].vlsr NE 0)
cgPlot,covlsr[ikp],v[ikp].mol[0].vlsr,psym=16,symsize=0.4,xr=[-50,150],$
       xtit='!u13!nCO V!dLSR!n  [km s!u-1!n]',$
       ytit='HCO!u+!n V!dLSR!n  [km s!u-1!n]',charsize=1.0

multiplot,/doyaxis

cgsig = cgSymbol('sigma')

plothist,covlsr[ikp]-v[ikp].mol[0].vlsr,bin=0.2,xr=[-12,12],/xst,$
         xtit='!u13!nCO V!dLSR!n - HCO!u+!n V!dLSR!n  [km s!u-1!n]',$
         charsize=1.0,ytit='N per bin',xarr,yarr
yf = mpfitpeak(xarr,yarr,A,nterms=3)
cgOplot,xarr,yf,color='orange'
al_legend,/top,/right,[cgsig+' = '+string(A[2],format="(F0.2)")+' km/s'],$
          charsize=0.8

myps,/done,/mp

;;===================================================================
;; Linewidths
myps,'../plots/hcop_thco_lw.eps',xsize=10
multiplot_xm,[2,1],xgap=0.04,mpcharsize=1.0
;; cgPlot,covlsr,v.mol[0].vlsr,psym=16,symsize=0.4,xr=[-50,150],$
;;        xtit='!u13!nCO V!dLSR!n  [km s!u-1!n]',$
;;        ytit='HCO!u+!n V!dLSR!n  [km s!u-1!n]',charsize=1.0

ikp = where(colw NE 0 AND v.mol[0].lw NE 0)
cgPlot,colw[ikp],v[ikp].mol[0].lw,psym=16,symsize=0.4,$;xr=[-50,150],$
       xtit='!u13!nCO V!dLSR!n FWHM  [km s!u-1!n]',$
       ytit='HCO!u+!n V!dLSR!n FWHM  [km s!u-1!n]',charsize=1.0
one2one,color='red',linestyle=3
multiplot,/doyaxis

plothist,colw[ikp]-v[ikp].mol[0].lw,bin=0.2,$;xr=[-12,12],/xst,$
         xtit='!u13!nCO V!dLSR!n FWHM - HCO!u+!n V!dLSR!n FWHM  [km s!u-1!n]',$
         charsize=1.0,ytit='N per bin',xarr,yarr
yf = mpfitpeak(xarr,yarr,A,nterms=3)
cgOplot,xarr,yf,color='orange'
al_legend,/top,/right,[cgsig+' = '+string(A[2],format="(F0.2)")+' km/s'],$
          charsize=0.8

myps,/done,/mp






iii = where(covlsr NE 0, niii)
print,niii,float(niii)/n


END


lv = readfits('/d1/BGPS/distance-omnibus/local/13co_lv/13COGAL_LV_plot.fits',hd)

pr = [-4,0]

skip = 1b
IF ~skip THEN BEGIN
   window,0                     ;,xsize=1200,ysize=800
   
   
   plotmap,lv,hd,/log,xr=[max(v.l),min(v.l)],ct=41,/notsquare,axsel=7,$
           charsize=1.0,xtickformat="(F0.1)",yr=[-30,130],range=pr,$
           outimg=crlv,outhdr=crhd,lonarr=larr,latarr=barr,xc=xc,yc=yc
   print,pr
   
   cgOplot,v.l,v.vlsr,psym=16,symsize=0.7,color='black'
ENDIF ELSE $
   plotmap,lv,hd,/log,xr=[max(v.l),min(v.l)],ct=41,/notsquare,axsel=7,$
           charsize=1.0,xtickformat="(F0.1)",yr=[-30,130],range=pr,$
           outimg=crlv,outhdr=crhd,lonarr=larr,latarr=barr,xc=xc,yc=yc,/noplot



;;======================
;; NEXT: Bin the On-Off GRS spectra into 0.02 degree bins, and repeat.



synlv = crlv*0.
synlvon = crlv*0.
synlvbdr = crlv*0.
synhd = crhd

bin = 0.02
yarr = histogram(v.l,loc=xarr,bin=bin,min=min(larr),reverse_indices=revind)
nx = n_elements(xarr)

;; window,2
;; cgPlot,xarr+0.01,yarr,psym=10







;; Loop over bins
FOR i=0,nx-1 DO BEGIN
   
   IF revind[i] NE revind[i+1] THEN BEGIN
      inds = revind[revind[i] : revind[i+1]-1]
      
      ni = n_elements(inds)
      ;; print,ni,yarr[i],ni-yarr[i]
      
      sp = double(yc*0.)
      spon = double(yc*0.)
      spbdr = double(yc*0.)
      
      ;; Loop over constituents to this bin
      FOR j=0,ni-1 DO BEGIN
         
         spbdr += interpol(grs[inds[j]].bdrspec,v_std,yc)
         spon += interpol(grs[inds[j]].onspec,v_std,yc)
         sp += interpol(grs[inds[j]].spectrum,v_std,yc)
         
         
      ENDFOR                    ; End loop over constituents
      sp /= ni
      
      ;; help,synlv[nx-1-i,*]
      synlv[nx-1-i,*]    = sp
      synlvon[nx-1-i,*]  = spon
      synlvbdr[nx-1-i,*] = spbdr
      
   ENDIF
ENDFOR                          ; End loop over bins

undefine,spr
spr = [-3,2]
window,2;,xsize=600,ysize=500
plotmap,synlvon,synhd,ct=41,axsel=7,/notsquare,charsize=1.0,xtickformat="(F0.1)",$
        range=spr,/log
cgOplot,v.l,v.vlsr,psym=16,symsize=0.7,color='black'
print,spr

undefine,spr
spr = [-5,0]
;; window,3;,xsize=600,ysize=500

myps,'./lv_grsmatch.eps',ct=41
plotmap,synlv,synhd,ct=41,axsel=7,/notsquare,charsize=1.0,xtickformat="(F0.1)",$
        range=spr,/log
cgOplot,v.l,v.vlsr,psym=16,symsize=0.4,color='black'
print,spr
myps,/done

undefine,spr
spr = [-3,2]
window,4;,xsize=600,ysize=500
plotmap,synlvbdr,synhd,ct=41,axsel=7,/notsquare,charsize=1.0,xtickformat="(F0.1)",$
        range=spr,/log
cgOplot,v.l,v.vlsr,psym=16,symsize=0.7,color='black'
print,spr








END

