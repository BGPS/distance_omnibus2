

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



;; cnums = minmax(cn)
cnums = minmax(s.cnum)
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

hco_one = ~v.mol[0].multiv

hcop   = v.mol[0].tmb * hco_one ; HCO+ Tmb
thco   = dblarr(n)              ; 13CO On-Off Tmb
covlsr = dblarr(n)              ; 13CO vlsr
colw   = dblarr(n)              ; 13CO linewidth
thcoon = dblarr(n)              ; 13CO On Tmb
mvrat  = dblarr(n)              ; Ratio of peaks for multi-velocity 13CO
nsegs  = bytarr(n)              ; Number of spectrum segments

FOR i=0,n-1 DO BEGIN            ; Loop through the N sources
   
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
      ;; print,[ns,minmax(v_std[sigind]),max(v_std[sigind])-min(v_std[sigind])]
      
      step = sigind[1:*]-sigind[0:*]
      ;; help,step,sigind
      
      
      ;;===================================================
      IF max(step) EQ 1 THEN BEGIN ; Single peak
         
         nseg = 1
         mvrat[i] = 0.

         ;;GOTO, skipsingle
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
         
         skipsingle:
      ENDIF ELSE BEGIN          ; Multipeak
         sigind = [sigind,n_elements(sp)]
         step = sigind[1:*]-sigind[0:*]
         
         
         ist = [where(step NE 1),n_elements(step)-1]
         ist = ist[uniq(ist,sort(ist))] ; Make unique
         
         nseg = n_elements(ist)
         print,'N Segments: ',nseg
         startind = 0
         ;; print,nseg,ist
         ;; print,sigind
         
         tpk = 0.
         
         tastar = fltarr(nseg)
         
         FOR jjj=0,nseg-1 DO BEGIN
            
            sigind2 = sigind[startind:ist[jjj]]
            
            estimates = [max(sp),v_std[median(sigind2)],2]
            yf = mpfitpeak(v_std,sp,A,NTERMS=3,$
                           estimates=estimates)
            
            ;;print,A[1],A[0]
            
            tastar[jjj] = max(sp[sigind2],jjjj)
            
            IF tastar[jjj] LT tpk THEN CONTINUE
            covlsr[i] = A[1]
            colw[i]   = A[2] * 2.355
            thco[i]   = max(sp[sigind2],jjjj)
            thcoon[i] = spon[jjjj]
            
            tpk = max(sp[sigind2],jjjj)
            
            ;; print,'Segment '+strtrim(jjj+1,2),$
            ;;       n_elements(step[startind:ist[jjj]])
            
            ;; print,step[startind:ist[jjj]]
            ;; print,sigind[startind:ist[jjj]]
            
            
            startind=ist[jjj]+1
         ENDFOR
         
         ;; FIND THE SMALLEST TPK/TASTAR RATIO GREATER THAN 1!
         
         tastar = tpk / tastar
         ind = where(tastar GT 1., nother)
         print,'Smallest ratio: ',min(tastar[ind])
         
         mvrat[i] = min(tastar[ind])
         
      ENDELSE
      
      nsegs[i] = nseg
   ENDIF
   
   IF ~skipplot THEN wait,0.5
ENDFOR
delta = thco*(mvrat - 1.d) / mvrat
keepmv = (mvrat EQ 0) OR (mvrat GE 1.5)

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
myps,'../plots/hcop_thco_vlsr.eps',xsize=10,ysize=4.5
multiplot_xm,[2,1],xgap=0.04,mpcharsize=1.0
;; cgPlot,covlsr,v.mol[0].vlsr,psym=16,symsize=0.4,xr=[-50,150],$
;;        xtit='!u13!nCO V!dLSR!n  [km s!u-1!n]',$
;;        ytit='HCO!u+!n V!dLSR!n  [km s!u-1!n]',charsize=1.0

ikp = where(covlsr NE 0 AND v.mol[0].vlsr NE 0 AND hco_one AND keepmv, nikp)
cgPlot,covlsr[ikp],v[ikp].mol[0].vlsr,psym=16,symsize=0.4,xr=[-50,150],$
       xtit='!u13!nCO V!dLSR!n  [km s!u-1!n]',$
       ytit='HCO!u+!n V!dLSR!n  [km s!u-1!n]',charsize=1.0
vline,/h,-5

iki = where( abs(covlsr[ikp]-v[ikp].mol[0].vlsr) GE 5, niki)
al_legend,/top,/left,box=0,charsize=0.8,$
          ['Fraction w/ > 5 km/s diff: '+$
           string(float(niki)/nikp,format="(F0.3)")]

iik = where(v[ikp].mol[0].vlsr LE -5)
cgOplot,covlsr[ikp[iik]],v[ikp[iik]].mol[0].vlsr,psym=16,symsize=0.4,$
        color='deep pink'

cgPlots,[-10,150],[-50,110],linestyle=3,color='blk4'
cgPlots,[-50,110],[-10,150],linestyle=3,color='blk4'

multiplot,/doyaxis

cgsig = cgSymbol('sigma')

plothist,covlsr[ikp]-v[ikp].mol[0].vlsr,bin=0.2,xr=[-12,12],/xst,$
         xtit='!u13!nCO V!dLSR!n - HCO!u+!n V!dLSR!n  [km s!u-1!n]',$
         charsize=1.0,ytit='N per bin',xarr,yarr
yf = mpfitpeak(xarr,yarr,A,nterms=3)
cgOplot,xarr,yf,color='cyan'
cgAxis,xaxis=0,/xst,xtickformat='blank_axis'
al_legend,/top,/right,[cgsig+' = '+string(A[2],format="(F0.2)")+' km/s'],$
          charsize=0.8,linestyle=0,color='cyan',linsize=0.5

myps,/done,/mp

;;===================================================================
;; Linewidths
myps,'../plots/hcop_thco_lw.eps',xsize=10
multiplot_xm,[2,1],xgap=0.04,mpcharsize=1.0
;; cgPlot,covlsr,v.mol[0].vlsr,psym=16,symsize=0.4,xr=[-50,150],$
;;        xtit='!u13!nCO V!dLSR!n  [km s!u-1!n]',$
;;        ytit='HCO!u+!n V!dLSR!n  [km s!u-1!n]',charsize=1.0

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


;;===================================================================
;; 13CO multi-vel peak ratio versus 13CO - HCO+ velocity diff
myps,'../plots/vdiff_vs_mvrat.eps',xsize=10,ysize=10
multiplot_xm,[2,2],gap=0.035,mpcharsize=1.0,/doyaxis,/doxaxis

ikp2 = where(v.mol[0].vlsr NE 0 AND mvrat GT 0.5 AND hco_one, nikp2)

cgPlot,abs(covlsr[ikp2]-v[ikp2].mol[0].vlsr),mvrat[ikp2],psym=16,symsize=0.4,$
       xtit='|!u13!nCO V!dLSR!n - HCO!u+!n V!dLSR!n|  [km s!u-1!n]',$
       charsize=1.0,xr=[-5,85],/xst,$
       ytit='!u13!nCO Multiple-Velocity Peak Ratio',yr=[0.5,10],/yst

vline,/h,2,color='cyan'
vline,/h,1.5,color='lime green'
vline,/h,1.25,color='crimson'

vline,40,color='blk4',linestyle=3

iik = where(v[ikp2].mol[0].vlsr LT -5, niik)
IF niik NE 0 THEN cgOplot,abs(covlsr[ikp2[iik]]-v[ikp2[iik]].mol[0].vlsr),$
                          mvrat[ikp2[iik]],psym=16,symsize=0.4,color='deep pink'

multiplot,/doyaxis,/doxaxis

plothist,mvrat,bin=0.2,xr=!y.crange,/xst,yst=8,charsize=1.0,$
         ytit='N',xtit='!u13!nCO Multiple-Velocity Peak Ratio'
vline,1

cgaxis,yaxis=1,yr=[0,1],/save,ytit='CDF',charsize=1.0,color='blk5'
mvi = where(mvrat GT 0.5)
mvcdf = mvrat[mvi]
cdf = ccdf(mvcdf,/cdf)
cgOplot,mvcdf,cdf,color='blk5'
vline,2,color='cyan'
vline,1.5,color='lime green'
vline,1.25,color='crimson'
vline,/h,interpol(cdf,mvcdf,2),color='cyan'
vline,/h,interpol(cdf,mvcdf,1.5),color='lime green'
vline,/h,interpol(cdf,mvcdf,1.25),color='crimson'

multiplot,/doyaxis,/doxaxis
;; DELTA

;; ikp2 = where(v.mol[0].vlsr NE 0 AND mvrat GT 0.5 AND hco_one, nikp2)

cgPlot,abs(covlsr[ikp2]-v[ikp2].mol[0].vlsr),delta[ikp2],psym=16,symsize=0.4,$
       xtit='|!u13!nCO V!dLSR!n - HCO!u+!n V!dLSR!n|  [km s!u-1!n]',$
       charsize=1.0,xr=[-5,85],/xst,$
       ytit='!u13!nCO Multiple-Velocity Peak Delta  [K]',yr=[-0.3,3.5],/yst

vline,/h,0.3,color='cyan'
;; vline,/h,1.5,color='cyan'
;; vline,/h,1.25,color='cyan'

vline,40,color='blk4',linestyle=3

iik = where(v[ikp2].mol[0].vlsr LT -5, niik)
IF niik NE 0 THEN cgOplot,abs(covlsr[ikp2[iik]]-v[ikp2[iik]].mol[0].vlsr),$
                          delta[ikp2[iik]],psym=16,symsize=0.4,color='deep pink'

multiplot,/doyaxis,/doxaxis

plothist,delta[where(finite(delta))],bin=0.05,xr=!y.crange,/xst,ytit='N',yst=8,$
         charsize=1.0,xtit='!u13!nCO Multiple-Velocity Peak Delta  [K]'


cgaxis,yaxis=1,yr=[0,1],/save,ytit='CDF',charsize=1.0,color='blk5'
mvi = where(mvrat GT 0.5)
mvcdf = delta[mvi]
cdf = ccdf(mvcdf,/cdf)
cgOplot,mvcdf,cdf,color='blk5'
vline,0.3,color='cyan'
;; vline,1.5,color='lime green'
;; vline,1.25,color='crimson'
vline,/h,interpol(cdf,mvcdf,0.3),color='cyan'
;; vline,/h,interpol(cdf,mvcdf,1.5),color='lime green'
;; vline,/h,interpol(cdf,mvcdf,1.25),color='crimson'

myps,/done,/mp


;;===================================================================
;; 13CO multi-vel peak ratio versus 13CO On-Off Tmb
myps,'../plots/mvrat_vs_thco_tmb.eps'

ikp3 = where( (v.mol[0].vlsr EQ 0 OR ~hco_one) AND mvrat GT 0.5,  nikp3)
print,'HCO+ Flag = 0: ',nikp3


cgPlot,thco[ikp2],mvrat[ikp2],psym=16,symsize=0.4,charsize=1.0,$
       xtit='!u13!nCO On-Off T!dA!u*!n  [K]',/nodata,yr=[0.5,10],/yst,$
       ytit='!u13!nCO Multiple-Velocity Peak Ratio'

xp = findgen(101)/100.*(!x.crange[1]-!x.crange[0])+!x.crange[0]
ikk = where(xp GT 0.3)
cgOplot,xp[ikk],1.d/(1.d - 0.3d/xp[ikk]),color='lime green',linestyle=4

vline,/h,color='blk5',linestyle=4,1.5
cgPlots,[1,!y.crange[1]]*0.3,[1,!y.crange[1]],color='crimson',linestyle=3
cgPlots,[0.3,!x.crange[1]],[1,1],color='crimson',linestyle=3

cgOplot,thco[ikp3],mvrat[ikp3],psym=16,symsize=0.4,color='black'
cgOplot,thco[ikp2],mvrat[ikp2],psym=16,symsize=0.4,color='cyan'

al_legend,/top,/left,psym=16,symsize=0.7,color=['black','cyan'],box=0,$
          ['No HCO!u+!n Velocity','Yes HCO!u+!n Velocity'],/clear

cgText,2.5,0.70,'Minimally Distinct',charsize=0.8,color='crimson',align=0.5
cgText,1.8,6.25,'Maximally Distinct',charsize=0.8,color='crimson',$
       orien=44,alig=0.5
cgText,0.4,6,cgSymbol('Delta')+'T = 0.3 K',charsize=0.8,align=0.5,$
       color='lime green',orien=-87


myps,/done


;;===================================================================
;; 13CO multi-vel peak DELTA versus 13CO On-Off Tmb
myps,'../plots/delta_vs_thco_tmb.eps'

ikp3 = where( (v.mol[0].vlsr EQ 0 OR ~hco_one) AND mvrat GT 0.5,  nikp3)
print,'HCO+ Flag = 0: ',nikp3

cgPlot,thco[ikp2],delta[ikp2],psym=16,symsize=0.4,charsize=1.0,$
       xtit='!u13!nCO On-Off T!dA!u*!n  [K]',/nodata,yr=[-0.3,3.5],/yst,$
       ytit='!u13!nCO Multiple-Velocity Peak Delta  [K]'

vline,/h,0.3,color='lime green',linestyle=4
xp = findgen(101)/100.*(!x.crange[1]-!x.crange[0])+!x.crange[0]
cgOplot,xp,xp/3.,color='blk5',linestyle=4

cgPlots,[0.3,!x.crange[1]],[0,0],color='crimson',linestyle=3
cgPlots,[0,!y.crange[1]]+0.3,[0,!y.crange[1]],color='crimson',linestyle=3

cgOplot,thco[ikp3],delta[ikp3],psym=16,symsize=0.4,color='black'
cgOplot,thco[ikp2],delta[ikp2],psym=16,symsize=0.4,color='cyan'

al_legend,/top,/left,psym=16,symsize=0.7,color=['black','cyan'],box=0,$
          ['No HCO!u+!n Velocity','Yes HCO!u+!n Velocity']

cgText,2.5,-0.20,'Minimally Distinct',charsize=0.8,color='crimson',align=0.5
cgText,1.8,1.61,'Maximally Distinct',charsize=0.8,color='crimson',$
       orien=37,alig=0.5
cgText,3,1.05,'MVRAT = 1.5',charsize=0.8,color='blk5',align=0.5,orien=15


myps,/done


iii = where(covlsr NE 0, niii)
print,niii,float(niii)/n







END
;;===========================================================================
;;===========================================================================
;;===========================================================================
;;===========================================================================
;;===========================================================================
;;===========================================================================






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

