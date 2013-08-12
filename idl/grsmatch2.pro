pro grsmatch2
  data_dir = '/srv/astro/erosolo/bgps/scutum/edmonton2013/'
  restore,file=data_dir+'v2.0_ds2_l030_13pca_bolocat.sav'
;  grs_data = readfits(data_dir+'grs-30-cube.fits',hd_grs)
  grs_data = readfits(data_dir+'MOS_029_contincluded.Tb.fit',hd_grs)
  bgps_map = readfits(data_dir+'v2.0_ds2_l030_13pca_map20.fits',hd_bgps)
  obj = readfits(data_dir+'v2.0_ds2_l030_13pca_labelmask.fits',hd_bgps_mask)

  IF n_elements(r) EQ 0 THEN r = 2
  elt = shift(dist(2*r+1,2*r+1),r,r) LE r
  
  ;; Check to see if this BGPS source lies within the GRS coverage
;  IF bgps.glon LT 14 OR bgps.glon GT 60 THEN BEGIN
;     message,'BGPS source lies outside GRS coverage...',/inf
;     RETURN, null_spec
;  ENDIF


  ;; Read in BGPS map and labelmap
;  data = readfits(map_dir+catroots, bolo_hd, /SILENT)
  
  bgps = bolocat_struct[20]

  
  ;; Determine where BOLOCAM points land in GRS data
  label = bgps.cloudnum
  objmask = obj eq label
  
  ind = where(objmask, ct)
  sz = size(bgps_map)
  wts = bgps_map[ind]
  x_orig = ind mod sz[1]
  y_orig = ind / sz[1]
  
  border = (dilate(objmask, elt)-objmask)*(obj eq 0)
  ind2 = where(border,ctbdr)
  x_bdr = ind2 mod sz[1]
  y_bdr = ind2 / sz[1]
  extast,hd_bgps,bolo_astrom
  xy2ad, x_orig, y_orig, bolo_astrom, glon, glat
  xy2ad, x_bdr, y_bdr, bolo_astrom, glonbdr, glatbdr
  
  wtorder = sort(glon)
  glon = glon[wtorder]
  glat = glat[wtorder]
  wts = wts[wtorder]
  lout = strcompress(string(round(glon)), /rem)
;  filename = 'grs-'+lout+'-cube.fits'
;  fn = grs_dir+filename
  
  wt2order = sort(glonbdr)
  glonbdr = glonbdr[wt2order]
  glatbdr = glatbdr[wt2order]
  lout2 = strcompress(string(round(glonbdr)), /rem)
  
;  runspec = fltarr(!MW.NVBINS)
;  bdrspec = fltarr(!MW.NVBINS)
  runspec = fltarr(1000)
  bdrspec = fltarr(1000)
  v_std = findgen(1000)-500
  ;; Loop through all BGPS pixels assigned to this source
  FOR j = 0, n_elements(ind)-1 DO BEGIN
;     IF ~ file_test(fn[j]) THEN CONTINUE
     
     ;; Pick GRS filename and load a new one if necessary
;     lastfn = fn[j]
     extast, hd_grs, astrom
     rdhd, hd_grs, s = h
     sz_grs = size(grs_data)
     ad2xy, glon[j], glat[j], astrom, x, y

     ;; Look up spectrum
     if x lt 0 or x gt sz_grs[1]-1 or $
        y lt 0 or y gt sz_grs[2]-1 then begin
        wts[j] = 0
        continue
     endif
     spectrum = interpolate(grs_data,replicate(x,sz_grs[3]),$
                            replicate(y,sz_grs[3]),indgen(sz_grs[3]))
     
     ;; Interpolate to common velocity scale
     specout = interpol(spectrum, h.v, v_std)
     specout[where(v_std lt min(h.v) or v_std gt max(h.v))] = !values.f_nan
     runspec = runspec+specout*wts[j]

  ENDFOR   ;; End of loop for ON-SOURCE spectra extraction
  
  
  ;; Loop through all BORDER pixels
  FOR j = 0,n_elements(ind2)-1 DO BEGIN
;     IF ~ file_test(fn2[j]) THEN CONTINUE
     
     ;; Pick GRS filename and load a new one if necessary
     ;; IF fn2[j] NE lastfn THEN BEGIN
     ;;    UNDEFINE,grs_data
     ;;    print,fn2[j]
     ;;    grs_data = readfits(fn2[j], hd)
     ;; ENDIF
;     lastfn = fn2[j]
     extast, hd_grs, astrom
     rdhd, hd_grs, s = h
     sz_grs = size(grs_data)
     ad2xy, glonbdr[j], glatbdr[j], astrom, x, y
     
     ;; Look up spectrum
     if x lt 0 or x gt sz_grs[1]-1 or $
        y lt 0 or y gt sz_grs[2]-1 then begin
        ctbdr = ctbdr - 1
        continue
     endif
     spectrum = interpolate(grs_data,replicate(x,sz_grs[3]),$
                            replicate(y,sz_grs[3]),indgen(sz_grs[3]))
     
     ;; Interpolate to common velocity scale
     specout = interpol(spectrum, h.v, v_std)
     specout[where(v_std lt min(h.v) or v_std gt max(h.v))] = !values.f_nan
     bdrspec = bdrspec+specout
     
  ENDFOR   ;; End of loop for OFF-SOURCE spectra extraction
  
  ;; Filter the spectra
  filter = savgol(17,17,0,6)
  
  specout = runspec/total(wts)
  spectrum = convol(specout,filter,/edge_trun,/nan) >0
  
  bdrspec = bdrspec/ctbdr
  bdrspec = convol(bdrspec,filter,/edge_trun,/nan) >0

  
  ;; Create the Final Spectrum     
  onspec = spectrum
  emmask = spectrum GT 4*mad(spectrum)
  spectrum = ((spectrum - bdrspec)>0)*emmask
  
  
  ;; Plot, if desired
  IF KEYWORD_SET( plot ) THEN BEGIN
     window,1
     plot,v_std,onspec,/xst,xtitle='V!DLSR!N (km/s)',$
          ytitle='!6T!DA!N (K)',linestyle=2
     oplot,v_std,bdrspec,ps=10,color='0000ff'x,linestyle=2
     oplot,v_std,spectrum,ps=10,color='ffff00'x

     v_tan = V0*(1. - sin(s.glon_peak * !dtor))
     plots,v_tan,!y.crange[0]
     plots,v_tan,!y.crange[1],/continue,linestyle=2
     print,!x.crange
     print,n_e(v_std),min(v_std),max(v_std)

  ENDIF
  
  ;; Square the spectrum  -- WHY?
;  spectrum = spectrum^2

;  IF total(spectrum) EQ 0 THEN RETURN,null_spec
  
;  RETURN,spectrum

  
  stop

  return
end
