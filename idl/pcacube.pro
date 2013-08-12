pro pcacube, data_in, size = size, linewidth = linewidth

  data_in = double(data_in)

  data = data_in;-total(data_in)/n_elements(data_in) 
  sz = size(data)
  bad_dataz = where(data ne data,ct) 
  if ct gt 0 then data[bad_dataz]=0.0
  matrix = reform(data, sz[1]*sz[2], sz[3])
  covar = transpose(matrix)#matrix/(sz[1]*sz[2])

  eigenval = eigenql(covar, eigenvectors = evec)
  eval_norm = eigenval/total(eigenval)

  proj = matrix#evec
  ims = reform(proj, sz[1], sz[2], sz[3])

  size = fltarr(sz[3])+!values.f_nan
  linewidth = fltarr(sz[3])+!values.f_nan

  noise_level = mad(data)
  noiseim = fltarr(sz[1], sz[2])
  ctr = 0
  for k = (3*sz[3]/4), sz[3]-1 do begin
    noiseim = noiseim+convolve(ims[*, *, k], /auto)/sz[1]/sz[2]/eigenval[k]
    ctr = ctr+1
  endfor
  noiseim = noiseim/ctr*noise_level^2

  maj = size
  min =size

  for z = 0, 3*sz[3]/4-1  do begin
    acor_image = convolve(ims[*, *, z], /auto)
    acor_image = acor_image/total(acor_image^2)*total(ims[*, *, z]^2)
    acor_image = acor_image-noiseim
    null = max(noiseim, maxind)
    null = max(acor_image)
    x0 = maxind mod sz[1]
    y0 = maxind / sz[1]
    xvec = acor_image[sz[1]/2:*, sz[2]/2]/NULL
    yvec = reform(acor_image[sz[1]/2, sz[2]/2:*]/NULL)
    
    l = label_region(xvec ge 0)
    maxind = max(where(l eq 1))+1
    maj[z] = interpol(findgen(maxind+1), $
                         xvec[findgen(maxind+1)], exp(-1))
    l = label_region(yvec ge 0)
    maxind = max(where(l eq 1))+1
    min[z] = interpol(findgen(maxind+1), $
                         yvec[findgen(maxind+1)], exp(-1))


    vec = evec[*, z]
    acor_vec = a_correlate(vec, findgen(sz[3]))
    l = label_region(acor_vec ge 0)
    maxind = max(where(l eq 1))+1
    linewidth[z] = interpol(findgen(maxind+1), $
                         acor_vec[findgen(maxind+1)], exp(-1))


    vec = evec[*, z]
    acor_vec = a_correlate(vec, findgen(sz[3]))
    l = label_region(acor_vec ge 0)
    maxind = max(where(l eq 1))+1
    linewidth[z] = interpol(findgen(maxind+1), $
                         acor_vec[findgen(maxind+1)], exp(-1))
  endfor

  size = sqrt(maj^2+min^2)


  return
end
