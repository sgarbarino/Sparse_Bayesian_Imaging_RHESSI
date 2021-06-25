
;+
; NAME:
;   hsi_vis_bayes_compute_peaks
;
; PURPOSE:
;   Compute the max of the distribution
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;
;-

function hsi_vis_bayes_compute_peaks, Nsamp, H

  peaks = make_array(3,Nsamp^2)
  npeaks = -1

  for j = 0,Nsamp-1 do begin
    for i = 0,Nsamp-1 do begin

      n = 0
      s = 0
      w = 0
      e = 0
      n_w = 0
      s_w = 0
      n_e = 0
      s_e = 0

      c = H[j,i]
      if (i gt 0) then begin
        n = H[j,i-1]
      endif

      if (i lt Nsamp-1) then begin
        s = H[j,i+1]
      endif

      if (j gt 0) then begin
        w = H[j-1,i]
      endif

      if (j lt Nsamp-1) then begin
        e = H[j+1,i]
      endif

      if (j gt 0) and (i gt 0) then begin
        n_w = H[j-1,i-1]
      endif

      if (j gt 0) and (i lt Nsamp-1) then begin
        s_w = H[j-1,i+1]
      endif

      if (j lt Nsamp-1) and (i gt 0) then begin
        n_e = H[j+1,i-1]
      endif

      if (j lt Nsamp-1) and (i lt Nsamp-1) then begin
        s_e = H[j+1,i+1]
      endif

      pos = [n,s,w,e,n_w,s_w,n_e,s_e]
      neighb = 0

      for k = 0, N_elements(pos)-1 do begin
        if (c ge pos[k]) then begin
          neighb = neighb + 1
        endif
      endfor

      if (neighb eq 8) then begin
        npeaks = npeaks+1
        peaks[0,npeaks] = j
        peaks[1,npeaks] = i
        peaks[2,npeaks] = c
      endif

    endfor
  endfor

  peaks = peaks[*,reverse(sort(peaks[2,*]))]
  print, peaks
  return, peaks

end
;----------------------------------------------------------------------------


;+
; NAME:
;   hsi_vis_bayes_data_stack
;
; PURPOSE:
;   Data stacking
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;
;-

function hsi_vis_bayes_data_stack, param,  num_sources,  sample, weights, H, types, weight_type=weight_type, pc_pe_pl=pc_pe_pl

  ; Data stacking
  mean_param = make_array(param.Nsamp,param.Nsamp,param.N_param)
  weight_tmp = make_array(param.Nsamp,param.Nsamp)
  weight_type = make_array(param.Nsamp,param.Nsamp,3)
  pc_pe_pl = make_array(param.Nsamp, param.Nsamp, 3)

  par_tmp=[0,1,4,5]

  for p = 0, param.N_particles-1 do begin
    for k = 0, num_sources[p]-1 do begin

      r = floor(param.Nsamp*(sample[k,1,p]-param.ysamp[0])/(param.ysamp[param.Nsamp-1]-param.ysamp[0])) ; find the bin where the k-th source is located
      s = floor(param.Nsamp*(sample[k,0,p]-param.xsamp[0])/(param.xsamp[param.Nsamp-1]-param.xsamp[0]))

      if (r ge 0) and (r lt param.Nsamp) and (s ge 0) and (s lt param.Nsamp) then begin ; if source is inside the fov

        mean_param[s,r,par_tmp] = mean_param[s,r,par_tmp] + sample[k,par_tmp,p]*weights[p]

        ; Compute the mean value (according to the shape) and the number of circles, ellipses, and loops
        if (sample[k,3,p] lt param.ecc_toll_min) then begin
          weight_type[s,r, 0] =  weight_type[s,r, 0] + weights[p]
          pc_pe_pl[s,r,0]=pc_pe_pl[s,r,0]+1
        endif else if ((sample[k,3,p] ge param.ecc_toll_min) and (sample[k,3,p] lt param.ecc_toll_max)) then begin
          if (abs(sample[k,6,p]) lt param.loop_toll_max) then begin
            weight_type[s,r, 1] =  weight_type[s,r, 1] + weights[p]
            pc_pe_pl[s,r,1]=pc_pe_pl[s,r,1]+1
          endif else begin
            weight_type[s,r, 2] =  weight_type[s,r, 2] + weights[p]
            pc_pe_pl[s,r,2]=pc_pe_pl[s,r,2]+1
          endelse
        endif else if (sample[k,3,p] ge param.ecc_toll_max) then begin
          if (abs(sample[k,6,p]) ge param.loop_toll_min) then begin
            weight_type[s,r, 2] =  weight_type[s,r, 2] + weights[p]
            pc_pe_pl[s,r,2]=pc_pe_pl[s,r,2]+1
          endif else begin
            weight_type[s,r, 1] =  weight_type[s,r, 1] + weights[p]
            pc_pe_pl[s,r,1]=pc_pe_pl[s,r,1]+1
          endelse
        endif

        ; For the pa, eccentricity and loop angle we do not use the mean value but the mode
        if (weights[p] ge weight_tmp[s,r]) then begin
          mean_param[s,r, [2,3,param.N_param-1]] = sample[k,[2,3,param.N_param-1],p]
          weight_tmp[s,r]=weights[p]
        endif

      endif

    endfor
  endfor


  for par = 0, 3 do begin
    mean_param[*,*,par_tmp[par]] = mean_param[*,*,par_tmp[par]]/H
  endfor

  return, mean_param
end
;----------------------------------------------------------------------------



;+
; NAME:
;   hsi_vis_bayes_reconstr_vis
;
; PURPOSE:
;  Find the solution and show the map of the event
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;
;-

function hsi_vis_bayes_reconstr_vis, param, distr_num_src,  H, jdum, mean_param,  shape,  visxyobs, visxyrec=visxyrec, sample, N_src_expected=N_src_expected, weight_type, pc_pe_pl, prop_cel_100=prop_cel_100, time_title
  COMMON uvdata, u,v, pa, mapcenter

  ; Determine the number of sources of the recovered image
  N_src_expected = where(distr_num_src eq max(distr_num_src))
  N_src_expected = N_src_expected[0]

  if (N_src_expected ge 1) then begin

    smcsources = replicate({hsi_vis_src_structure},N_src_expected) ; initialization of the optimal particle

    peaks=hsi_vis_bayes_compute_peaks(param.Nsamp, H) ; find the local maximum of the distribution

    ; Initialization
    visxyrec = dblarr(N_elements(jdum))
    prob_cel=make_array(N_src_expected,3)
    prop_cel_100=prob_cel

    for sorg = 0, N_src_expected-1 do begin
      pa=0
      ecc=0
      la=0

      ; Determine the index where there is the maximum of the distribution
      j = peaks[0,sorg]
      i = peaks[1,sorg]

      ; According to the shape, determine the eccentricity, pa, and loop angle
      if (weight_type[j,i,1] ge weight_type[j,i,0]) and (weight_type[j,i,1] ge weight_type[j,i,2]) then begin
        pa = mean_param[j,i,2]
        ecc = mean_param[j,i,3]
      endif else if (weight_type[j,i,2] ge weight_type[j,i,0]) and (weight_type[j,i,2] ge weight_type[j,i,1]) then begin
        pa = mean_param[j,i,2]
        ecc = mean_param[j,i,3]
        la = mean_param[j,i,6]
      endif

      ; A posteriori evaluation of the shape of the source (it could be deleted since it should be not necessary anymore)
      if (ecc lt param.ecc_toll_min) then begin
        type = 0
      endif $
      else if ((ecc ge param.ecc_toll_min) and (ecc lt param.ecc_toll_max)) then begin
        if (abs(la) lt param.loop_toll_max) then begin
          type = 1
        endif else begin
          type = 2
        endelse
      endif $
      else if (ecc ge param.ecc_toll_max) then begin
        if (abs(la) ge param.loop_toll_min) then begin
          type = 2
        endif else begin
          type = 1
        endelse
      endif

      ; Assign the parameter to the reconstructed source
      smcsources[sorg].srctype=shape[type]
      smcsources[sorg].albedo_ratio = 0
      smcsources[sorg].srcheight = 0
      smcsources[sorg].srcx = mean_param[j,i,0]
      smcsources[sorg].srcy = mean_param[j,i,1]
      smcsources[sorg].srcpa = pa
      smcsources[sorg].eccen = ecc
      smcsources[sorg].srcfwhm = mean_param[j,i,4]
      smcsources[sorg].srcflux = mean_param[j,i,5]
      smcsources[sorg].loop_angle = la

      ; compute the visibilities
      srcparm   = hsi_vis_fwdfit_structure2array(smcsources[sorg], mapcenter)
      visxyrec  += hsi_vis_fwdfit_func(jdum, srcparm)

      ; Compute the number of ellipses, circles and loops (not required)
      prob_cel[sorg,*]=pc_pe_pl[j,i,*]
      prop_cel_tot=prob_cel[sorg,0]+prob_cel[sorg,1]+prob_cel[sorg,2] ;sum(prob_cel[sorg,*])
      prop_cel_100[sorg,*]=prob_cel[sorg,*]*100/ prop_cel_tot


    endfor

    ; Show the reconstructed image

    ;HSI_VIS_SOURCE2MAP, srcstrin, mapcenter, data,pixel=param.pixel_size, mapsize = param.fov/param.pixel_size ; _extra={pixel:pixel, mapsize : fov/pixel}
    ;map_asmc    = MAKE_MAP(data, xc=mapcenter[0], yc=mapcenter[1], dx=param.pixel_size, dy=param.pixel_size, time=anytim(time_title,/ecs))
     
    
     map_asmc = hsi_vis_bayes_showmap(param.fov,param.pixel_size, mapcenter, smcsources)

     window, 1
     plot_map, map_asmc, title= anytim(time_title,/ecs)


    ; Display the result
    print, ' '
    print, 'Estimated source(s) with the Bayesian approach:'
    HSI_VIS_FWDFIT_PRINT, smcsources

  endif

  return, smcsources
end
;----------------------------------------------------------------------------



;+
; NAME:
;   hsi_vis_bayes_sort_part
;
; PURPOSE:
;   Sort particles and their weights
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;
;-

function hsi_vis_bayes_sort_part, param, n_src, smcsources, sample, Nsources, weights=weights, wei_s=wei_s
  COMMON uvdata, u,v, pa, mapcenter

  j=0
  z_n=0
  c= make_array( double(n_src), param.N_param, param.N_particles, /double)
  wei=make_array(param.N_particles, /double)
  for i_p=0, param.N_particles-1 do begin
    if Nsources[i_p] eq n_src then begin
      c[*,*,j]= sample[0:n_src-1,*,i_p]
      wei[j]=weights[i_p]
      j=j+1
    endif
  endfor

  N_part_cond=j

  c_cond=c[*,*,0:N_part_cond-1]
  wei_s=wei[0:N_part_cond-1]


  par_s = make_array(n_src, param.N_param, N_part_cond, /double)



  if n_src eq 1 then begin

    par_s=c_cond


  endif else if n_src ge 2 then begin
    p_xyrf=[0,1]

    smp=make_array(param.N_param,double(n_src), /double)

    for i_sr=0, n_src-1 do begin

      smp[ 0,i_sr] = smcsources[i_sr].srcx
      smp[ 1,i_sr] = smcsources[i_sr].srcy
      smp[ 4,i_sr] = smcsources[i_sr].srcfwhm
      smp[ 5,i_sr] = smcsources[i_sr].srcflux
      smp[ 2,i_sr] = smcsources[i_sr].srcpa
      smp[ 3,i_sr] = smcsources[i_sr].eccen
      smp[ 6,i_sr] = smcsources[i_sr].loop_angle

    endfor

    dist_sm = make_array(n_src,n_src,   /double)

    for i_p= 0, N_part_cond-1 do begin
      non=0
      dist_sm = make_array(n_src,n_src,   /double)
      for j1=0, n_src-1 do begin
        for j2=0, n_src-1 do begin
          dist_sm[j1,j2] = norm(c_cond[j1,p_xyrf,i_p]-smp[p_xyrf,j2],lnorm=2)
        endfor
      endfor

      i_tmp=min(dist_sm, dimension=2,ii)
      iii= ii/n_src

      u_iii=iii(UNIQ(iii, SORT(iii))) ; uniq only with adjacent values

      if n_elements(iii) eq n_elements(u_iii) then begin
        for jj=0, n_src-1 do begin
          par_s[iii(jj),*,i_p] = c_cond[jj,*,i_p]
        endfor
      endif else begin
        z_n=z_n+1
        v_index=make_array(n_src, /double)
        d_index=make_array(n_src, /double)

        for jj=0, n_src-1 do begin
          q=where(iii[jj] eq iii, q_n)
          if q_n eq 1 then begin
            par_s[iii(jj),*,i_p] = c_cond[jj,*,i_p]
            v_index[iii(jj)]=1
            d_index[jj]=1
          endif
        endfor
        d_ind=where(d_index eq 0)
        dist_sm1=dist_sm[d_ind, *]
        v_ind=where(v_index eq 0)
        dist_sm2=dist_sm1[ *,v_ind]
        i_tmp2=min(dist_sm2,ww)
        ww2=ww mod (n_src - total(v_index))
        ww1=ww/(n_src-total(v_index))

        par_s[v_ind[ww1],*,i_p] = c_cond[d_ind(ww2),*,i_p]
        v_index[v_ind(ww1)]=0
        d_index[d_ind(ww2)]=0
        par_s[where(v_index eq 0),*,i_p] = c_cond[where(d_index eq 0),*,i_p]


      endelse



    endfor
  endif

  ;par_s=par_s[*,*,0:N_part_cond-z_n-1]


  return, par_s

end
;-------------------------------------------------------



;+
; NAME:
;   hsi_vis_bayes_func
;
; PURPOSE:
;   Compute posterior distribution of the number of sources and (if set) compute the
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;
;-

function hsi_vis_bayes_func, param, Nsources,  weights,fl,  sample, H, types, jdum,  shape,  visxyobs, visxyrec=visxyrec, N_src_expected=N_src_expected, smcsources, wei_s=wei_s, par_s=par_s, time_title

  ; Compute the posterior distribution of the number of sources:
  distr_num_src = make_array(param.N_max_sources+1)
  for num = 0,param.N_max_sources do begin
    for src = 0,param.N_particles-1 do begin
      if (Nsources[src] eq num) then begin
        distr_num_src[num] = distr_num_src[num] + weights[src]
      endif
    endfor
  endfor


  if fl eq 1 then begin ; when the algorithm converges, compute the recovered sources and the parameters of each particle.

    mean_param = hsi_vis_bayes_data_stack(param, Nsources,  sample,  weights, H, types, weight_type=weight_type, pc_pe_pl=pc_pe_pl)
    smcsources = hsi_vis_bayes_reconstr_vis(param, distr_num_src, H, jdum, mean_param,  shape,  visxyobs, visxyrec=visxyrec,  sample, N_src_expected=N_src_expected, weight_type, pc_pe_pl,  prop_cel_100=prop_cel_100, time_title)
    par_s = hsi_vis_bayes_sort_part(param, N_src_expected, smcsources, sample,  Nsources, weights=weights, wei_s=wei_s)

  endif

  return, distr_num_src
end
;----------------------------------------------------------------------------
