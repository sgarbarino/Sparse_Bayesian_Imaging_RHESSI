;+
; NAME:
;   stx_vis_bayes_prior 
;
; PURPOSE:
;   Draw the parameter of the sources from the prior distribution
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;-

function vis_bayes_prior, param, prior=prior, types=types, cerberus=cerberus 
  COMMON uvdata, u,v, pa, mapcenter

  ; Initialization 
  prior = make_array(param.N_max_sources, param.N_param, param.N_particles) 
  types = make_array(param.N_max_sources, param.N_particles,/integer)

  ; Generate the number of sources for each sample
  if cerberus eq 1 then Nsources = long(2 + make_array(param.N_particles)) 
  if cerberus eq 0 then Nsources = long(randomu(seed, param.N_particles, poisson=param.lam))
  
  ; If the number of sources is larger than param.N_max_sources then set it to the maximum
  index_tmp=where(Nsources gt param.N_max_sources, count_tmp)
  if count_tmp gt 0 then Nsources[index_tmp]= param.N_max_sources

  ; Draw from the prior the parameters for each particles
  for p = 0, param.N_particles-1 do begin

    if (Nsources[p] ge 1) then begin
      
      ;; Common parameters amongst all structures:
      ; position:
      xx = mapcenter[0]+(randomu(seed,Nsources[p])-0.5)*param.fov 
      yy = mapcenter[1]+(randomu(seed,Nsources[p])-0.5)*param.fov
      if cerberus eq 1 then begin
        xx[1]=-xx[0]
        yy[1]=-yy[0]
      endif
      ; fwhm and flux:
      fw = randomu(seed,Nsources[p])*param.maxfwhm
      fl = randomu(seed,Nsources [p])*param.maxflux
      
      ;; Initialise parameters specific to just some structures:
      ; source angle, loop angle, eccentricity and type
      angle = make_array(Nsources[p])
      loopang = make_array(Nsources[p])
      ecc = make_array(Nsources[p])
      type_tmp = make_array(Nsources[p],/integer)   
      
      ; Compute angle, loop, eccentricity.
      rand=randomu(seed, Nsources[p])

      ; Disk:
      index_tmp=where(rand le param.cdf_prior_types[0], count_tmp)
      if count_tmp gt 0 then type_tmp[index_tmp]=0 

      ; Elipse:
      index_tmp=where((rand gt param.cdf_prior_types[0]) and (rand le param.cdf_prior_types[1]), count_tmp)
      if count_tmp gt 0 then begin
        type_tmp[index_tmp]=1
        angle[index_tmp] = randomu(seed, count_tmp)*param.pa_max
        ecc[index_tmp] = randomu(seed, count_tmp)*0.9+0.1
      endif

      ; Loop:
      index_tmp=where((rand gt param.cdf_prior_types[1]) and (rand le param.cdf_prior_types[2]), count_tmp)
      if count_tmp gt 0 then begin
        type_tmp[index_tmp]=2
        angle[index_tmp] = randomu(seed, count_tmp)*param.pa_max
        ecc[index_tmp] = randomu(seed, count_tmp)*0.9+0.1
        for sorg = 0, count_tmp-1 do begin
          tmp_max_ecc = vis_bayes_compute_la(ecc[index_tmp[sorg]])
          loopang[index_tmp[sorg]] = (randomu(seed)*tmp_max_ecc*2-tmp_max_ecc) 
        endfor
      endif
 
      prior[0:Nsources[p]-1,*,p] = [[xx], [yy], [angle], [ecc], [fw], [fl], [loopang]]
      types[0:Nsources[p]-1,p] = type_tmp
      
    endif
  endfor

  return, Nsources

end