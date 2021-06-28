;+
; NAME:
;    COMBIGEN
;
; PURPOSE:
;    Generates all possible combinations n-choose-k.
;
; CATEGORY:
;    Math
;
; CALLING SEQUENCE:
;    Result = COMBIGEN(N, K)
;
; INPUTS:
;    N:    Maximum number.
;
;    K:    Number of elements in each combination.
;
; OUTPUTS:
;    Returns a M x K array of all K-length combinations of numbers from 0 to N-1.
;
; EXAMPLE:
;    Generate all combinations 5-choose-3:
;
;    IDL> print, combigen(5,3)
;           0       0       0       0       0       0       1       1       1       2
;           1       1       1       2       2       3       2       2       3       3
;           2       3       4       3       4       4       3       4       4       4
;
; MODIFICATION HISTORY:
;    Written by Jeremy Bailin
;    1 April 2011   Initial writing.
;-
function combigen, n, k

  if n lt 2 then message, 'N must be greater than 1.'
  if k gt n then message, 'K must be less than or equal to N.'

  possible_prev = indgen(n)
  for vi=1l, k-1 do begin
    possible_next = indgen(n)
    ; only really possible when possible_next > possible_prev
    prevsize = size(possible_prev, /dimen)
    possible_prev = rebin(reform(possible_prev,1,prevsize), [n,prevsize], /sample)
    possible_next = rebin(possible_next, [n,prevsize], /sample)
    good = where(possible_next gt possible_prev, ngood)

    if vi eq 1 then begin
      ; set up array that answers get stored in
      result = [[possible_prev[good]],[possible_next[good]]]
    endif else begin
      ; add to result
      good_ai = array_indices(possible_prev, good)
      result = [[result[good_ai[1,*],*]], [possible_next[good]]]
    endelse

    possible_prev = possible_next[good]
  endfor

  return, result

end


;+
; NAME:
;   vis_bayes_accept_source
;
; PURPOSE:
;   Accept or reject the new particle
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;
;-

function vis_bayes_accept_source, likeratio, Nsources, new_types, proposta,p,  new_sample=new_sample, sample=sample, types=types, dato_new, old_vis=old_vis, phase = phase


  if (randomu(seed) le likeratio) then begin
    types[0:Nsources[p]-1,p] = new_types
    sample[*,*,p] = proposta
    new_sample[*,*,p] = proposta
    old_vis = dato_new
  endif

  return, Nsources
end
;----------------------------------------------------------------------------



;+
; NAME:
;   vis_bayes_birth_source
;
; PURPOSE:
;   Perform the birth move
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;
;-

function vis_bayes_birth_source, param, visxyobs, expon, p, Nsources=Nsources, sample,  n_circles, n_ellipses, n_loops, types=types, new_sample=new_sample, old_vis=old_vis, phase=phase

  COMMON uvdata, u,v, pa, mapcenter,alb_apply_index
  COMMON srcshape, shape
  ; Determine the shape and the parameters of the proposed new source.
  newx = mapcenter[0]+(randomu(seed)-0.5)*param.fov
  newy = mapcenter[1]+(randomu(seed)-0.5)*param.fov
  newfw = randomu(seed)*param.maxfwhm
  newfl = randomu(seed)*param.maxflux
  rand = randomu(seed)
  
  if (rand le param.cdf_prior_types[0]) then begin
    type = 0
    newsrc = transpose([newx,newy,0,0,newfw,newfl,0])
  endif $
  else if (rand gt param.cdf_prior_types[0]) and (rand le param.cdf_prior_types[1]) then begin
    type = 1
    newsrc = transpose([newx,newy,randomu(seed)*param.pa_max,randomu(seed)*0.9+0.1,newfw,newfl,0])
  endif $
  else begin
    type = 2
    newecc = randomu(seed)*0.9+0.1
    tmp_max_ecc = vis_bayes_compute_la(newecc)
    newloop = (randomu(seed)*tmp_max_ecc*2-tmp_max_ecc)
    newsrc = transpose([newx,newy,randomu(seed)*param.pa_max, newecc,newfw,newfl,newloop])
  endelse

  if (Nsources[p] ge 1) then begin
    proposed_sample = [sample[0:Nsources[p]-1,*,p], newsrc]
    types_proposed_sample = [types[0:Nsources[p]-1,p],type]
  endif $
  else if (Nsources[p] eq 0) then begin
    proposed_sample = newsrc
    types_proposed_sample = type
  endif

  ; compute the visibilities of the proposed source
  data2 = vis_bayes_compute_vis(types_proposed_sample, proposed_sample, Nsources[p]+1, phase)
  ; compute the loglikelihood ratio
  logratio = vis_bayes_compute_logratio(old_vis, data2, param.sigma_noise, visxyobs, expon)
  jumpratio = (param.lam/(Nsources[p]+1))*param.ratio_birth*exp(logratio)
  case type of
    0: jumpratio = jumpratio/(n_circles + 1)
    1: jumpratio = jumpratio/(n_ellipses + 1)
    2: jumpratio = jumpratio/(n_loops + 1)
  endcase

  ; Decide either to accept or not the proposed source.
  if (randomu(seed) le jumpratio) then begin
    sample[Nsources[p],*,p] = newsrc
    new_sample[Nsources[p],*,p] = newsrc
    types[Nsources[p],p] = type
    Nsources[p] = Nsources[p] + 1
    old_vis = data2
  endif

  return, sample
end
;----------------------------------------------------------------------------



;+
; NAME:
;   vis_bayes_death_source
;
; PURPOSE:
;   Perform the death move
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;
;-

function vis_bayes_death_source, param, visxyobs, expon,  p, Nsources=Nsources, sample,  n_circles, n_ellipses, n_loops, types=types, new_sample=new_sample, old_vis=old_vis, phase=phase

  COMMON uvdata, u,v, pa, mapcenter,alb_apply_index
  COMMON srcshape, shape

  if (Nsources[p] gt 1) then begin ; if there are more than one source
    ; Propose a source to be eliminated
    JJ = [(RANDOMU(SEED,/LONG) MOD NSOURCES[P])]
    ARRAY=SAMPLE[0:NSOURCES[P]-1,*,P]
    DIMS = SIZE(ARRAY, /DIMENSIONS)
    proposed_sample = array[Where(~Histogram(jj, MIN=0, MAX=dims[0]-1)),*] ;removecols(sample[0:Nsources[p]-1,*,p],jj)
    array = types[0:Nsources[p]-1,p]
    dims = Size(array, /DIMENSIONS)
    types_proposed_sample = array[Where(~Histogram(jj, MIN=0, MAX=dims[0]-1))] ; removeind(types[0:Nsources[p]-1,p],jj)

    ; Compute the corresponding visibilities and the loglikelihood ratio.
    data2 =  vis_bayes_compute_vis(types_proposed_sample,proposed_sample, Nsources[p]-1, phase)
    logratio = vis_bayes_compute_logratio(old_vis, data2, param.sigma_noise, visxyobs, expon)
    jumpratio = (Nsources[p]/param.lam)*param.ratio_death*exp(logratio)

    case types[jj,p] of
      0: jumpratio = jumpratio*n_circles
      1: jumpratio = jumpratio*n_ellipses
      2: jumpratio = jumpratio*n_loops
    endcase

    ; Accept or reject
    if (randomu(seed) le jumpratio) then begin
      sample[*,*,p] = [proposed_sample,make_array(param.N_max_sources-Nsources[p]+1,param.N_param)]
      new_sample[*,*,p] = sample[*,*,p]
      types[*,p] = [types_proposed_sample,make_array(param.N_max_sources-Nsources[p]+1,/integer)]
      Nsources[p] = Nsources[p] - 1
      old_vis = data2

    endif

  endif $
  else if (Nsources[p] eq 1) then begin ; if there is only one source:
    proposed_sample = make_array(param.N_param)
    types_proposed_sample = 0
    ; Compute the corresponding visibilities and the loglikelihood ratio.
    data2 = vis_bayes_compute_vis(types_proposed_sample, proposed_sample, Nsources[p]-1, phase)
    logratio = vis_bayes_compute_logratio(old_vis, data2, param.sigma_noise, visxyobs, expon)
    jumpratio = (Nsources[p]/param.lam)*param.ratio_death*exp(logratio)


    ; Accept or reject
    if (randomu(seed) le jumpratio) then begin
      sample[0,*,p] = proposed_sample
      new_sample[0,*,p] = proposed_sample
      types[0,p] = types_proposed_sample
      Nsources[p] = Nsources[p] - 1
      old_vis = data2

    endif

  endif

  return, sample

end
;----------------------------------------------------------------------------



;+
; NAME:
;   vis_bayes_change_sources
;
; PURPOSE:
;   Perform the source changing move: ellipse <-> cirlce and ellipse <-> loop
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;
;-

function vis_bayes_change_sources, param, Nsources, p,  sample=sample, types=types,  visxyobs, expon, new_sample=new_sample, old_vis=old_vis, phase=phase



  ; CIRCLE <---> ELLIPSE
  circles = where(types[0:Nsources[p]-1,p] eq 0, n_circles)
  ellipses = where(types[0:Nsources[p]-1,p] eq 1, n_ellipses)

  ; Circle --> Ellipse
  if (n_circles ge 1) then begin
    ; Set the new three parameters: eccentricity, pa, and loop angle.
    ind_src = (randomu(seed,/long) mod n_circles)
    src = circles[ind_src]
    new_types = types[0:Nsources[p]-1,p]
    new_types[src] = 1 ; circle is transformed into an ellipse
    proposed_sample = sample[*,*,p] ; use the same first four parameters of the circle
    proposed_sample[src,2] = randomu(seed)*param.pa_max ; ellipse angle computed with the uniform.
    rand = randomu(seed)
    proposed_sample[src,3] = 1-sqrt(1-rand)  ; ecc is computed by using the inverse transform of g(u) = 2(1-u)
    ; Compute the visibilities and the loglikelihood ratio
    new_vis = vis_bayes_compute_vis(new_types,proposed_sample,Nsources[p], phase)
    logratio = vis_bayes_compute_logratio(old_vis, new_vis, param.sigma_noise, visxyobs, expon)
    likeratio = exp(logratio)*(n_circles/((n_ellipses+1)*(2*(1-rand))))*(param.prior_types[1]/param.prior_types[0])
    ; Accept or reject the proposal
    Nsources= vis_bayes_accept_source(likeratio, Nsources, new_types, proposed_sample, p, new_sample=new_sample, sample=sample, types=types, new_vis, old_vis=old_vis, phase=phase)

  endif

  circles = where(types[0:Nsources[p]-1,p] eq 0, n_circles)
  ellipses = where(types[0:Nsources[p]-1,p] eq 1, n_ellipses)

  ; Ellipse --> Circle
  if (n_ellipses ge 1) then begin
    ; choose the ellipse to be converted into a circle
    ind_src = (randomu(seed,/long) mod n_ellipses)
    src = ellipses[ind_src]
    new_types = types[0:Nsources[p]-1,p]
    new_types[src] = 0 ; ellipse -> circle
    proposed_sample = sample[*,*,p]
    proposed_sample[src,2] = 0
    proposed_sample[src,3] = 0
    ; Compute the visibilities and the loglikelihood ratio
    new_vis = vis_bayes_compute_vis(new_types,proposed_sample,Nsources[p], phase)
    logratio = vis_bayes_compute_logratio(old_vis, new_vis, param.sigma_noise, visxyobs, expon)
    likeratio = exp(logratio)*(n_ellipses*(2*(1-sample[src,3,p]))/(n_circles+1))*(param.prior_types[0]/param.prior_types[1])
    ; Accept or reject the proposal
    Nsources= vis_bayes_accept_source(likeratio, Nsources, new_types, proposed_sample, p, new_sample=new_sample, sample=sample, types=types, new_vis, old_vis=old_vis , phase=phase)

  endif

  ; ELLIPSE <---> LOOP
  ; ------------------
  loop = where(types[0:Nsources[p]-1,p] eq 2, n_loops)
  ellipses = where(types[0:Nsources[p]-1,p] eq 1, n_ellipses)

  ; Ellipse -> loop
  if (n_ellipses ge 1) then begin
    ; Select the ellipse to be converted into a loop and compute the loop angle
    ind_src = (randomu(seed,/long) mod n_ellipses)
    src = ellipses[ind_src]
    new_types = types[0:Nsources[p]-1,p]
    new_types[src] = 2 ; ellipse -> loop
    proposed_sample = sample[*,*,p] ; use the first 6 parameters of the ellipse
    temp=randomu(seed)
    la_tmp=vis_bayes_compute_la(proposed_sample[src,3])
    proposed_sample[src,6] = (la_tmp-sqrt(la_tmp^2-(abs(temp*la_tmp*2-la_tmp))^2))*(abs(temp-0.5)/(temp-0.5))
    ; Compute the visibilities and the loglikelihood ratio
    new_vis = vis_bayes_compute_vis(new_types,proposed_sample,Nsources[p], phase)
    logratio=vis_bayes_compute_logratio(old_vis, new_vis, param.sigma_noise, visxyobs, expon)
    likeratio = exp(logratio)*(n_ellipses/(n_loops+1))*(param.prior_types[2]/param.prior_types[1])
    ; Accept or reject the proposal
    Nsources= vis_bayes_accept_source(likeratio, Nsources, new_types, proposed_sample, p, new_sample=new_sample, sample=sample, types=types, new_vis, old_vis=old_vis, phase=phase)

  endif


  loop = where(types[0:Nsources[p]-1,p] eq 2, n_loops)
  ellipses = where(types[0:Nsources[p]-1,p] eq 1, n_ellipses)

  if (n_loops ge 1) then begin ; LOOP --> ELLIPSE
    ; Select the loop to be converted into an ellipse
    ind_src = (randomu(seed,/long) mod n_loops)
    src = loop[ind_src]
    new_types = types[0:Nsources[p]-1,p]
    new_types[src] = 1 ; loop -> ellipse
    proposed_sample = sample[*,*,p] ; use the first six parameters of the loop
    proposed_sample[src,6] = 0
    ; Compute the visibilities and the loglikelihood ratio
    new_vis = vis_bayes_compute_vis(new_types,proposed_sample,Nsources[p], phase)
    logratio=vis_bayes_compute_logratio(old_vis, new_vis, param.sigma_noise, visxyobs, expon)
    likeratio = exp(logratio)*(n_loops/(n_ellipses+1))*(param.prior_types[1]/param.prior_types[2])
    ; Accept or reject the proposal
    Nsources= vis_bayes_accept_source(likeratio, Nsources, new_types, proposed_sample, p, new_sample=new_sample, sample=sample, types=types, new_vis, old_vis=old_vis, phase=phase)


  endif
  return, Nsources

end
;----------------------------------------------------------------------------



;+
; NAME:
;   vis_bayes_parameters_split_ellipse
;
; PURPOSE:
;   Compute the parameters of the two proposed ellipses
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;
;-

function vis_bayes_parameters_split_ellipse, n_ellipses, ellipses, sample, ind_src,p, Nsources, pr=pr

  ; Pick an ellipse (ell) to be split into two ellipses
  ind_src = (randomu(seed,/long) mod n_ellipses)
  ell = ellipses[ind_src]

  ; Determine the parameters of the two new proposed ellipses
  sinpa = sin(sample[ell,2,p]*!DTOR)
  cospa = cos(sample[ell,2,p]*!DTOR)
  semiaxis_max = sample[ell,4,p]/(3*(1-sample[ell,3,p]^2)^0.25)
  dx = semiaxis_max*sinpa
  dy = semiaxis_max*cospa

  z0 = randomu(seed)*dx/2 + (3./4)*dx
  z1 = randomu(seed)*dy/2 + (3./4)*dy
  z2 = randomu(seed)*sample[ell,2,p]/2 + (1./4)*sample[ell,2,p]
  z3 = randomu(seed)*sample[ell,3,p]/4 + (1./4)*sample[ell,3,p]/2
  z4 = randomu(seed)*sample[ell,4,p]/4 + (1./4)*sample[ell,4,p]/2
  z5 = randomu(seed)*sample[ell,5,p]/4 + (1./4)*sample[ell,5,p]/2
  pr = (abs(dx)/2)*(abs(dy)/2)*(sample[ell,2,p]/2)*(sample[ell,3,p]/4)*(sample[ell,4,p]/4)*(sample[ell,5,p]/4)

  x1 = sample[ell,0,p] + z0
  y1 = sample[ell,1,p] - z1
  pa1 = sample[ell,2,p] - z2
  e1 = sample[ell,3,p]/2 - z3
  fw1 = sample[ell,4,p]/2 - z4
  fl1 = sample[ell,5,p]/2 - z5
  s1 = transpose([x1,y1,pa1,e1,fw1,fl1,0])

  x2 = sample[ell,0,p] - z0
  y2 = sample[ell,1,p] + z1
  pa2 = sample[ell,2,p] + z2
  e2 = sample[ell,3,p]/2 + z3
  fw2 = sample[ell,4,p]/2 + z4
  fl2 = sample[ell,5,p]/2 + z5
  s2 = transpose([x2,y2,pa2,e2,fw2,fl2,0])

  proposed_sample = sample[*,*,p]
  proposed_sample[ell,*] = s1
  proposed_sample[Nsources[p],*] = s2

  return, proposed_sample
end
;----------------------------------------------------------------------------



;+
; NAME:
;   vis_bayes_split_ellipse
;
; PURPOSE:
;   Perform the split of an ellipse
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;
;-

function  vis_bayes_split_ellipse, param, Nsources,  p,  visxyobs, expon, new_sample=new_sample, sample=sample, types=types, old_vis=old_vis, phase=phase

  if (Nsources[p] lt param.N_max_sources) then begin ; if there are less than the num. max of sources

    ellipses = where(types[0:Nsources[p]-1,p] eq 1, n_ellipses)

    if (n_ellipses ge 1 ) then begin ; if there is at least an ellipse
      ; Compute the parameters of the two new proposed ellipses
      proposed_sample= vis_bayes_parameters_split_ellipse(n_ellipses, ellipses, sample,ind_src,p, Nsources, pr=pr)
      new_types = [types[0:Nsources[p]-1,p],1]

      ; Compute the visibilities and the loglikelihood ratio
      new_vis = vis_bayes_compute_vis(new_types,proposed_sample,Nsources[p]+1, phase)
      logratio=vis_bayes_compute_logratio(old_vis, new_vis, param.sigma_noise, visxyobs, expon)
      likeratio = (param.lam*n_ellipses/(Nsources[p]+1))*param.c_split*pr*exp(logratio)

      ; Accept or reject
      if (randomu(seed) le likeratio) then begin
        types[Nsources[p],p] = 1
        sample[*,*,p] = proposed_sample
        new_sample[*,*,p] = proposed_sample
        Nsources[p] = Nsources[p] + 1
        old_vis = new_vis
      endif

    endif
  endif

  return, Nsources
end
;----------------------------------------------------------------------------



;+
; NAME:
;   vis_bayes_merge_ellipse
;
; PURPOSE:
;   Perform the merge between two ellipses
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;
;-

function vis_bayes_merge_ellipse, param, p,  visxyobs, expon, sample=sample, Nsources, types=types, new_sample=new_sample, old_vis=old_vis, phase=phase

  ellipses = where(types[0:Nsources[p]-1,p] eq 1, n_ellipses)

  if (n_ellipses ge 2) then begin ; if there are at least two ellipses

    numcombi = long(factorial(n_ellipses)/(factorial(2)*factorial(n_ellipses-2)))
    combi = combigen(n_ellipses,2)
    merge_before = 0

    for couple = 0, numcombi-1 do begin
      if (keyword_set(merge_before) eq 0) then begin  ; if the merge move has not been done before

        ; Check if the two ellipses are close enough to be merged
        ind1 = combi[couple,0]
        ind2 = combi[couple,1]
        ell1 = ellipses[ind1]
        ell2 = ellipses[ind2]

        x1 = sample[ell1,0,p]
        y1 = sample[ell1,1,p]
        x2 = sample[ell2,0,p]
        y2 = sample[ell2,1,p]

        dist = distance_measure([[x1,y1],[x2,y2]])
        ecc_merge = sample[ell1,3,p]+sample[ell2,3,p]
        fw_merge = sample[ell1,4,p]+sample[ell2,4,p]
        merge_toll_max = 5*fw_merge/(6*(1-ecc_merge^2)^0.25)
        merge_toll_min = fw_merge/(2*(1-ecc_merge^2)^0.25)

        if ((dist ge merge_toll_min) and (dist le merge_toll_max)) then begin

          ; Compute the parameters of the proposed ellipse
          merge_before = 1
          x0 = (x1+x2)/2
          y0 = (y1+y2)/2
          pa0 = (sample[ell1,2,p]+sample[ell2,2,p])/2
          e0 = ecc_merge
          fw0 = fw_merge
          fl0 = sample[ell1,5,p]+sample[ell2,5,p]

          s0 = transpose([x0,y0,pa0,e0,fw0,fl0,0])

          sinpa = sin(pa0*!DTOR)
          cospa = cos(pa0*!DTOR)
          max_semiaxis = fw0/(3*(1-e0^2)^0.25)
          dx = max_semiaxis*sinpa
          dy = max_semiaxis*cospa

          pr = (abs(dx)/2)*(abs(dy)/2)*(pa0/2)*(e0/4)*(fw0/4)*(fl0/4)

          proposed_sample = sample[0:Nsources[p]-1,*,p]
          proposed_sample[ell1,*] = s0
          dims = Size(proposed_sample, /DIMENSIONS)
          proposed_sample = proposed_sample[Where(~Histogram([ell2], MIN=0, MAX=dims[0]-1)),*]; removecols(proposed_sample,ell2)
          new_types = types[0:Nsources[p]-1,p]
          dims = Size(new_types, /DIMENSIONS)
          new_types = new_types[Where(~Histogram([ell2], MIN=0, MAX=dims[0]-1))]  ;removeind(new_types,ell2)

          ; Compute the visibilities and the loglikelihood ratio
          new_vis = vis_bayes_compute_vis(new_types,proposed_sample,Nsources[p]-1, phase)
          logratio=vis_bayes_compute_logratio(old_vis, new_vis, param.sigma_noise, visxyobs, expon)
          likeratio = (Nsources[p]/(param.lam*(n_ellipses-1)))*(1./(param.c_split*pr))*exp(logratio)

          ; Accept or reject
          if (randomu(seed) le likeratio) then begin
            types[0:Nsources[p]-2,p] = new_types
            types[Nsources[p]-1,p] = 0
            sample[0:Nsources[p]-2,*,p] = proposed_sample
            sample[Nsources[p]-1,*,p] = make_array(param.N_param)
            new_sample[*,*,p] = sample[*,*,p]
            Nsources[p] = Nsources[p] - 1
            old_vis = new_vis

          endif

        endif
      endif
    endfor
  endif

  return, Nsources
end
;----------------------------------------------------------------------------



;+
; NAME:
;   vis_bayes_update_parameters
;
; PURPOSE:
;   Update (slightly perturb) the parameters of each source
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;
;-

function vis_bayes_update_parameters, param, Nsources, types, p,  sample, visxyobs, expon, new_sample=new_sample,  old_vis=old_vis, phase=phase, cerberus=cerberus
  COMMON uvdata, u,v, pa, mapcenter,alb_apply_index

  ; parameter update move of the Metropolis - Hastings algorithm

  sig=param.sigmas
  for k = 0, Nsources[p]-1 do begin

    ; Adaptively determine the sigmas to perturb the parameters
    if (types[k,p] ne 0) then begin
      sig[2] = param.sigmas_max[2]
      new_sigma_ecc = param.sigmas_max[3]*(1-sample[k,3,p])    ; linear
      ;new_sigma_ecc = sigma_ecc*exp(-sample[k,3,p])*(1-sample[k,3,p]) ; exponential
      ;new_sigma_ecc = sigma_ecc*(1-sample[k,3,p])/(1+sample[k,3,p])   ; hyperbolic
      sig[3] = min([new_sigma_ecc,sample[k,3,p]/4])
    endif
    sig[4] = min([param.sigmas_max[4],sample[k,4,p]/4])
    sig[5] = min([param.sigmas_max[5],(param.maxflux-sample[k,5,p])/4,sample[k,5,p]/4])
    if (types[k,p] eq 2) then begin
      sig[6] = min([param.sigmas_max[6],abs((180-sample[k,6,p])/4),abs(sample[k,6,p]/4)])
    endif


    for i = 0, param.N_param-1 do begin

      if  (types[k,p] eq 2) $
        or $
        ( (types[k,p] eq 1) and (i ne 6) )  $
        or $
        ( (types[k,p] eq 0) and (i ne 2) and (i ne 3) and (i ne 6) ) $
        then begin
        new_sample[k,i,p] = randomn(seed)*sig[i] + sample[k,i,p]  ; update the parameters with the sigmas computed before

        if cerberus eq 1 then begin
          if (i le 1) then begin    
            if (k eq 0) then begin 
              new_sample[1,i,p] = -new_sample[0,i,p]
            endif else begin
              new_sample[0,i,p] = -new_sample[1,i,p]        
            endelse
          endif
        endif
        
        if (i eq 6) and (new_sample[k,i,p] ge param.la_max) then new_sample[k,i,p]=param.la_max ; project the loop angle into the correct interval
        if (i eq 6) and (new_sample[k,i,p] le -param.la_max) then new_sample[k,i,p]=-param.la_max

        if (i eq 4)  and (new_sample[k,i,p] le 2) then new_sample[k,i,p]=2. ; lower bound for the fwhm

        ; Compute the visibilities and the loglikelihood ratio
        new_vis = vis_bayes_compute_vis(types[0:Nsources[p]-1,p],new_sample[0:Nsources[p]-1,*,p], Nsources[p], phase)
        logratio=vis_bayes_compute_logratio(old_vis, new_vis, param.sigma_noise, visxyobs, expon)
        likeratio = exp(logratio)

        ; Accept or reject
        if (randomu(seed) le likeratio) then begin
          sample[k,i,p] = new_sample[k,i,p]
          if cerberus eq 1 then begin
            if (i le 1) then begin
              if (k eq 0) then begin
                sample[1,i,p] = - sample[0,i,p]
              endif else begin
                  sample[0,i,p] = - sample[1,i,p]
              endelse
            endif
          endif
          old_vis = new_vis

        endif $
        else begin
          new_sample[k,i,p] = sample[k,i,p]
        endelse

      endif


    endfor

    ; Re-assign the shape flag (types)
    if (new_sample[k,3,p] lt param.ecc_toll_min) then begin
      types[k,p] = 0
    endif $
    else if ((new_sample[k,3,p] ge param.ecc_toll_min) and (new_sample[k,3,p] lt param.ecc_toll_max)) then begin
      if (abs(new_sample[k,6,p]) lt param.loop_toll_max) then begin
        types[k,p] = 1
      endif else begin
        types[k,p] = 2
      endelse
    endif $
    else if (new_sample[k,3,p] ge param.ecc_toll_max) then begin
      if (abs(new_sample[k,6,p]) ge param.loop_toll_min) then begin
        types[k,p] = 2
      endif else begin
        types[k,p] = 1
      endelse
    endif



  endfor


  return, sample
end
;----------------------------------------------------------------------------



;+
; NAME:
;   vis_bayes_compute_incremental_weights
;
; PURPOSE:
;   Compute incremental weights
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;
;-

function vis_bayes_compute_incremental_weights, param, visxyobs, old_nsource, old_types, old_sample, p, incr, log_weights, like_unit=like_unit,  new_log_weights, phase=phase

  d=param.dmax
  vis_tmp = dblarr(N_elements(visxyobs))
  if (old_nsource ge 1) then begin ; The (j+1)-th weight depends only on the old_sample (j-th sample).
    vis_tmp = vis_bayes_compute_vis(old_types,old_sample[0:old_nsource-1,*,p], old_nsource, phase)
  endif

  like_unit[p] = -(norm((visxyobs - vis_tmp)/param.sigma_noise))^2
  incr[p] = (d/2.)*like_unit[p]

  new_log_weights[p] = log_weights[p] + incr[p]

  return, new_log_weights
end
;----------------------------------------------------------------------------



;+
; NAME:
;   vis_bayes_birth_death_sources_post
;
; PURPOSE:
;   Call the birth and death functions
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;
;-

function vis_bayes_birth_death_sources, param, visxyobs, expon, p, Nsources, sample=sample,  types=types, new_sample=new_sample, old_vis=old_vis, phase=phase

  ; find the number of circles, ellipses, loops and cerberi
  if (Nsources[p] ge 1) then begin
    cerchi = where(types[0:Nsources[p]-1,p] eq 0, n_circles)
    ellissi = where(types[0:Nsources[p]-1,p] eq 1, n_ellipses)
    loop = where(types[0:Nsources[p]-1,p] eq 1, n_loops)
    cerberi = where(types[0:Nsources[p]-1,p] eq 1, n_cerberi)
  endif else begin ; Do we need it?
    n_circles = 0
    n_ellipses = 0
    n_loops = 0
    n_cerberi = 0
  endelse

  q = randomu(seed)

  ; birth move
  if (q ge param.Q_birth) and (Nsources[p] lt param.N_max_sources) then begin
    sample=vis_bayes_birth_source(param, visxyobs, expon, p, Nsources=Nsources, sample, n_circles, n_ellipses, n_loops, types=types, new_sample=new_sample, old_vis=old_vis, phase=phase)
  endif $
    ; death move
  else if (q le param.Q_death) and (Nsources[p] ge 1) then begin
    sample=vis_bayes_death_source(param, visxyobs, expon, p, Nsources=Nsources, sample,  n_circles, n_ellipses, n_loops, types=types, new_sample=new_sample, old_vis=old_vis, phase=phase)
  endif

  return, Nsources

end
;----------------------------------------------------------------------------



;+
; NAME:
;   vis_bayes_moves
;
; PURPOSE:
;   Perform the Monte Carlo moves (birth, death, merge, split and update) and compute the weights
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;
;-

function vis_bayes_moves, param, visxyobs, expon,  p, Nsources=Nsources, sample=sample, types=types, new_sample=new_sample, $
  old_vis=old_vis, phase=phase, cerberus=cerberus, old_types, old_sample, incr, log_weights, like_unit=like_unit,  new_log_weights
   
  COMMON uvdata, u,v, pa, mapcenter, alb_apply_index

  for p = 0, param.N_particles-1 do begin
    ; Save for each particle p, the number of sources and the visibilities at the j-th iteration.
    old_nsource = Nsources[p]
    if (Nsources[p] ge 1) then begin
      old_vis = vis_bayes_compute_vis(types[0:Nsources[p]-1,p],sample[0:Nsources[p]-1,*,p], Nsources[p], phase)
    endif $
    else if (Nsources[p] eq 0) then begin
      old_vis = make_array(2*size(pa,/dimension), /double)
      if (phase eq 0) then old_vis_complex = complex(old_vis[0:N_ELEMENTS(old_vis)/2-1],old_vis[N_ELEMENTS(old_vis)/2:*])
      if (phase eq 0) then old_vis = abs(old_vis_complex)
    endif

    if (cerberus eq 1) then begin
      
      ; Update move:
      sample = vis_bayes_update_parameters(param, Nsources, types, p, sample, visxyobs, expon, new_sample=new_sample, old_vis=old_vis, phase=phase, cerberus=cerberus)
      ; Compute incremental weights
      new_log_weights = vis_bayes_compute_incremental_weights(param,visxyobs, old_nsource, old_types, old_sample, p, incr, log_weights, like_unit=like_unit,  new_log_weights, phase=phase)
    
    endif else begin
      
      ; Perform the birth-death move:
      Nsources = vis_bayes_birth_death_sources(param, visxyobs, expon,  p, Nsources, sample=sample, types=types, new_sample=new_sample, old_vis=old_vis, phase=phase)
      if (Nsources[p] ge 1) then begin ; if the particle has at least one source, perform the other moves.
        ; Change move:
        Nsources = vis_bayes_change_sources(param, Nsources, p, sample=sample, types=types,  visxyobs, expon,  new_sample=new_sample, old_vis=old_vis, phase=phase)
        ; Split move:
        Nsources = vis_bayes_split_ellipse(param, Nsources, p, visxyobs, expon, new_sample=new_sample, sample=sample, types=types, old_vis=old_vis, phase=phase)
        ; Merge move:
        Nsources = vis_bayes_merge_ellipse(param, p, visxyobs, expon, sample=sample, Nsources, types=types, new_sample=new_sample, old_vis=old_vis, phase=phase)
        ; Update move:
        sample = vis_bayes_update_parameters(param, Nsources, types, p, sample, visxyobs, expon, new_sample=new_sample, old_vis=old_vis, phase=phase, cerberus=cerberus)
        ; Compute incremental weights
        new_log_weights = vis_bayes_compute_incremental_weights(param,visxyobs, old_nsource, old_types, old_sample, p, incr, log_weights, like_unit=like_unit,  new_log_weights, phase=phase)
      endif

    endelse
  endfor

  return, new_log_weights

end
