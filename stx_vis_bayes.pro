;+
;
; NAME:
;   stx_vis_bayes
;
; PURPOSE:
;   Adaptive sequential Monte Carlo method to estimate the flare parameters from the visibilities.
; The method approximates the posterior distribution, as given by Bayes theorem, for an a priori unknown
;   number of objects in the image; then individual source estimates are computed.
;
; CALLING SEQUENCE:
;   vis = ...
;   pixel_size = 1.
;   fov = 64.
;   map = stx_vis_bayes(vis, fname, pixel_size, fov, lam, pC, pE, pL, N_particles, autoshape)
;
; INPUTS:
;   vis: input visibility structure in standard format
;   fname: name of the file where the output of the code will be saved. The file name should end with ".sav"
;
; KEYWORDS:
;   PIXEL_SIZE: size in arcsec (default set to 1)
;   FOV: field of view (in arcsec) of the image (default set to 64)
;   LAM: expected value of the number of objects in the image (default is 1)
;   PC: prior probability of having a circle (default set to 1/2)
;   PE: prior probability of having an ellipse (default set to 1/4)
;   PL: prior probability of having a loop(default set to 1/4)
;   N_PARTICLES: number of Monte Carlo samples (or "particles", default is 5000)
;   AUTOSHAPE: if set:
;                 (PC, PE, PL) = (1/4, 1/4, 1/2) if max(energy in keV)<=15,
;                 (PC, PE, PL) = (1/2, 1/4, 1/4) otherwise
;   FINAL_PLOTS: returns the map of the solution

;
; OUTPUTS:
;   map_asmc: image map in the structure format provided by the routine make_map.pro
;   fname.sav: the file includes the variables used in the post-processing step
;
; RETURNS:
;   bayes_sol: structure that includes the recovered particles (par_s), their weights (wei_s), and the reconstructed sources (smcsources)
;
;
; RESTRICTIONS:
;   -
;
; PAPER:
; - Sciacchitano, F., Sorrentino, A., Emslie, A. G., Massone, A. M., & Piana, M. (2018). Identification of Multiple Hard X-Ray Sources in Solar Flares: A Bayesian Analysis of the 2002 February 20 Event. The Astrophysical Journal, 862(1), 68.
; - Sciacchitano, F., Lugaro, S., & Sorrentino, A. (2019). Sparse Bayesian Imaging of Solar Flares. SIAM Journal on Imaging Sciences, 12(1), 319-343.
;
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;-


function stx_vis_bayes, vis, fname, PIXEL_SIZE=pixel_size, FOV=fov, LAM=lam, PC=pC, PE=pE, PLoop=pL, N_PARTICLES=N_particles, $
  AUTOSHAPE=autoshape, PLOTHIST=plothist, phase = phase, cerberus=cerberus

  COMMON uvdata, u,v, pa, mapcenter, alb_apply_index
  COMMON srcshape, shape

  ; If not set use the following values:
  default, pixel_size, 1.
  default, fov, 64.
  default, lam, 1.
  default, pC, 0.5
  default, pE, 0.5
  default, pLoop, 0.25
  default, N_particles, 5000.
  default, autoshape, 1
  default, final_plots, 1
  default, plothist, 1
  default, phase, 1
  default, cerberus, 0
  
  alb_apply_index = -1 ; added May 2018

  ; Visibilities:
  visin        = vis
  nvis         = N_ELEMENTS(visin)
  npt          = 2*nvis
  jdum         = FINDGEN(npt)
  ; If phase, I collect complex visibility as usual
  if (phase eq 1) then begin
    visx         = FLOAT(visin.obsvis) ; real part
    visy         = IMAGINARY(visin.obsvis) ; imaginary part
    visxyobs     = [visx, visy] ; data
    ampobs       = ABS(visin.obsvis) ; amplitude
    maxflux      = 2 * MAX(ABS(visx))
  endif
  ; if no phase, I only collect amplitude data
  if (phase eq 0) then begin
    ampobs       = FLOAT(visin.obsvis) ; amplitude
    visxyobs     = ampobs ; data
    maxflux      = 2 * max(visxyobs)
  endif
  mapcenter    = visin[0].xyoffset
  visxyobsorig = visxyobs
  visxyexp     = dblarr(N_elements(visxyobs))
  u            = visin.u
  v            = visin.v

  ; Position angle
  twopi = 2.*!pi
  pa = ((atan(v, u) + twopi) mod twopi) * !radeg

  ; Adjust the input visibility errors to include a systematic term.
  syserr          = 0.05
  sigamp    = SQRT(visin.sigamp^2  + syserr^2 * ampobs^2)
  if (phase eq 1) then sigma_noise     = [sigamp, sigamp]
  if (phase eq 0) then sigma_noise     = sigamp/2

  ; If not set by the user, set the prior probability of having circles, ellipses, loops and cerberus
  if keyword_set(autoshape) then begin
    pC=.5
    pLoop=.25
    pE=.5
  endif

  Nsamp = 50.0
  ds = fov/Nsamp
  xsamp = mapcenter[0]-fov/2. + dindgen(Nsamp)*ds ;  x_0 = mapcenter[0]-fov/2
  ysamp = mapcenter[1]+fov/2. - dindgen(Nsamp)*ds ;  y_0 = mapcenter[1]+fov/2

  param = {N_particles:double(N_particles), fov:double(fov), lam:double(lam), $
    N_max_sources:'', c_split:2./(60.*fov^3), Q_death:1.0/3., Q_birth:2.0/3., $
    ratio_birth:'', ratio_death:'', $
    prior_types:[pC,pE,pLoop],$
    maxflux:maxflux, maxfwhm:fov/3.*2., $
    cdf_prior_types:dblarr(3), $ 
    pixel_size: pixel_size, $
    ecc_toll_min:0.1,  ecc_toll_max:0.4,  loop_toll_max : 45.,  loop_toll_min:10., flux_toll_min:maxflux/10,$
    N_param:7, $ ;;7
    pa_max:180., la_max:179., dmax:1e-1, dmin:1e-5, $
    toll:1e-5, $
    imin:0.9, imax:0.99, $
    sigma_pos:'', sigma_fwhm : 5.0 ,sigma_flux:'', sigma_noise:sigma_noise, sigma_pa:10., sigma_ecc:0.1, sigma_loop:10., $
    sigmas:dblarr(7), sigmas_max:dblarr(7), $
    Nsamp:Nsamp, ds:ds, xsamp:xsamp, ysamp:ysamp, $
    nmaxiter:1500 }

  if cerberus eq 1 then begin
    param.N_max_sources = 2
  endif else begin
    if phase eq 0 then begin
      param.N_max_sources = 1 
    endif else begin
      param.N_max_sources = 5
    endelse
  endelse
      
  
  param.ratio_birth = param.Q_death/(1-param.Q_birth)
  param.ratio_death = 1.0/param.ratio_birth
  param.cdf_prior_types = cum_sum(param.prior_types)
  param.sigma_pos = param.fov/30.
  param.sigma_flux =  0.1 * param.maxflux

  shape = ['circle','ellipse','loop']
  param.sigmas = [param.sigma_pos, param.sigma_pos, param.sigma_pa, param.sigma_ecc, param.sigma_fwhm, param.sigma_flux, param.sigma_loop]
  param.sigmas_max = param.sigmas

  ; Prior and weights;
  Nsources = vis_bayes_prior(param,  prior=prior, types=types, cerberus=cerberus)
  sample = prior
  new_sample = sample

  ; Initialization
  new_log_weights = make_array(param.N_particles,/double)
  like_unit = make_array(param.N_particles,/double)
  incr = make_array(param.N_particles,/double)
  log_weights = make_array(param.N_particles,/double,value=alog(double(1)/param.N_particles))
  num_sources=make_array(param.N_max_sources+1, param.nmaxiter, /double)
  j = 0
  expon= 0
  print_esp = 0

  ; Compute the initial ESS
  ess = vis_bayes_ess(weights=weights, log_weights)

  progbar, progobj, /init

  while (expon le 1) do begin

    j = j+1

    ; Compute the number of sources
    num_sources [*, j-1] = vis_bayes_func(param, Nsources, weights, 0, phase=phase, cerberus=cerberus)

    ; Save values of the j-th iteration
    old_sample = sample
    old_types = types

    ; Compute the ESS
    ess = vis_bayes_resampling(param, sample=sample, types=types, weights=weights, ess, new_sample=new_sample, log_weights=log_weights, Nsources=Nsources)

    ; Perform the Monte Carlo moves
    new_log_weights = vis_bayes_moves(param, visxyobs, expon,  p, Nsources=Nsources, sample=sample, types=types, new_sample=new_sample, old_vis=old_vis, phase=phase, cerberus=cerberus, old_types, old_sample, incr, log_weights, like_unit=like_unit, new_log_weights)

    ; Compute the ESS and the exponent
    ess_new = vis_bayes_ess(weights=weights, new_log_weights)
    ess = vis_bayes_adapt_exponent(param, expon=expon, ess_new=ess_new, ess, like_unit, new_log_weights=new_log_weights,  log_weights, print_esp=print_esp , weights=weights,Nsources,  sample, H=H)

    progbar,  progobj, /update, percent=(expon*100.)
    progbar,  progobj, cancel=cancelled

  end
  
  ; Compute the number of sources and find the optimal particle (should be modified for phase data as mean/max is no optimal anymore)

  distr_num_src = vis_bayes_func(param, Nsources,  weights,1,  sample, H, types, jdum, shape,  visxyobs, visxyrec=visxyrec, N_src_expected=N_src_expected, smcsources, wei_s=wei_s, par_s=par_s, vis[0].time_range[0], phase=phase, cerberus=cerberus)

  progbar,  progobj, /destroy
  
 ; 1) (for all) return position angle in interval [0,180] 
 dim = size(par_s)
 for k = 0, dim[1]-1 do begin
    par_s[k,2,*] = ((par_s[k,2,*] mod (180)) + 180) mod (180) 
 endfor


  ;2) if cerberus, separate peaks for fwhm and flux and return mean as estimate. 
; if cerberus eq 1 then begin
;   tmp_fw = [reform(par_s[0,4,*]),reform(par_s[1,4,*])] ; for fwhm
;   tmp_fl = [reform(par_s[0,5,*]),reform(par_s[1,5,*])] ; for flux
;
;   tmp_comb = transpose([[tmp_fw], [tmp_fl]])
;   tmp_comb = transpose([tmp_fl])
;
;   ;weights = CLUST_WTS(tmp_comb, N_CLUSTERS = 2)  
;   ;result = reform(CLUSTER(tmp_comb, weights, N_CLUSTERS = 2))
;   result = kmeans_cluster(tmp_comb, k = 2, initial = 'random') 
;  
;   fw_0 = tmp_fw[where(result.clusters eq 0)]
;   fw_1 = tmp_fw[where(result.clusters eq 1)]
;   fl_0 = tmp_fl[where(result.clusters eq 0)]
;   fl_1 = tmp_fl[where(result.clusters eq 1)]
;
;    ;iplot, tmp_comb(where(result.clusters)), symbol = '*'
;    
;   ; re-assign par_s
;     par_s[0,4, *] = fw_0
;     par_s[1,4, *] = fw_1
;     par_s[0,5, *] = fl_0
;     par_s[1,5, *] = fl_1
;     
;   ; re-assign the flux/fwhm at the smcsources    
;   smcsources[0].srcfwhm = mean(fw_0)
;   smcsources[1].srcfwhm = mean(fw_1)
;   smcsources[0].srcflux = mean(fl_0)
;   smcsources[1].srcflux = mean(fl_1)
;    
;   endif
   
 ; 3) if no phase and no cerberus, return MAP estimates for pa and ecc
   if phase eq 0 and cerberus eq 0 then begin
     dimensions=size(par_s)
     tmp_pa = par_s(0,2,0:dimensions[3]-1)
     maximum = max(tmp_pa)
     minimum = min(tmp_pa)
     if maximum eq minimum then begin
       maximum = minimum +1
     endif
     pt_pa = min(tmp_pa) + (maximum-minimum) * findgen(11)/11
     pdf_pa = fltarr(11)
     FOR i = 0, dimensions[3]-1 DO BEGIN
       idx_pa = floor((par_s(0,2,i)-minimum)/(maximum-minimum)*10.9998)
       pdf_pa[idx_pa] = pdf_pa[idx_pa] + wei_s(i)
     endfor
     smcsources[0].srcpa = pt_pa(where(pdf_pa eq max(pdf_pa)))

     tmp_ecc = par_s(0,3,0:dimensions[3]-1)
     maximum = max(tmp_ecc)
     minimum = min(tmp_ecc)
     if maximum eq minimum then begin
       maximum = minimum +1
     endif
     pt_ecc = minimum + (maximum-minimum) * findgen(11)/11
     pdf_ecc = fltarr(11)
     FOR i = 0, dimensions[3]-1 DO BEGIN
       idx_ecc = floor((par_s(0,3,i)-minimum)/(maximum-minimum)*10.9998)
       pdf_ecc[idx_ecc] = pdf_ecc[idx_ecc] + wei_s(i)
     endfor
     smcsources[0].eccen = pt_ecc(where(pdf_ecc eq max(pdf_ecc)))
     
   endif

; Create the structure to be returned:

bayes_sol={par_s:par_s,  wei_s:wei_s, smcsources:smcsources}

save, par_s, N_src_expected, smcsources,  Nsources,  wei_s, distr_num_src, FILENAME = fname


; Display the result
 print, ' '
 print, 'Estimated source(s) with the Bayesian approach:'
 
 if phase eq 0 and cerberus eq 0 then begin  ; show mean and std without positions
  hsi_VIS_FWDFIT_PRINT_nopos, smcsources
  dim = size(par_s)
  FOR n = 0, dim[1]-1 DO BEGIN
    ; add std
    temp        = [ stddev(par_s[n,5,*]),  $
      stddev(par_s[n,4,*]),  stddev(par_s[n,3,*]), stddev(par_s[n,2,*]), stddev(par_s[n,6,*]) ]
    PRINT,n+1 , '(std)', temp, '--', '--', FORMAT="(I5, A13, F12.2, 1F10.2, F10.3, 2F12.1, A13, A13)"
  ENDFOR
 
 endif else if cerberus eq 1 then begin  ; re-compute standard deviation of fwhm and flux according to new histograms
   hsi_VIS_FWDFIT_PRINT, smcsources
   dim = size(par_s)
   ; add std
    

   temp        = [ stddev(par_s[0,5,*]), $ ;fl_0
     stddev(par_s[0,0,*]),  stddev(par_s[0,1,*]),  $
      stddev(par_s[0,4,*]),  $   ;fw_0
      stddev(par_s[0,3,*]), stddev(par_s[0,2,*]), stddev(par_s[0,6,*])]
   PRINT, 1, '(std)', temp, '--', '--', FORMAT="(I5, A13, F12.2, 2F10.2, F10.3, 1F12.1, F12.3, F12.1, A13, A13)"
    
   temp        = [ stddev(par_s[1,5,*]), $;fl_1
      stddev(par_s[1,0,*]),  stddev(par_s[1,1,*]),  $
      stddev(par_s[1,4,*]), $ ;fw_1
       stddev(par_s[1,3,*]), stddev(par_s[1,2,*]), stddev(par_s[1,6,*])]
   PRINT, 2, '(std)', temp, '--', '--', FORMAT="(I5, A13, F12.2, 2F10.2, F10.3, 1F12.1, F12.3, F12.1, A13, A13)"
    
 endif else begin ; standard mean (std)
   hsi_VIS_FWDFIT_PRINT, smcsources
   dim = size(par_s)
   ; add std
   FOR n = 0, dim[1]-1 DO BEGIN
     temp        = [ stddev(par_s[n,5,*]),  stddev(par_s[n,0,*]),  stddev(par_s[n,1,*]),  $
       stddev(par_s[n,4,*]),  stddev(par_s[n,3,*]), stddev(par_s[n,2,*]), stddev(par_s[n,6,*])]
     PRINT, n+1, '(std)', temp, '--', '--', FORMAT="(I5, A13, F12.2, 2F10.2, F10.3, 1F12.1, F12.3, F12.1, A13, A13)"
   ENDFOR
 endelse

  if plothist EQ 1 then begin
    if phase eq 1 then begin
      res = vis_bayes_histograms(par_s, wei_s)
    endif else begin
      if cerberus eq 1 then begin
        res = vis_bayes_histograms_cerberus(par_s, wei_s)    
      endif else begin
      res = vis_bayes_histograms_nopos(par_s, wei_s)      
      endelse
    endelse
    
  endif
    

  return, bayes_sol
end

