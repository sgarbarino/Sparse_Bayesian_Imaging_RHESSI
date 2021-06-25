;+
;
; NAME:
;   hsi_vis_bayes 
;
; PURPOSE:
;   Adaptive sequential Monte Carlo method to estimate the flare parameters from the visibilities.
;	The method approximates the posterior distribution, as given by Bayes theorem, for an a priori unknown
;		number of objects in the image; then individual source estimates are computed.
;
; CALLING SEQUENCE:
;   vis = ...
;   pixel_size = 1.
;   fov = 64.
;   map = hsi_vis_bayes(vis, fname, pixel_size, fov, lam, pC, pE, pL, N_particles, autoshape)
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


function hsi_vis_bayes, vis, fname, PIXEL_SIZE=pixel_size, FOV=fov, LAM=lam, PC=pC, PE=pE, PpL=pL, N_PARTICLES=N_particles, AUTOSHAPE=autoshape, PLOTHIST=plothist

  COMMON uvdata, u,v, pa, mapcenter,alb_apply_index
  COMMON srcshape, shape 


; If not set use the following values:
 default, pixel_size, 1.
 default, fov, 64.
 default, lam, 1.
 default, pC, 0.5
 default, pE, 0.25
 default, pL, 0.25
 default, N_particles, 5000.
 default, autoshape, 1
 default, final_plots, 1
 default, plothist, 1
 
   alb_apply_index = -1 ; added May 2018
   

 ; Visibilities:
  visin        = vis
  nvis         = N_ELEMENTS(visin)
  npt          = 2*nvis 
  jdum         = FINDGEN(npt)
  visx         = FLOAT(visin.obsvis)
  visy         = IMAGINARY(visin.obsvis)
  visxyobs     = [visx, visy]
  mapcenter    = visin[0].xyoffset 
  visxyobsorig = visxyobs
  visxyexp     = dblarr(N_elements(visxyobs))
  dummy        = hsi_vis_select (visin, PAOUT=paout)
  u            = visin.u
  v            = visin.v
  pa           = paout
  maxflux      = max(abs(visx))

  ; Adjust the input visibility errors to include a systematic term.
  syserr          = 0.05
  visamp          = ABS(visin.obsvis)                              ; input amplitudes
  visin.sigamp    = SQRT(visin.sigamp^2  + syserr^2 * visamp^2)
  sigma_noise     = [visin.sigamp, visin.sigamp] 

; If not set by the user, set the prior probability of having circles, ellipses, and loops (depending on the energy range)
 if keyword_set(autoshape) then begin
    if vis[0].erange[1] le 15. then begin
      pC=0.25
      pL=0.5
    endif else if vis[0].erange[1] gt 15. then begin
      pC=0.5
      pL=0.25
    endif
 endif

  Nsamp = 50.0
  ds = fov/Nsamp 
  xsamp = mapcenter[0]-fov/2. + dindgen(Nsamp)*ds ;x_0 = mapcenter[0]-fov/2
  ysamp = mapcenter[1]+fov/2. - dindgen(Nsamp)*ds ;  y_0 = mapcenter[1]+fov/2
  param = {N_particles:double(N_particles), fov:double(fov), lam:double(lam), N_max_sources:5., c_split:2./(60.*fov^3), Q_death:1.0/3., Q_birth:2.0/3.,  ratio_birth:'', $
           ratio_death:'', prior_types:[pC,pE,pL], maxflux:maxflux,  maxfwhm:fov/3.*2., cdf_prior_types:dblarr(3), pixel_size: pixel_size, $
          ecc_toll_min:0.1,  ecc_toll_max:0.4,  loop_toll_max : 45.,  loop_toll_min:10.,  N_param:7,  pa_max:360.,  la_max:179., $
          dmax:1e-1, dmin:1e-5, toll:1e-5, imin:0.9, imax:0.99, sigma_pos:'', sigma_fwhm : 5.0 ,sigma_flux:'', sigma_noise:sigma_noise, $
          sigma_pa:10.,  sigma_ecc:0.1,  sigma_loop:10., sigmas:dblarr(7),  sigmas_max:dblarr(7), Nsamp:Nsamp, ds:ds, xsamp:xsamp, ysamp:ysamp, nmaxiter:1500 }
  param.ratio_birth = param.Q_death/(1-param.Q_birth)
  param.ratio_death = 1.0/param.ratio_birth
  param.cdf_prior_types = cum_sum(param.prior_types)
  param.sigma_pos = param.fov/30.
  param.sigma_flux =  0.1 * param.maxflux
  shape = ['circle','ellipse','loop']
   param.sigmas = [param.sigma_pos, param.sigma_pos, param.sigma_pa, param.sigma_ecc, param.sigma_fwhm, param.sigma_flux,param.sigma_loop]
    param.sigmas_max = param.sigmas


  ; Prior and weights;
  Nsources=hsi_vis_bayes_prior(param,  prior=prior, types=types)
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
  ess = hsi_vis_bayes_ess(weights=weights, log_weights)

  progbar, progobj, /init

  while (expon le 1) do begin 
 
    j = j+1
    
    ; Compute the number of sources
    num_sources [*, j-1] = hsi_vis_bayes_func(param, Nsources, weights,0)
    
    ; Save values of the j-th iteration
    old_sample = sample  
    old_types = types

    ; Compute the ESS 
    ess=hsi_vis_bayes_resampling(param, sample=sample, types=types, weights=weights, ess, new_sample=new_sample, log_weights=log_weights, Nsources=Nsources)

    ; Perform the Monte Carlo moves
    new_log_weights=hsi_vis_bayes_moves(param, visxyobs, expon,  p, Nsources=Nsources, sample=sample, types=types, new_sample=new_sample, old_vis=old_vis, old_types, old_sample, incr, log_weights, like_unit=like_unit,  new_log_weights)

    ; Compute the ESS and the exponent 
    ess_new = hsi_vis_bayes_ess(weights=weights, new_log_weights)
    ess = hsi_vis_bayes_adapt_exponent(param, expon=expon, ess_new=ess_new, ess, like_unit, new_log_weights=new_log_weights,  log_weights, print_esp=print_esp , weights=weights,Nsources,  sample, H=H)
    
    progbar,  progobj, /update, percent=(expon*100.)
    progbar,  progobj, cancel=cancelled

  end 
   ; Compute the number of sources, find the optimal particle and show the map.
     distr_num_src = hsi_vis_bayes_func(param, Nsources,  weights,1,  sample, H, types, jdum,  shape,  visxyobs, visxyrec=visxyrec, N_src_expected=N_src_expected, smcsources, wei_s=wei_s, par_s=par_s, vis[0].trange[0])
  
;  b=asmc_hist(par_s, bins, wei_s, N_src_expected, distr_num_src)
;

  progbar,  progobj, /destroy
  ; Create the structure to be returned:
  bayes_sol={par_s:par_s,  wei_s:wei_s, smcsources:smcsources}
    
  save, par_s,N_src_expected, smcsources,  Nsources,  wei_s, distr_num_src, FILENAME = fname
  
  if plothist EQ 1 then begin
  res = hsi_vis_bayes_histograms(par_s, wei_s)
  endif
  
  return, bayes_sol
end

