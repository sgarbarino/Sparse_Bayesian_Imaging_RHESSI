;+
;
; NAME:
;   stx_amp_bayes
;
; PURPOSE:
;   Adaptive sequential Monte Carlo method to estimate the flare parameters from visibility amplitudes.
;   The method approximates the posterior distribution, as given by Bayes theorem.
;
; CALLING SEQUENCE:
;   bayes_sol = stx_vis_bayes(type, ampobs, sigamp)
;
; INPUTS:
;   type: parametric shape to use for the forward fitting method
;         - 'circle' : Gaussian circular source
;         - 'ellipse': Gaussian elliptical source
;         - 'multi'  : double Gaussian circular source
;   
;   ampobs: array containing the values of the observed visibility amplitudes
;   sigamp: array containing the values of the errors on the observed visibility amplitudes
;
; KEYWORDS:
;   FOV: field of view (in arcsec) of the image (default set to 64)
;   N_PARTICLES: number of Monte Carlo samples (or "particles", default is 5000)
;   PLOTHIST: if set to 1, the histograms of the probability distributions are plotted
;
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
; - Massa, P. et al. (2021). Imaging from STIX visibility amplitudes
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;   March 2021 S.Garbarino adapted to visibility amplitudes 
;
; CONTACT:
;   garbarino [at] dima.unige.it
;   massa.p [at] dima.unige.it


function stx_amp_bayes, type, ampobs, sigamp, FOV=fov, N_PARTICLES=N_particles, PLOTHIST=plothist
  
  COMMON uvdata, u, v, pa, mapcenter, alb_apply_index
  COMMON srcshape, shape

  ; If not set use the following values:
  default, pixel_size, 1.
  default, fov, 64.
  default, lam, 1.
  default, N_particles, 5000.
  default, autoshape, 0
  default, final_plots, 1
  default, plothist, 1
  default, phase, 0
  
  alb_apply_index = -1 ; added May 2018
  
  case type of
    
    'circle': begin
              pE = 0
              pC = 1
              pL = 0
              cerberus=0
              end
    
    'ellipse': begin
               pE = 1
               pC = 0
               pL = 0
               cerberus=0
               end
      
     'multi': begin
              pE = 0
              pC = 0
              pL = 0
              cerberus=1
              end
    
  endcase
  

  ; Amplitudes:
  namp         = N_ELEMENTS(ampobs)
  npt          = 2*namp
  jdum         = FINDGEN(npt)
  
  visxyobs     = ampobs ; data
  maxflux      = 2 * max(visxyobs)

  mapcenter    = [0., 0.]
  visxyobsorig = visxyobs
  visxyexp     = dblarr(N_elements(visxyobs))

  ; Position angle
  twopi = 2.*!pi
  pa = ((atan(v, u) + twopi) mod twopi) * !radeg

  sigma_noise     = sigamp/2

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
    prior_types:[pC,pE,pL],$
    maxflux:maxflux, maxfwhm:fov/3.*2., $
    cdf_prior_types:dblarr(3), $ 
    ;pixel_size: pixel_size, $
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

  distr_num_src = vis_bayes_func(param, Nsources,  weights,1,  sample, H, types, jdum, shape,  visxyobs, visxyrec=visxyrec, N_src_expected=N_src_expected, smcsources, wei_s=wei_s, par_s=par_s, phase=phase, cerberus=cerberus)

  progbar,  progobj, /destroy
  
 ; 1) (for all) return position angle in interval [0,180] 
 dim = size(par_s)
 for k = 0, dim[1]-1 do begin
    par_s[k,2,*] = ((par_s[k,2,*] mod (180)) + 180) mod (180) 
 endfor

 ;2) if cerberus, remove the ambiguity related to the simmetry of the problem by 
 ; forcing assignment of one of the two possible (symmetric) configurations
  if cerberus eq 1 then begin

    ; Identify the two sources via clustering on the fluxes
    tmp_par = par_s[*,[5],*]  
    par = reform([(tmp_par[0,*,*]), (tmp_par[1,*,*])])   
    result = kmeans_cluster(par, k = 2, initial = 'random')  
    ;iplot, par,  LINESTYLE = 6, SYM_INDEX = 4
    
    par_complete = [reform(par_s[0,*,*]), reform(par_s[1,*,*])]
    config_0 = par_complete[*,where(result.clusters eq 0)]
    config_1 = par_complete[*,where(result.clusters eq 1)]
    
    ; Switch flux and fwhm in one of the two sources  
    tmp = config_1
    config_1[5,*] = tmp[12,*]
    config_1[12,*] = tmp[5,*]
    config_1[4,*] = tmp[11,*]
    config_1[11,*] = tmp[4,*]
    
    ; Merge the two configurations
    config = transpose([transpose(config_0), transpose(config_1)])
    
    ; Re-assign par_s
    par_s[0,*,*] = config[0:6,*]
    par_s[1,*,*] = config[7:13,*]
 
    ; Re-assign smcsources
    smcsources[0].srcx = mean(par_s[0,0,*])
    smcsources[0].srcy = mean(par_s[0,1,*])
    smcsources[0].srcflux = mean(par_s[0,5,*])
    smcsources[0].srcfwhm = mean(par_s[0,4,*])
    
    smcsources[1].srcx = mean(par_s[1,0,*])
    smcsources[1].srcy = mean(par_s[1,1,*])
    smcsources[1].srcflux = mean(par_s[1,5,*])
    smcsources[1].srcfwhm = mean(par_s[1,4,*])
 
    endif



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

;save, par_s, N_src_expected, smcsources,  Nsources,  wei_s, distr_num_src, FILENAME = fname

; Display the result
 print, ' '
 print, 'Estimated source(s) with the Bayesian approach:'
 
 case type of
  
   'circle': begin

             PRINT
             PRINT, 'COMPONENT    TYPE          FLUX         FWHM  '
             PRINT, '                         cts/s/keV     arcsec '
             PRINT
        
             temp        = [ smcsources[0].srcflux, smcsources[0].srcfwhm]
             PRINT, 1, smcsources[0].srctype, temp, FORMAT="(I5, A13, F13.2, 1F13.1)"
             temp        = [ stddev(par_s[0,5,*]),stddev(par_s[0,4,*])]
             PRINT, ' ', '(std)', temp, FORMAT="(A7, A11, F13.2, 1F13.1)"
             PRINT, ' '
        
             if plothist EQ 1 then begin
        
               res = vis_bayes_histograms_nopos(par_s, wei_s, type, dimensions = [600, 300])
        
             endif
        
             end
  
 
 'ellipse': begin
            
            PRINT
            PRINT, 'COMPONENT    TYPE          FLUX         FWHM      Eccentricity      Angle           [FWHM min]       [FWHM max]'
            PRINT, '                         cts/s/keV     arcsec                        deg             [arcsec]         [arcsec]'
            PRINT
        
            fwhm_min = smcsources[0].srcfwhm * (1. - smcsources[0].eccen^2)^0.25
            fwhm_max = smcsources[0].srcfwhm / (1. - smcsources[0].eccen^2)^0.25
            temp        = [ smcsources[0].srcflux, smcsources[0].srcfwhm,  smcsources[0].eccen, smcsources[0].srcpa + 90., fwhm_min, fwhm_max]
            PRINT, 1, smcsources[0].srctype, temp, FORMAT="(I5, A13, F13.2, 1F13.1, F12.1, 2F17.1, 2F19.1, 2F16.1)"
            temp        = [ stddev(par_s[0,5,*]),stddev(par_s[0,4,*]), stddev(par_s[0,3,*]), stddev(par_s[0,2,*])]
            PRINT, ' ', '(std)', temp, FORMAT="(A7, A11, F13.2, 1F13.1, F12.1, 2F17.1)"
            PRINT, ' '
            
            if plothist EQ 1 then begin
              
            par_s[*,2,*] += 90.
            par_s[*,2,*] = par_s[*,2,*] mod 180.
            res = vis_bayes_histograms_nopos(par_s, wei_s, type)
                

            endif 
              
            end 
  
   'multi': begin
    
            PRINT
            PRINT, 'COMPONENT  PROFILE      FLUX     X(+ve W)  Y(+ve N)    FWHM   '
            PRINT, '                      ph/cm2/s    arcsec    arcsec    arcsec  '
            PRINT
            
            nsrc = N_ELEMENTS(smcsources)
            FOR n = 0, nsrc-1 DO BEGIN
            temp        = [ smcsources[n].srcflux,  smcsources[n].srcx,  smcsources[n].srcy,  smcsources[n].srcfwhm]
            PRINT, n+1, smcsources[n].srctype, temp, FORMAT="(I5, A13, F12.2, 3F10.2, F10.3, 2F12.1)"
            
            temp        = [ stddev(par_s[n,5,*]),  stddev(par_s[n,0,*]),  stddev(par_s[n,1,*]),  stddev(par_s[n,4,*])]
            PRINT, n+1, '(std)', temp, FORMAT="(I5, A13, F12.2, 2F10.2, F10.3, 1F12.1)"
            
            ENDFOR
    
            if plothist EQ 1 then begin

            res = vis_bayes_histograms_cerberus(par_s, wei_s)

            endif
    
            end
  
 endcase

  return, bayes_sol
end

