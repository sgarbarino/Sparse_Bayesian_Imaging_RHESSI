;+
;
; NAME:
;   stx_amp_pso
;
; PURPOSE:
;   Forward fitting method from visibility amplitudes based on Particle Swarm Optimization
;
; CALLING SEQUENCE:
;   stx_amp_pso, type, ampobs, sigamp, u, v, n_free
;
; INPUTS:
;   type: parametric shape to use for the forward fitting method
;         - 'circle' : Gaussian circular source
;         - 'ellipse': Gaussian elliptical source
;         - 'multi'  : double Gaussian circular source
;   
;   ampobs: array containing the values of the observed visibility amplitudes
;   sigamp: array containing the values of the errors on the observed visibility amplitudes
;   u: u coordinates of the sampling frequencies
;   v: v coordinates of the sampling frequencies
;   n_free: degrees of freedom (difference between the number of visibility amplitudes and the number of
;           parameters of the source shape)
;
; KEYWORDS:
;   SwarmSize: number of particles used in PSO (default is 100)
;   TolFun: tolerance for the stopping criterion (default is 1e-6)
;   maxiter: maximum number of iterations of PSO (defult is the product between of the numbers of parameters and the number of particles)
;   uncertainty: set to 1 for the computation of the parameters uncertainty (confidence strip approach)
;   silent: set to 1 for avoiding the print of the retrieved parameters
;
; HISTORY: January 2021, Massa P., Perracchione E. created
;
; CONTACT:
;   massa.p [at] dima.unige.it
;   perracchione [at] dima.unige.it

function stx_amp_pso, type, ampobs, sigamp, u, v, n_free, $
                     SwarmSize = SwarmSize, TolFun = TolFun, maxiter = maxiter, uncertainty = uncertainty, $
                     silent = silent

default, SwarmSize, 100.
default, TolFun, 1e-06
default, silent, 0

fun_name = 'amp_fwdfit_func'

aampobs = transpose(cmreplicate(ampobs, SwarmSize))
ssigamp = transpose(cmreplicate(sigamp, SwarmSize))
extra = {type: type, $
         ampobs: aampobs, $
         sigamp: ssigamp, $
         u: u, $
         v: v, $
         n_free: n_free}

estimate_flux = max(ampobs)

if type eq 'circle' then begin

  lb = [0.1*estimate_flux, 1.]
  ub = [1.5*estimate_flux, 100.]
  Nvars = n_elements(lb)

  if ~keyword_set(maxiter) then maxiter = Nvars*SwarmSize

  optim_f = swarmintelligence(fun_name, Nvars, lb, ub, SwarmSize, TolFun, maxiter, extra = extra)
  xopt = optim_f.xopt
  
  srcstr = {amp_src_structure}
  srcstr.srctype ='circle'

  fitsigmas = {amp_src_structure}
  fitsigmas.srctype ='std.dev'

  srcstr.srcflux         = xopt[0]
  srcstr.srcfwhm_max     = xopt[1]
  srcstr.srcfwhm_min     = xopt[1]

  if keyword_set(uncertainty) then begin
    
    print, ' '
    print, 'Uncertainty: '
    print, ' 
    
    ntry = 20
    namp = N_ELEMENTS(ampobs)

    trial_results = fltarr(2, ntry)
    for n=0,ntry-1 do begin
      testerror = RANDOMN(iseed, namp)          ; nvis element vector normally distributed with sigma = 1
      amptest   = ampobs + testerror * sigamp
      amptest = transpose(cmreplicate(amptest, SwarmSize))
      
      extra = {type: type, $
        ampobs: amptest, $
        sigamp: ssigamp, $
        u: u, $
        v: v, $
        n_free: n_free}
      
      optim_f = swarmintelligence(fun_name, Nvars, lb, ub, SwarmSize, TolFun, maxiter, extra = extra)
      xopt = optim_f.xopt

      trial_results[*,n]  = xopt

    endfor

    std_dev_par = stddev(trial_results, dimension=2)

    fitsigmas.srcflux         = std_dev_par[0]
    fitsigmas.srcfwhm_max     = std_dev_par[1]
    fitsigmas.srcfwhm_min     = std_dev_par[1]

  endif
  
endif


if type eq 'ellipse' then begin
  
  lb = [0.1*estimate_flux, 1., 0., -5.]
  ub = [1.5*estimate_flux, 100., 1., 5.]
  Nvars = n_elements(lb)

  if ~keyword_set(maxiter) then maxiter = Nvars*SwarmSize

  optim_f = swarmintelligence(fun_name, Nvars, lb, ub, SwarmSize, TolFun, maxiter, extra = extra)
  xopt = optim_f.xopt

  srcstr = {amp_src_structure}
  srcstr.srctype ='ellipse'

  fitsigmas = {amp_src_structure}
  fitsigmas.srctype ='std.dev'

  srcstr.srcflux = xopt[0]

  ecmsr = REFORM(SQRT(xopt[2]^2 + xopt[3]^2))
  eccen = SQRT(1 - EXP(-2*ecmsr))

  IF ecmsr GT 0 THEN srcstr.srcpa = reform(ATAN(xopt[3], xopt[2]) * !RADEG) + 90.

  srcstr.srcfwhm_min = xopt[1] * (1-eccen^2)^0.25
  srcstr.srcfwhm_max = xopt[1] / (1-eccen^2)^0.25
  

  if keyword_set(uncertainty) then begin

    print, ' '
    print, 'Uncertainty: '
    print, '

    ntry = 20
    namp = N_ELEMENTS(ampobs)
    
    trial_results = fltarr(Nvars, ntry)
    for n=0,ntry-1 do begin
      testerror = RANDOMN(iseed, namp)          ; nvis element vector normally distributed with sigma = 1
      amptest   = ampobs + testerror * sigamp
      amptest   = transpose(cmreplicate(amptest, SwarmSize))
      
      extra = {type: type, $
        ampobs: amptest, $
        sigamp: ssigamp, $
        u: u, $
        v: v, $
        n_free: n_free}

      optim_f = swarmintelligence(fun_name, Nvars, lb, ub, SwarmSize, TolFun, maxiter, extra = extra)
      xopt = optim_f.xopt
      

      ecmsr = REFORM(SQRT(xopt[2]^2 + xopt[3]^2))
      eccen = SQRT(1 - EXP(-2*ecmsr))
        
        IF ecmsr GT 0 THEN trial_results[3,n] = reform(ATAN(xopt[3], xopt[2]) * !RADEG) + 90.
        
        trial_results[0,n]  = xopt[0]
        trial_results[1,n]  = xopt[1] / (1-eccen^2)^0.25
        trial_results[2,n]  = xopt[1] * (1-eccen^2)^0.25
        
      endfor

      fitsigmas.srcflux         = stddev(trial_results[0, *])
      fitsigmas.srcfwhm_max     = stddev(trial_results[1, *])
      fitsigmas.srcfwhm_min     = stddev(trial_results[2, *])
      avsrcpa                   = ATAN(TOTAL(SIN(trial_results[3, *] * !DTOR)), $
                                  TOTAL(COS(trial_results[3, *] * !DTOR))) * !RADEG
      groupedpa                   = (810 + avsrcpa - trial_results[3, *]) MOD 180. 
      fitsigmas.srcpa          = STDDEV(groupedpa)

  endif

endif


if type EQ 'multi' then begin
  
  lb = [0.,  0.1*estimate_flux, 1e-3,  0.1*estimate_flux, 0., 0.]
  ub = [100., 1.5*estimate_flux, 1., 1.5*estimate_flux, 30., 180.]
  Nvars = n_elements(lb)
  
  if ~keyword_set(maxiter) then maxiter = Nvars*SwarmSize
  
  Nruns = 20
  xx_opt = []
  f = fltarr(Nruns)
  
  for i = 0,Nruns-1 do begin
    optim_f = swarmintelligence(fun_name, Nvars, lb, ub, SwarmSize, TolFun, maxiter, extra = extra)
    f[i] = optim_f.fopt
    xx_opt = [[xx_opt],optim_f.xopt]
  endfor
  
    dummy = min(f,location)
    xopt = xx_opt(location,*)
    
    srcstr = {amp_src_structure}
    srcstr.srctype ='ellipse'
    srcstr = amp_fwdfit_bifurcate(srcstr)
    
    fitsigmas = {amp_src_structure}
    fitsigmas.srctype ='std.dev'
    fitsigmas = amp_fwdfit_bifurcate(fitsigmas)
    
    srcstr[0].srcflux     = xopt[1]
    srcstr[0].srcfwhm_max     = xopt[0]
    srcstr[0].srcfwhm_min     = xopt[0]
    srcstr[0].srcx     = xopt[4] * cos(xopt[5] * !dtor)
    srcstr[0].srcy     = xopt[4] * sin(xopt[5] * !dtor)

    srcstr[1].srcflux     = xopt[3]
    srcstr[1].srcfwhm_max     = xopt[0]*xopt[2]
    srcstr[1].srcfwhm_min     = xopt[0]*xopt[2]
    srcstr[1].srcx     = -xopt[4] * cos(xopt[5] * !dtor)
    srcstr[1].srcy     = -xopt[4] * sin(xopt[5] * !dtor)
    
    
    if keyword_set(uncertainty) then begin

    print, ' '
    print, 'Uncertainty: '
    print, '

    ntry = 20
    namp = N_ELEMENTS(ampobs)

    trial_results = fltarr(Nvars, ntry)
    for n=0,ntry-1 do begin
      nn = n
      testerror  = RANDOMN(nn, namp)
      amptest    = ampobs + testerror * sigamp
      amptest   = transpose(cmreplicate(amptest, SwarmSize))
      
      extra = {type: type, $
        ampobs: amptest, $
        sigamp: ssigamp, $
        u: u, $
        v: v, $
        n_free: n_free}
        
      xx_opt = []
      f = fltarr(Nruns)

      for i = 0,Nruns-1 do begin
        optim_f = swarmintelligence(fun_name, Nvars, lb, ub, SwarmSize, TolFun, maxiter, extra = extra)
        f[i] = optim_f.fopt
        xx_opt = [[xx_opt],optim_f.xopt]
      endfor

      dummy = min(f,location)
      xopt = xx_opt(location,*)
        
      trial_results[0, n] = xopt[1]
      trial_results[1, n] = xopt[0]
      trial_results[2, n] = xopt[4] * cos(xopt[5] * !dtor)
      trial_results[3, n] = xopt[4] * sin(xopt[5] * !dtor)

      trial_results[4, n] = xopt[3]
      trial_results[5, n] = xopt[0]*xopt[2]
  
     endfor
     
     fitsigmas[0].srcflux     = stddev(trial_results[0, *])
     fitsigmas[0].srcfwhm_max     = stddev(trial_results[1,*])
     fitsigmas[0].srcfwhm_min     = stddev(trial_results[1,*])
     fitsigmas[0].srcx     = stddev(trial_results[2,*])
     fitsigmas[0].srcy     = stddev(trial_results[3,*])

     fitsigmas[1].srcflux     = stddev(trial_results[4,*])
     fitsigmas[1].srcfwhm_max     = stddev(trial_results[5,*])
     fitsigmas[1].srcfwhm_min     = stddev(trial_results[5,*])
     fitsigmas[1].srcx     = stddev(trial_results[2,*])
     fitsigmas[1].srcy     = stddev(trial_results[3,*])
        
     endif

endif



if ~keyword_set(silent) then begin
  
PRINT
PRINT, 'COMPONENT    TYPE          FLUX       FWHM MAX    FWHM MIN      Angle     X loc      Y loc'
PRINT, '                         cts/s/keV     arcsec      arcsec        deg      arcsec     arcsec'
PRINT
nsrc = N_ELEMENTS(srcstr)
FOR n = 0, nsrc-1 DO BEGIN
  temp        = [ srcstr[n].srcflux,srcstr[n].srcfwhm_max,  srcstr[n].srcfwhm_min, srcstr[n].srcpa, srcstr[n].srcx, srcstr[n].srcy]
  PRINT, n+1, srcstr[n].srctype, temp, FORMAT="(I5, A13, F13.2, 1F13.1, F12.1, 2F11.1, F11.1, 2F12.1)"
  temp        = [ fitsigmas[n].srcflux,fitsigmas[n].srcfwhm_max, fitsigmas[n].srcfwhm_min, fitsigmas[n].srcpa, fitsigmas[n].srcx, fitsigmas[n].srcy]
  PRINT, ' ', '(std)', temp, FORMAT="(A7, A11, F13.2, 1F13.1, F12.1, 2F11.1, F11.1, 2F12.1)"
  PRINT, ' '
end

endif

return, {srcstr: srcstr, fitsigmas: fitsigmas}

end