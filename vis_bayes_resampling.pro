;+
; NAME:
;   stx_vis_bayes_resampling
;
; PURPOSE:
;   Perform the systematic resampling if the ESS (Effective sample Size) is small enough
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;
;-

function vis_bayes_resampling, param, sample=sample, types=types, weights=weights, ess, new_sample=new_sample, log_weights=log_weights, Nsources=Nsources   ; N_particles   
      
       if (ess le double(param.N_particles)/2) then begin
         
             cum_sum_weights = cum_sum(weights)
             unif = make_array(param.N_particles,/double,value=double(1)/param.N_particles)
             unif[0] = double(0)
             z = randomu(seed)/param.N_particles 
             cum_sum_unif = z + cum_sum(unif)                                             
                                                   
             ccc=sample   
             Nsources_old=Nsources 
                                                
             for p = 0, param.N_particles-1 do begin
         
                 index = min(where(cum_sum_weights ge cum_sum_unif[p])) 
                 
                 if (index ge 0) then begin
                    sample[*,*,p] = ccc[*,*,index]
                    new_sample[*,*,p] = sample[*,*,p]
                    types[*,p] = types[*,index]
                    Nsources[p]=Nsources_old[index]
                 endif $
                 else if (index eq -1) then begin
                    print, 'weights too small!'
                 endif                                                 
             endfor      
          
          log_weights = make_array(param.N_particles,/double,value=alog(double(1)/param.N_particles)) 
          ess = vis_bayes_ess(weights=weights, log_weights)
                   
       endif
          
   return, ess
end
