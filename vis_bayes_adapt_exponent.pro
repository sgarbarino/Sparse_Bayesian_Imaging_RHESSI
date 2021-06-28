;+
; NAME:
;   vis_bayes_compute_d
;
; PURPOSE:
;   Compute d with the bisection method
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;-

function vis_bayes_compute_d, param, expon, ess_new, ess, d, da, db

  if (ess_new/ess le param.imin) then begin ; d decreases
    db = d
    d = max([(da + db)/2,  param.dmin])
  endif $
  else if (ess_new/ess ge  param.imax) then begin ; d increases
    da = d
    d = min([(da + db)/2,  param.dmax])
  endif

  if (expon + d gt 1) then begin
    d = max([1 - expon,  param.dmin])
  endif

  return, d
end
;----------------------------------------------------------------------------



;+
; NAME:
;   vis_bayes_adapt_exponent
;
; PURPOSE:
;   Compute the exponent gamma, the ESS and do the data stack
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;-

function vis_bayes_adapt_exponent, param, expon=expon, ess_new=ess_new, ess,  like_unit, new_log_weights=new_log_weights,  log_weights, print_esp=print_esp, weights=weights,Nsources,  sample, H=H
        
          da = param.dmin    
          db = param.dmax 
          d = param.dmax
          if (expon + d gt 1) then begin
             d = max([1 - expon, param.dmin])
          endif
      
          while (ess_new/ess ge param.imax) or (ess_new/ess le param.imin) do begin  
                                
                d = vis_bayes_compute_d(param, expon, ess_new, ess, d, da, db)  ;compute b by the bisection method
                
                incremento = (d/2.)*like_unit
                new_log_weights = log_weights+incremento
                
                ess_new = vis_bayes_ess(weights=weights, new_log_weights) 

                if (abs(d-param.dmin) le param.toll) or (abs(d-param.dmax) le param.toll) or (abs(expon+d-1) le param.toll) then begin
                   break 
                endif  
          end 
      
          expon = expon + d
          ess = ess_new
          log_weights = new_log_weights
          
          print_esp = [print_esp,expon] ;  
          
          
           ; Data stack:

           H = make_array(param.Nsamp,param.Nsamp)

           for p = 0, param.N_particles-1 do begin
             for k = 0, Nsources[p]-1 do begin
               ; find where the particle is located in the grid
               r = floor(param.Nsamp*(sample[k,1,p]-param.ysamp[0])/(param.ysamp[param.Nsamp-1]-param.ysamp[0]))
               s = floor(param.Nsamp*(sample[k,0,p]-param.xsamp[0])/(param.xsamp[param.Nsamp-1]-param.xsamp[0]))

               ; if the indexes are inside the image space
               if (r ge 0) and (r lt param.Nsamp) and (s ge 0) and (s lt param.Nsamp) then begin
                 H[s,r] = H[s,r] + weights[p]
               endif

             endfor
           endfor
           
          
    return,ess       
end
