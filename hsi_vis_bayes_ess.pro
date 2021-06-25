;+
; NAME:
;   hsi_vis_bayes_ess
;
; PURPOSE:
;   Compute the ESS (Effective Sample Size) value 
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;
;-

function hsi_vis_bayes_ess, weights=weights, log_weights

             w = max(log_weights)
             logsum = w + alog(total(exp(log_weights-w)))
             weights = exp(log_weights-logsum)
             ess = double(1.0)/total(weights*weights)  
             
      return, ess
end
