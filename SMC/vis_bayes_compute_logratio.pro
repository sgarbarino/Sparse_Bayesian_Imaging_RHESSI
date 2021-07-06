;+
; NAME:
;   vis_bayes_compute_logratio
;
; PURPOSE:
;   Compute the ratio of two log likelihood functions.
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;
;-

function vis_bayes_compute_logratio, dato1, dato2, sigma_noise, visxyobs, esponente


  loglike1 = -(1.0/2.) * (norm((visxyobs - dato1)/sigma_noise))^2
  loglike2 = -(1.0/2.) * (norm((visxyobs - dato2)/sigma_noise))^2

  logratio = esponente*(loglike2 - loglike1)

  return, logratio
  end
