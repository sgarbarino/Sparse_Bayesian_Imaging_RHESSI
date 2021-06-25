;+
; NAME:
;   hsi_vis_bayes_compute_la
;
; PURPOSE:
; Return the maximum value of the loop angle, given the eccentricity value
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;
;-

function hsi_vis_bayes_compute_la, ecc
         if (ecc le 0.5) then begin
         la=170
         endif else if (ecc ge 0.5) and (ecc le 0.55) then begin
         la=165
         endif else if (ecc ge 0.55) and (ecc le 0.65) then begin
         la=150
         ;endif else if (ecc ge 0.6) and (ecc le 0.65) then begin
         ;la=150
         endif else if (ecc ge 0.65) and (ecc le 0.7) then begin 
         la=140
         endif else if (ecc ge 0.7) and (ecc le 0.75) then begin
         la=130
         endif else if (ecc ge 0.75) and (ecc le 0.8) then begin
         la=120
         endif else if (ecc ge 0.8) and (ecc le 0.85) then begin
         la=105
         endif else if (ecc ge 0.85) and (ecc le 0.9) then begin
         la=85
         endif else if (ecc ge 0.9) and (ecc le 0.95) then begin 
         la=60
         endif else begin  
         la=30
         endelse
         return, la
end
