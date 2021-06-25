;+
; NAME:
;   hsi_vis_bayes_histograms
;
; PURPOSE:
;   Plot histograms of each variable (excluding loop angle), for each estimated source
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;-

function hsi_vis_bayes_histograms, par_s, wei_s

dimensions = size(par_s)

nums = dimensions(1)
numpart = dimensions(3)

FOR k = 0, nums-1 do begin
  FOR ipar = 0, 5 do begin
    minimum = min(par_s(k,ipar,0:numpart-1)) 
    maximum = max(par_s(k,ipar,0:numpart-1)) 
    punti = minimum + (maximum-minimum) * findgen(11)/11
    ; calcoliamo noi la pdf
    pdf = fltarr(11)
    FOR i = 0, numpart-1 DO BEGIN
      index = floor((par_s(k,ipar,i)-minimum)/(maximum-minimum)*10.999)
    ;  print, index
      pdf[index] = pdf[index] + wei_s(i)
    endfor
    
    if ipar EQ 0 then begin 
    phisto = BARPLOT(punti, pdf, LAYOUT = [3,2,1], $
      AXIS_STYLE=1, TITLE='X location')
    endif else begin 
      phisto = BARPLOT(punti, pdf, /CURRENT, LAYOUT=[3,2,ipar+1], $
      AXIS_STYLE=1)    
    endelse
    case ipar OF
      1: phisto.title = 'Y location'     
      2: phisto.title = 'Angle'          
      3: phisto.title = 'Eccentricity'      
      4: phisto.title = 'Width'
      5: phisto.title = 'Flux'
    else: 
    endcase
   endfor
endfor
end