;+
; NAME:
;   vis_bayes_histograms
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

function vis_bayes_histograms_nopos, par_s, wei_s, type, dimensions = dim

  dimensions = size(par_s)

  nums = dimensions(1)
  numpart = dimensions(3)
  
  
  case type of
    
    'circle': begin

      FOR k = 0, nums-1 do begin

        FOR ipar = 4, 5 do begin
          minimum = min(par_s(k,ipar,0:numpart-1))
          maximum = max(par_s(k,ipar,0:numpart-1))
          if maximum eq minimum then maximum = minimum +1
          punti = minimum + (maximum-minimum) * findgen(11)/11

          ; calcoliamo noi la pdf
          pdf = fltarr(11)
          FOR i = 0, numpart-1 DO BEGIN
            index = floor((par_s(k,ipar,i)-minimum)/(maximum-minimum)*10.9998)
            ;  print, index
            pdf[index] = pdf[index] + wei_s(i)
          endfor

          if ipar EQ 4 then begin
            phisto = BARPLOT(punti, pdf, LAYOUT = [2,1,1], $
              AXIS_STYLE=1, TITLE='FWHM', dimensions = dim)
          endif else begin
            phisto = BARPLOT(punti, pdf, /CURRENT, LAYOUT=[2,1,2], $
              AXIS_STYLE=1, TITLE='Flux', dimensions = dim)
          endelse
        endfor
      endfor

    end 

  
  'ellipse': begin
  
  FOR k = 0, nums-1 do begin  
    
    FOR ipar = 2, 5 do begin
      minimum = min(par_s(k,ipar,0:numpart-1))
      maximum = max(par_s(k,ipar,0:numpart-1))
      if maximum eq minimum then maximum = minimum +1
      punti = minimum + (maximum-minimum) * findgen(11)/11

      ; calcoliamo noi la pdf
      pdf = fltarr(11)
      FOR i = 0, numpart-1 DO BEGIN
        index = floor((par_s(k,ipar,i)-minimum)/(maximum-minimum)*10.9998)
        ;  print, index
        pdf[index] = pdf[index] + wei_s(i)
      endfor

      if ipar EQ 2 then begin
        phisto = BARPLOT(punti, pdf, LAYOUT = [2,2,1], $
          AXIS_STYLE=1, TITLE='Angle', dimensions = dim)
      endif else begin
        phisto = BARPLOT(punti, pdf, /CURRENT, LAYOUT=[2,2,ipar-1], $
          AXIS_STYLE=1, dimensions = dim)
      endelse
      case ipar OF
        3: phisto.title = 'Eccentricity'
        4: phisto.title = 'FWHM'
        5: phisto.title = 'Flux'
        else:
      endcase
    endfor
  endfor
  
  end
  
  
  endcase


  
end