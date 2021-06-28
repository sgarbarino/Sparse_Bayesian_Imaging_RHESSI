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

function vis_bayes_histograms_positions, par_s, wei_s

  dimensions = size(par_s)

  nums = dimensions(1) ; number of sources
  numpart = dimensions(3)

  binsize = 50


  FOR k = 0, nums-1 do begin
    print, k
    FOR ipar = 0, 1 do begin

      ; calcoliamo noi la pdf
      pdf = fltarr(binsize)

      FOR i = 0, numpart-1 DO BEGIN
        minimum = min(par_s(k,ipar,0:numpart-1))
        maximum = max(par_s(k,ipar,0:numpart-1))
        if maximum eq minimum then  maximum = minimum +1
        punti = minimum + (maximum-minimum) * findgen(binsize)/binsize
        index = floor((par_s(k,ipar,i)-minimum)/(maximum-minimum)*(binsize-0.0001))
        pdf[index] = pdf[index] + wei_s(i)
      endfor


      if ipar EQ 0 then begin
        phisto = BARPLOT(punti, pdf, LAYOUT = [1,2,1], $
          AXIS_STYLE=1, TITLE='X location', xtitle = 'arcsec', ytitle = 'pdf', FONT_SIZE=12, $
          XRANGE = [min([mean(par_s[k,ipar,*])+3,mean(par_s[k,ipar,*])-3]), max([mean(par_s[k,ipar,*])+3,mean(par_s[k,ipar,*])-3])])
          oplot, punti
      endif else if ipar eq 1 then begin
        phisto = BARPLOT(punti, pdf, LAYOUT=[1,2,2], $
          AXIS_STYLE=1,/CURRENT, TITLE = 'Y location', xtitle = 'arcsec', ytitle = 'pdf', FONT_SIZE=12, $
          XRANGE = [min([mean(par_s[k,ipar,*])+3,mean(par_s[k,ipar,*])-3]), max([mean(par_s[k,ipar,*])+3,mean(par_s[k,ipar,*])-3])]) 
      
      endif

    endfor
  endfor
end