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

function vis_bayes_histograms_cerberus, par_s, wei_s, fw_0, fw_1, fl_0, fl_1

dimensions = size(par_s)

nums = dimensions(1)
numpart = dimensions(3)

binsize = 50

FOR k = 0, nums-1 do begin

  FOR ipar = 0, 5 do begin
    
    ; calcoliamo noi la pdf 
    pdf = fltarr(binsize)
    
;    if ipar eq 4 then begin  
;        if k eq 0 then begin 
;          fw = fw_0
;          dim = size(fw_0)
;          len = dim[1]
;        endif else begin
;          fw = fw_1
;          dim = size(fw_1)
;          len = dim[1]
;        endelse
;      
;        FOR i = 0, len-1 DO BEGIN
;          minimum = min(fw)
;          maximum = max(fw)
;          if maximum eq minimum then  maximum = minimum +1
;          punti = minimum + (maximum-minimum) * findgen(binsize)/binsize
;          index = floor((fw(i)-minimum)/(maximum-minimum)*(binsize-0.0001))
;          if i gt numpart-1 then begin
;            pdf[index] = pdf[index]+ wei_s(0)
;          endif else begin
;            pdf[index] = pdf[index]+ wei_s(i) 
;          endelse
;        endfor
    
;    endif else 
;    if ipar eq 5 then begin ; va nella riga sopra nel caso
;        if k eq 0 then begin
;          fl = fl_0
;          dim = size(fl_0)
;          len = dim[1]
;        endif else begin
;          fl = fl_1
;          dim = size(fl_1)
;          len = dim[1]
;        endelse
;        FOR i = 0, len-1 DO BEGIN
;          minimum = min(fl)
;          maximum = max(fl)
;          if maximum eq minimum then  maximum = minimum +1
;          punti = minimum + (maximum-minimum) * findgen(binsize)/binsize
;          index = floor((fl(i)-minimum)/(maximum-minimum)*(binsize-0.0001))          
;          if i gt numpart-1 then begin
;              pdf[index] = pdf[index]+ wei_s(0)
;          endif else begin
;            pdf[index] = pdf[index]+ wei_s(i) 
;          endelse
;        endfor    
;    endif else begin
      FOR i = 0, numpart-1 DO BEGIN
        minimum = min(par_s(k,ipar,0:numpart-1))
        maximum = max(par_s(k,ipar,0:numpart-1))
        if maximum eq minimum then  maximum = minimum +1
        punti = minimum + (maximum-minimum) * findgen(binsize)/binsize
        index = floor((par_s(k,ipar,i)-minimum)/(maximum-minimum)*(binsize-0.0001))
        pdf[index] = pdf[index] + wei_s(i)
      endfor

;    endelse
    
    if ipar EQ 0 then begin 
      phisto = BARPLOT(punti, pdf, LAYOUT = [2,2,1], $
      AXIS_STYLE=1, TITLE='X location');, XRANGE = [-10, -4]);[4, 10]
    endif else if ipar eq 1 then begin 
      phisto = BARPLOT(punti, pdf, LAYOUT=[2,2,2], $
      AXIS_STYLE=1,/CURRENT, TITLE = 'Y location');, XRANGE = [0, 8])    ;[-8, 0]
 
;    endif else if ipar eq 4 then begin
;      phisto = barplot(histogram(new_fw[k,*]), LAYOUT = [2,2,3], $
;      AXIS_STYLE = 1,/CURRENT, TITLE = 'Width')
;    endif else if ipar eq 5 then begin
;      phisto = barplot(histogram(new_fl[k,*]), LAYOUT = [2,2,4], $
;      AXIS_STYLE = 1,/CURRENT, TITLE = 'Flux')
;    endif
    
    endif else if ipar eq 4 then begin
      phisto = barplot(punti, pdf, LAYOUT = [2,2,3], $
      AXIS_STYLE = 1,/CURRENT, TITLE = 'Width');, XRANGE = [0, 50])
    endif else if ipar eq 5 then begin
      phisto = barplot(punti, pdf, LAYOUT = [2,2,4], $
      AXIS_STYLE = 1,/CURRENT, TITLE = 'Flux')
    endif
    

   endfor
endfor
end