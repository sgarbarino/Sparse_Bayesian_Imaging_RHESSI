
pro stx_plot_fit, im_map, ampobs, u, v, sigamp, n_free, det_used, imsize, pixel, wwindow=wwindow

default, wwindow, 0

F = vis_map2vis_matrix(u, v, imsize, pixel)
vispred = F # im_map[*]

ampobsmap = abs(vispred)

if det_used eq 30 then begin
  ind_d = [8L, 10L, 15L, $
           9L, 16L, 14L, $
           6L, 26L, 0L,  $
           22L, 4L, 20L, $
           5L, 27L, 1L,  $
           12L, 24L, 28L,$ 
           21L, 7L, 25L, $
           18L, 23L, 3L, $
           13L, 11L, 29L,$
           2L, 17L, 19L]
  xx = (findgen(30))/3. + 1.2
endif
                                          
if det_used eq 24 then begin

  ind_d = [26L, 6L, 0L,  $
          4L, 22L, 20L, $
          27L, 5L, 1L,  $
          12L, 28L, 24L, $
          21L, 25L, 7L, $
          18L, 3L, 23L, $
          29L, 11L, 13L,$
          19L, 17L, 2L]
  
  xx = (findgen(30))/3. + 1.2
  xx = xx[6:29]
endif
        
chi2 = total((ampobsmap[ind_d] - ampobs[ind_d])^2./sigamp[ind_d]^2.)/n_free

charsize = 1.2
leg_size = 1.4
thick = 1.6

linecolors, /quiet

title = 'CHI2: ' + num2str(chi2)
units = 'counts s!U-1!n keV!U-1!n'
xtitle = 'Label'

window, wwindow

plot, xx, ampobs[ind_d], /nodata, xrange=[1.,11.], /xst, xtickinterval=1, xminor=-1, $
  title=title, xtitle=xtitle, ytitle=units, yrange=yrange, charsize=charsize, thick=thick, _extra=_extra

; draw vertical dotted lines at each detector boundary
for i=1,10 do oplot, i+[0,0], !y.crange, linestyle=1


errplot, xx, (ampobs[ind_d]-sigamp[ind_d] > !y.crange[0]), ampobs[ind_d]+sigamp[ind_d] < !y.crange[1], $
  width=0, thick=thick, COLOR=7
oplot, xx, ampobs[ind_d], psym=7, thick=thick
oplot, xx, ampobsmap[ind_d], psym=4, col=2, thick=thick

leg_text = ['Observed', 'Error on Observed', 'From Image']
leg_color = [255, 7,2]
leg_style = [0, 0, 0]
leg_sym = [7, -3, 4]
ssw_legend, leg_text, psym=leg_sym, color=leg_color, linest=leg_style, box=0, charsize=leg_size, thick=thick, /left


end