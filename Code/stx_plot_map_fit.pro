pro stx_plot_map_fit, amp_fwdfit_map, ampobs, u, v, sigamp, n_free, imsize, pixel, wwindow=wwindow, title=title

default, window, 0
default, title, ''

im_map = amp_fwdfit_map.data

F = vis_map2vis_matrix(u, v, imsize, pixel)
vispred = F # im_map[*]

ampobsmap = abs(vispred)


xx = (findgen(30))/3. + 1.2
xx = xx[6:29]

chi2 = total((ampobsmap - ampobs)^2./sigamp^2.)/n_free

charsize = 1.5
leg_size = 1.5
thick = 1.8
symsize = 1.8


units = 'counts s!U-1!n keV!U-1!n'
xtitle = 'Detector label'

window, wwindow, xsize=1200, ysize=500

loadct, 5

;set_viewport,0.1, 0.1, 0.4, 0.98
plot_map, amp_fwdfit_map, /cbar, title = title, multi=[1,2], $
          position=[0.1, 0.1, 0.4, 0.85], cb_position=[0.1, 0.95, 0.4, 0.98], charsize=charsize

linecolors, /quiet

set_viewport,0.5, 0.95, 0.1, 0.85


plot, xx, ampobs, /nodata, xrange=[1.,11.], /xst, xtickinterval=1, xminor=-1, $
  title='VISIBILITY AMPLITUDE FIT - CHI2: ' + trim(chi2, '(f12.2)'), $
  xtitle=xtitle, ytitle=units, yrange=yrange, charsize=charsize, thick=thick, /noe

; draw vertical dotted lines at each detector boundary
for i=1,10 do oplot, i+[0,0], !y.crange, linestyle=1


errplot, xx, (ampobs - sigamp > !y.crange[0]), (ampobs + sigamp < !y.crange[1]), $
  width=0, thick=thick, COLOR=7
oplot, xx, ampobs, psym=7, thick=thick, symsize=symsize
oplot, xx, ampobsmap, psym=4, col=2, thick=thick, symsize=symsize


leg_text = ['Observed', 'Error on Observed', 'From Image']
leg_color = [255, 7,2]
leg_style = [0, 0, 0]
leg_sym = [7, -3, 4]
ssw_legend, leg_text, psym=leg_sym, color=leg_color, linest=leg_style, box=0, charsize=leg_size, thick=thick, /left

!p.position = [0, 0, 0, 0]

end