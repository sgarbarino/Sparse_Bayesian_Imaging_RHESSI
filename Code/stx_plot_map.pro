pro stx_plot_map, im_map, title, wwindow=wwindow

  default, wwindow, 0

  xsize=820
  ysize=850
  csize = 2.

  loadct, 5, /silent
  hsi_linecolors

  window, wwindow

  plot_map, im_map, /cbar, /no_timestamp, /equal, $
    bottom=13, title=title;, cb_position=[.16, .95, .95, .97], charsize=csize, position=[.16, .1, .95, .85]


end 