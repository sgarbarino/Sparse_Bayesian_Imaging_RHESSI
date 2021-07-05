; add the path of the folder that contains the code and the data
folder = '/Users/admin/Documents/PhD/Papers/STIX Imaging group/STIX phase calibration/'

add_path, folder + '/Codes for ps plots/'
add_path, folder + '/Calibration/'

xy_offset = [0., 0.]

ph_src = stx_sim_flare(pixel_data=pixel_data, $
  src_shape = 'point', $
  src_xcen = xy_offset[0], $
  src_ycen = xy_offset[1], $
  src_flux = 100000.)

pixel_data = stx_pixel_sums(pixel_data, 1)
subc_str = stx_construct_subcollimator()
vis = stx_visgen(pixel_data, subc_str)
u = vis.u
v = vis.v

phase = atan(imaginary(vis.obsvis), float(vis.obsvis)) * !radeg

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; PLOT BACKPROJECTION LINES

fov      = 64.     ;size of the FOV of the plots (in arcsec)
labels   = [2,3] ;labels of the collimators selected for a superimposed plot
                    ;(they must be between 1 and 10)
                    
cf_location = [xy_offset[0],xy_offset[1], 5, 5] ;The first two entries are the x,y coordinates of the flaring source estimated by the CFL
                                                  ;The second two entries are the related uncertainties
                                                  ;(ATTENTION: 50 and 50 are arbitrary values used for this demo! The correct ones should be used)

; Plot the backprojection lines (after phase correction)
stix_bp_lines, phase, u, v, xyoffset = xy_offset, fov = fov, labels = labels, ps_folder = folder + '/Figures';, $
               ;cf_location = cf_location


vis_ind = 8               
stx_vis_bpmap, vis[vis_ind], BP_FOV=64., PIXEL=0.25, MAP=map

loadct, 5

print_options, /port
popen, folder + 'Figures/bp_vis_11_1a.ps', units='cm', xsize=8, ysize=6

;set_viewport,0.15,0.5,0.15,0.5

plot_map, make_map(map, dx=0.25, dy=0.25), /notitle, position=[0.15,0.15,0.48,0.6], multi=[1,2]

set_viewport,0.6,0.95,0.15,0.6

xyoffset = xy_offset
res30=1./(2.*sqrt(u^2. + v^2.))
phase = phase / 180. * !pi ; phase in radians

; Discretization of the x-axis (used in the plots)
xx = (findgen(257)/256 - 0.5)*fov + xyoffset[0]

plot, xx, -u[vis_ind]/v[vis_ind]*xx + phase[vis_ind]/(2.*!pi*v[vis_ind]), $
  /xst, yrange = [xyoffset[1] - fov/2., xyoffset[1] + fov/2.], /yst, /nodata, /isotropic, title = this_title, charsize=charsize, $
  xtitle = 'X (arcsec)', ytitle = 'Y (arcsec)', thick=thick

; Computation of the number of lines needed for covering the FOV
if -u[vis_ind]/v[vis_ind] ge 0 then begin
  tmp = round(v[vis_ind] * (xyoffset[1] - fov/2. + u[vis_ind]/v[vis_ind] * (xyoffset[0] + fov/2)) $
    - phase[vis_ind]/(2.*!pi))
endif else begin
  tmp = round(v[vis_ind] * (xyoffset[1] - fov/2. + u[vis_ind]/v[vis_ind] * (xyoffset[0] - fov/2)) $
    - phase[vis_ind]/(2.*!pi))
endelse
n_k = round(fov * sqrt(2)/res30[vis_ind]/2) + 1
k = signum(v[vis_ind]) * findgen(n_k) + tmp


for i=0, n_k-1 do begin ; For loop for plotting the lines

  if keyword_set(ps_folder) then begin
    oplot, xx, -u[vis_ind]/v[vis_ind]*xx + phase[vis_ind]/(2.*!pi*v[vis_ind]) $
      + k[i]/v[vis_ind], thick=thick
  endif else begin
    oplot, xx, -u[vis_ind]/v[vis_ind]*xx + phase[vis_ind]/(2.*!pi*v[vis_ind]) $
      + k[i]/v[vis_ind], thick=thick
  endelse

endfor


pclose

!p.position = [0, 0, 0, 0]         
               

end