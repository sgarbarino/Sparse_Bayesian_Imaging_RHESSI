;PRO main
; example code for using stx_vis_bayes with synthetic data generator

search_network, /enable
while !D.WINDOW ne -1 do wdelete
loadct, 5, /silent
COMMON uvdata, u,v, pa, mapcenter, alb_apply_index_orig
COMMON srcshape, shape


ph_src = stx_sim_flare(pixel_data=pixel_data, $
    src_shape = 'gaussian', $  
    src_xcen = 7, $
    src_ycen = -3, $
    src_flux = 1000., $
    src_fwhm_wd = 10., $
    src_fwhm_ht = 10., $
    src_phi = 0.)
   
pixel_data = stx_pixel_sums(pixel_data, 1)
subc_str = stx_construct_subcollimator()
vis = stx_visgen(pixel_data, subc_str)
l = 0.22
h = 0.92
M1 = l*h*4/(!pi^3.)*sin(!pi/4.)
vis.obsvis = vis.obsvis/(4. * M1)
vis.sigamp = vis.sigamp/(4. * M1)


ph_src_2 = stx_sim_flare(pixel_data=pixel_data_2, $
  src_shape = 'gaussian', $
  src_xcen = -7, $
  src_ycen = 3, $
  src_flux = 2000., $
  src_fwhm_wd = 10., $
  src_fwhm_ht = 10., $
  src_phi = 0.)

pixel_data_2 = stx_pixel_sums(pixel_data_2, 1)
subc_str_2 = stx_construct_subcollimator()
vis_2 = stx_visgen(pixel_data_2, subc_str_2)
vis_2.obsvis = vis_2.obsvis/(4. * M1)
vis_2.sigamp = vis_2.sigamp/(4. * M1)

vis_tot =  stx_visgen(pixel_data, subc_str)
vis_tot.obsvis = vis.obsvis + vis_2.obsvis
vis_tot.sigamp = vis.sigamp + vis_2.sigamp

;vis_tot = vis
;cerberus = 0
cerberus = 1


; RECONSTRUCTION

; Parameters:
pixel_size = 1.
fov=128.0

; Mean of the Poisson prior for the number of sources in a map
lamb = 1.;0.5

; Prior probabilities on the shape
; C = circle, E = ellipse, L = loop
pC = 1;1./2
pE = 0;1./4
pL = 0;1./4

; Number of Monte Carlo samples
N_particles = 5000 ;[recommended: between 100 and 10,000 (the more the better)]

; name of the file where the results will be saved
fname= '/Users/saragarbarino/Documents/WORK/MIDA/labs/astro/STIX/SMC/results/Bayes_.sav'


; full reconstruction with phases 
;phase = 1
;asmc_sources = stx_vis_bayes(vis_tot, fname, PIXEL_SIZE=pixel_size, FOV=fov, LAM=lam, PC=pC, PE=pE, PpL=pL, N_PARTICLES=N_particles, AUTOSHAPE=0, PLOtHIST=1, phase = phase)

; reconstruction from amplitude data
phase = 0

vis_aux = vis_tot
vis_aux.obsvis = abs(vis_tot.obsvis)

asmc_sources = stx_vis_bayes(vis_aux, fname, PIXEL_SIZE=pixel_size, FOV=fov, LAM=lam, PC=pC, PE=pE, PLoop=pL, $
  N_PARTICLES=N_particles, AUTOSHAPE=0, PLOtHIST=1, phase = phase, cerberus = cerberus)

map_asmc = asmc_bayes_showmap_stix(fov,pixel_size, vis_aux[0].xyoffset, asmc_sources.smcsources, phase=phase, cerberus=cerberus)
loadct, 5
window, 1
plot_map, map_asmc ;,title= anytim(time_title,/ecs)

  ;plot_map_stix, map_asmc, bla, savefolder + '/reconstruction.png'


END
