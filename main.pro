
; example code for using hsi_vis_bayes:

search_network, /enable
while !D.WINDOW ne -1 do wdelete 
loadct, 5, /silent 
COMMON uvdata, u,v, pa, mapcenter, alb_apply_index_orig
COMMON srcshape, shape 


;mapcenter= [918.0, 260.0];[914.0, 255.0]
time_interval=['2002-Feb-20 11:06:10', '2002-Feb-20 11:06:24' ];['2002-Feb-20 11:06:02', '2002-Feb-20 11:06:34' ]
epsmin=20;22
epsmax=30;26
detectors=[0,0,1,1,1,1,1,1,1]

ov=hsi_visibility()
ov->set, DET_INDEX_MASK=detectors,$
  xyoffset=mapcenter, $
  im_time_interval=time_interval,$
  phz_n_roll_bins_min=6, $
  phz_n_roll_bins_max=64,$
  im_energy_binning=[epsmin,epsmax], $
  vis_conjugate=0, $ 
  vis_normalize=1, $ 
  vis_edit=1, $ 
  use_phz_stacker= 1L, $
  modpat_skip= 4

vis=ov->getdata()


; Parameters:
pixel_size = 1.
fov=64.0  ;

; Mean of the Poisson prior for the number of sources in a map
lamb = 1.;0.5

; Prior probabilities on the shape 
; C = circle, E = ellipse, L = loop
pC = 1./2
pE = 1./4
pL = 1./4

; Number of Monte Carlo samples 
N_particles = 50 ;[recommended: between 100 and 10,000 (the more the better)]

; name of the file where the results will be saved

fname= 'Bayes_.sav'
asmc_sources= hsi_vis_bayes(vis, fname, PIXEL_SIZE=pixel_size, FOV=fov, LAM=lam, PC=pC, PE=pE, PpL=pL, N_PARTICLES=N_particles, AUTOSHAPE=0, PLOTHIST=1)

end
