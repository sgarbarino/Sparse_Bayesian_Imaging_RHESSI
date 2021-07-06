
folder = '/Users/admin/Documents/PhD/Papers/STIX Imaging group/STIX demonstration/'

;code_folder = folder + 'STIX/code_Paolo_ampiezze/'
;add_path, code_folder

;data_folder = folder + 'idl_SSW_data/'
;cd, data_folder

add_path, folder + 'SMC/'

savefolder = folder + 'results/'
if savefolder eq 0 then FILE_MKDIR, savefolder

fname = savefolder + '.sav'

cerberus = 0
phase = 0

;;;; LOAD DATA

;; if November - flare peak, flare location from FSI

;fff=file_search(data_folder + 'all_L1/solo_L1A*.fits', /FOLD_CASE)
;;
;;;path_sci_file = fff(106)
;;;path_bkg_file = data_folder + 'all_L1/solo_L1A_stix-sci-xray-l1-1270153232_20201122T200008-20201122T213004_006928_V01.fits'
;;;time_range = ['19-Nov-2020 03:43:20.794','19-Nov-2020 03:46:00.794']
;;;position = [1000.,-300.]
;energy_range = [7,12]

;path_sci_file = fff(96)
;path_bkg_file = data_folder + 'stix_L1_fsi/solo_L1A_stix-sci-xray-l1-1270153232_20201122T200008-20201122T213004_006947_V01.fits'
;time_range = ['18-Nov-2020 15:29:30.000','18-Nov-2020 15:30:30.000']
;energy_range = [7,12]

;path_sci_file=fff(20)
;path_bkg_file='stix_L1_fsi/solo_L1A_stix-sci-xray-l1-1270153232_20201122T200008-20201122T213004_006947_V01.fits'
;time_range = ['16-Nov-2020 08:16:30','16-Nov-2020 08:17:30']
;energy_range = [7,12]

;path_sci_file = data_folder + 'all_L1/solo_L1A_stix-sci-xray-l1-1267819536_20201118T053624-20201118T055415_006558_V01.fits'
;path_bkg_file = data_folder + 'stix_L1_fsi/solo_L1A_stix-sci-xray-l1-1270153232_20201122T200008-20201122T213004_006947_V01.fits'
;time_range = ['18-Nov-2020 05:45:30','18-Nov-2020 05:46:15']
;time_range = ['18-Nov-2020 05:46:30','18-Nov-2020 05:46:31']

;energy_range = [7,12]
;position = [1000.,-350]

;day = 18
;month = 11
;year = 2020
;hour = 05

;minute1 = 46
;second1 = 30
;minute2 = 46
;second2 = 31

;minute1 = 45
;second1 = 30
;minute2 = 46
;second2 = 15


;;; if june - is NOT in all_L1 folder and reads as below
;path_sci_file = data_folder + 'solo_L1_stix-sci-xray-l1-1178428688_20200607T213708-20200607T215208_V01.fits'
;path_bkg_file = data_folder + 'solo_L1_stix-sci-xray-l1-1178451984_20200607T225959-20200607T235900_V01.fits'
time_range = ['7-Jun-2020 21:38:49', '7-Jun-2020 21:45:00']
energy_range = [6,10]
day = 7
month = 6
year = 2020
hour = 21
minute1 = 38
second1 = 49
minute2 = 45
second2 = 00


;data = stix_display_pixels_trans(path_sci_file,path_bkg_file, anytim(time_range), xy_flare = position, energy_range, silent=1)
;data = stix_display_pixels(path_sci_file,path_bkg_file, anytim(time_range),  energy_range, silent=1)

restore, folder + 'Data/data_7Jun2020_214149_214249_6_10keV.dat';'data_18Nov2020_054530_054615_21_70keV.dat'

;;; DETECTOR USED

subc_str = stx_construct_subcollimator()
ddet_idx = where(subc_str.label ne 'cfl' and $
  subc_str.label ne 'bkg' and $
  subc_str.label ne '1a' and $
  subc_str.label ne '1b' and $
  subc_str.label ne '1c' and $
  subc_str.label ne '2a' and $
  subc_str.label ne '2b' and $
  subc_str.label ne '2c')


;;; CONSTRUCT VISIBILITIES
dummy_vis = stx_construct_visibility(subc_str[ddet_idx], obsvis = data.amp_both[ddet_idx], sigamp = data.amp_both_error[ddet_idx], energy_range = energy_range);, time_range = time_range )

time_range_min = JULDAY(month, day, year, hour, minute1, second1)
time_range_max = JULDAY(month, day, year, hour, minute2, second2)
struct_add_field, dummy_vis, 'time_min', time_range_max
struct_add_field, dummy_vis, 'time_max', time_range_max

; Careful with obsvis as it is stored in a complex number, but we have not REAL/IMAG; instead we have AMP/NOTHING

;;; RECONSTRUCTION 

; Parameters:
pixel_size = 1.
fov=128.0

; Mean of the Poisson prior for the number of sources in a map
lamb = 1.;0.5

; Prior probabilities on the shape
; C = circle, E = ellipse, L = loop
;if cerberus eq 0 then begin
;  pC = 1
;  pE = 0
;  pL = 0
;endif else begin
;  pC = 1
;  pE = 0
;  pL = 0
;endelse
; Number of Monte Carlo samples
N_particles = 1000 ;[recommended: between 100 and 10,000 (the more the better)]

; run
pE = 1
pC = 0
pL = 0
vis_aux = dummy_vis
vis_aux.obsvis = abs(dummy_vis.obsvis)
asmc_sources = stx_vis_bayes(vis_aux, fname, PIXEL_SIZE=pixel_size, FOV=fov, LAM=lam, PE=pE, PC=pC, PLoop=pL, $
 N_PARTICLES=N_particles, AUTOSHAPE=0, PLOTHIST=1, phase = phase, cerberus = cerberus)


;
;smcsources[0].srcx = -smcsources[0].srcx
;smcsources[0].srcy = -smcsources[0].srcy
;smcsources[1].srcx = -smcsources[1].srcx
;smcsources[1].srcy = -smcsources[1].srcy
map_asmc = asmc_bayes_showmap_stix(fov,pixel_size, vis_aux[0].xyoffset, asmc_sources.smcsources, phase=phase, cerberus=cerberus)
;loadct, 5
;window, 1
;plot_map, map_asmc ;,title= anytim(time_title,/ecs)
;plot_map_stix, map_asmc, bla, savefolder + '/reconstruction.png'
;map2fits, map_asmc, savefolder +'/reconstruction.fits'
;
;
; ;;; CONSTRUCT VISIBILITIES for fit
;subc_str_fit = stx_construct_subcollimator()
;ddet_idx_fit = where(subc_str.label ne 'cfl' and  subc_str.label ne 'bkg')
;syserr = 0.05
;ampobs = data.amp_both[ddet_idx_fit]
;sigamp = data.amp_both_error[ddet_idx_fit]
;sigamp = SQRT(sigamp^2  + syserr^2 * ampobs^2)
;dummy_vis_fit = stx_construct_visibility(subc_str[ddet_idx_fit])
; 
;dim_vis = size(dummy_vis.u) 
;if cerberus eq 1 then n_par_fit = 6
;if cerberus eq 0 then n_par_fit = 4
;plot_stix_fit, map_asmc.data, data.amp_both[ddet_idx_fit], dummy_vis_fit.u, dummy_vis_fit.v, $
;   sigamp, dim_vis[1]-n_par_fit, dim_vis[1], [fov, fov], [pixel_size, pixel_size], savefolder + '/fit.png'
;
;vis_bayes_histograms_cerberus(par_s, wei_s)
;
;
;
;
;plot, par_s[0,0,*]
;par_s[0,0, 500]
;idx_0 = where((par_s[0,0,*] le -50) or (par_s[0,0,*] ge 50))
;dim_0 = size(idx_0)
;for i = 0, dim_0[1]-1 do par_s[0,0, idx_0[i]] = par_s[0,0, 500]
;idx_0 = where((par_s[0,1,*] le -50) or (par_s[0,1,*] ge 50))
;dim_0 = size(idx_0)
;for i = 0, dim_0[1]-1 do par_s[0,1, idx_0[i]] = par_s[0,1, 500]
;idx_0 = where((par_s[1,0,*] le -50) or (par_s[1,0,*] ge 50))
;dim_0 = size(idx_0)
;for i = 0, dim_0[1]-1 do par_s[1,0, idx_0[i]] = par_s[1,0, 500]
;idx_0 = where((par_s[1,1,*] le -50) or (par_s[1,1,*] ge 50))
;dim_0 = size(idx_0)
;for i = 0, dim_0[1]-1 do par_s[1,1, idx_0[i]] = par_s[1,1, 500]
;plot, par_s[0,0,*]
;
;;par_s = par_s[*,*,0:1500]
;;wei_s = wei_s[0:1500]
;;;; clustering a posteriori if we have two possible configurations
;; what if we try clustering position and just after flux
;
;tmp_par = par_s[*,[5],*] ; positions and flux (no fwhm)
;;par = [reform(tmp_par[0,*,*]), reform(tmp_par[1,*,*])]
;par = reform([(tmp_par[0,*,*]), (tmp_par[1,*,*])])
;
;weights = CLUST_WTS(par, N_CLUSTERS = 2)
;result = CLUSTER(par, weights, N_CLUSTERS = 2)
;
;iplot, par,  LINESTYLE = 6, SYM_INDEX = 4
;h2d = HIST_2D(reform(par[0,*]), reform(par[1,*]),bin1=.02, bin2=.02)
;;ct = COLORTABLE(0, /REVERSE)
;;IMAGE(h2d);, RGB_TABLE=ct)
;
;; first and second configurations (now with all parameters):
;par_complete = [reform(par_s[0,*,*]), reform(par_s[1,*,*])]
;config_0 = par_complete[*,where(reform(result) eq 0)]
;config_1 = par_complete[*,where(reform(result) eq 1)]
;config_2 = par_complete[*,where(reform(result) eq 2)]
;
;; weights of the 2 configurations:
;wei_0 = wei_s[where(reform(result) eq 0)]
;wei_1 = wei_s[where(reform(result) eq 1)]
;wei_2 = wei_s[where(reform(result) eq 2)]
; 
;; re-shape in standard [N_SOURCES,N_PAR,N_PARTICLES] shape
;dim0=size(where(reform(result) eq 0))
;dim1=size(where(reform(result) eq 1))
;dim2=size(where(reform(result) eq 2))
;
;par_s_clustered = make_array(2,7,1000);
;par_s_clustered[0,*,0:dim0[1]-1] = config_0[0:6,*]
;par_s_clustered[1,*,0:dim0[1]-1] = config_0[7:13,*]
;
;par_s_clustered[0,*,dim0[1]:dim0[1]+dim1[1]-1] = config_1[0:6,*]
;par_s_clustered[1,*,dim0[1]:dim0[1]+dim1[1]-1] = config_1[7:13,*]
;
;par_s_clustered[0,*,dim0[1]+dim1[1]:dim0[1]+dim1[1]+dim2[1]-1] = config_2[0:6,*]
;par_s_clustered[1,*,dim0[1]+dim1[1]:dim0[1]+dim1[1]+dim2[1]-1] = config_2[7:13,*]
;
;wei_s_clustered = [wei_0,wei_1, wei_2]
;wei_s_clustered = [wei_0,wei_1];, wei_2]
;
;
;save, config_0, config_1, config_2, par_s_clustered, wei_0, wei_1, wei_2, filename=savefolder + "file_for_histo.sav"
;save, config_0, config_1, par_s_clustered, wei_0, wei_1, $
;   filename=savefolder+"/file_for_histo.sav"
;
;
;
;
;vis_bayes_histograms_cerberus(par_s_clustered[0,*,0:dim0[1]-1], wei_0)
;vis_bayes_histograms_cerberus(par_s_clustered[1,*,0:dim0[1]-1], wei_0)
;
;vis_bayes_histograms_cerberus(par_s_clustered[0,*,dim0[1]:dim0[1]+dim1[1]-1], wei_1)
;vis_bayes_histograms_cerberus(par_s_clustered[1,*,dim0[1]:dim0[1]+dim1[1]-1], wei_1)
;
;vis_bayes_histograms_cerberus(par_s_clustered[0,*,dim0[1]+dim1[1]:dim0[1]+dim1[1]+dim2[1]-1], wei_2)
;vis_bayes_histograms_cerberus(par_s_clustered[1,*,dim0[1]+dim1[1]:dim0[1]+dim1[1]+dim2[1]-1], wei_2)
;
;
;; mean values and std values
;mean_values_0_0 = mean(config_0[0:6,*], dimension = 2)
;std_values_0_0 = stddev(config_0[0:6,*], dimension = 2)
;mean_values_0_1 = mean(config_0[7:13,*], dimension = 2)
;std_values_0_1 = stddev(config_0[7:13,*], dimension = 2)
;
;mean_values_1_0 = mean(config_1[0:6,*], dimension = 2)
;std_values_1_0 = stddev(config_1[0:6,*], dimension = 2)
;mean_values_1_1 = mean(config_1[7:13,*], dimension = 2)
;std_values_1_1 = stddev(config_1[7:13,*], dimension = 2)
;
;mean_values_2_0 = mean(config_2[0:6,*], dimension = 2)
;std_values_2_0 = stddev(config_2[0:6,*], dimension = 2)
;mean_values_2_1 = mean(config_2[7:13,*], dimension = 2)
;std_values_2_1 = stddev(config_2[7:13,*], dimension = 2)
;
;source_0 = reform(par_s[0,*,*])
;source_1 = reform(par_s[1,*,*])
;
;print, 'Original'
;print, mean(source_0, dimension = 2)
;print, stddev(source_0, dimension = 2)
;print, mean(source_1, dimension = 2)
;print, stddev(source_1, dimension = 2)
;
;print, 'Config 0'
;print, mean_values_0_0
;print, std_values_0_0
;print, mean_values_0_1
;print, std_values_0_1
;
;print, 'Config 1'
;print, mean_values_1_0
;print, std_values_1_0
;print, mean_values_1_1
;print, std_values_1_1
;
;print, 'Config 2'
;print, mean_values_2_0
;print, std_values_2_0
;print, mean_values_2_1
;print, std_values_2_1
;
;new_smcsources_0 = smcsources
;new_smcsources_0[0].srcx = -mean_values_0_0[0]
;new_smcsources_0[0].srcy = -mean_values_0_0[1]
;new_smcsources_0[0].SRCFWHM = mean_values_0_0[4]
;new_smcsources_0[0].SRCFLUX = mean_values_0_0[5]
;
;new_smcsources_0[1].srcx = -mean_values_0_1[0]
;new_smcsources_0[1].srcy = -mean_values_0_1[1]
;new_smcsources_0[1].SRCFWHM = mean_values_0_1[4]
;new_smcsources_0[1].SRCFLUX = mean_values_0_1[5]
;
;new_smcsources_1 = smcsources
;new_smcsources_1[0].srcx =  mean_values_1_0[0]
;new_smcsources_1[0].srcy =  mean_values_1_0[1]
;new_smcsources_1[0].SRCFWHM = mean_values_1_0[4]
;new_smcsources_1[0].SRCFLUX = mean_values_1_0[5]
;
;new_smcsources_1[1].srcx =  mean_values_1_1[0]
;new_smcsources_1[1].srcy =  mean_values_1_1[1]
;new_smcsources_1[1].SRCFWHM = mean_values_1_1[4]
;new_smcsources_1[1].SRCFLUX = mean_values_1_1[5]
;
;new_smcsources_2 = smcsources
;new_smcsources_2[0].srcx = - mean_values_2_0[0]
;new_smcsources_2[0].srcy = - mean_values_2_0[1]
;new_smcsources_2[0].SRCFWHM = mean_values_2_0[4]
;new_smcsources_2[0].SRCFLUX = mean_values_2_0[5]
;
;new_smcsources_2[1].srcx = - mean_values_2_1[0]
;new_smcsources_2[1].srcy = - mean_values_2_1[1]
;new_smcsources_2[1].SRCFWHM = mean_values_2_1[4]
;new_smcsources_2[1].SRCFLUX = mean_values_2_1[5]
;
;map_asmc_config = asmc_bayes_showmap_stix(fov,pixel_size, vis_aux[0].xyoffset, new_smcsources_1, phase=phase, cerberus=cerberus)
;loadct, 5
;window, 1
;plot_map, map_asmc_config ;,title= anytim(time_title,/ecs)
;
;plot_map_stix, map_asmc_config, bla, savefolder + '/reconstruction_config1.png'
;map2fits, map_asmc_config, savefolder +'/reconstruction_config1.fits'
;
;plot_stix_fit, map_asmc_config.data, data.amp_both[ddet_idx_fit], dummy_vis_fit.u, dummy_vis_fit.v, $
;  sigamp, dim_vis[1]-n_par_fit, dim_vis[1], [fov, fov], [pixel_size, pixel_size], savefolder + '/fit_config1.png'


u = dummy_vis.u
v = dummy_vis.v
n_free = 24. - 4.
det_used = 24.


stx_plot_fit, map_asmc.data, data.amp_both[ddet_idx], u, v, data.amp_both_error[ddet_idx], n_free, det_used, size(map_asmc.data, /dim), [pixel_size,pixel_size], wwindow=1
stx_plot_map, map_asmc, 'SMC', wwindow=0


end


