
folder = '/Users/saragarbarino/Documents/WORK/MIDA/labs/astro/'

code_folder = folder + 'STIX/code_Paolo_ampiezze/'
add_path, code_folder

data_folder = folder + 'idl_SSW_data/'
cd, data_folder

savefolder = folder + '/STIX/SMC/results/18NovE3870_nocluster'
if savefolder eq 0 then FILE_MKDIR, savefolder

fname = savefolder + '.sav'

cerberus = 1
phase = 0

;;;; LOAD DATA

restore,  savefolder+'/18NovE3870_nocluster.sav'
restore,  savefolder+'/file_for_histo.sav'


path_sci_file = data_folder + 'all_L1/solo_L1A_stix-sci-xray-l1-1267819536_20201118T053624-20201118T055415_006558_V01.fits'
path_bkg_file = data_folder + 'stix_L1_fsi/solo_L1A_stix-sci-xray-l1-1270153232_20201122T200008-20201122T213004_006947_V01.fits'
time_range = ['18-Nov-2020 05:45:30','18-Nov-2020 05:46:15']

energy_range = [34,70]
position = [1000.,-350]

day = 18
month = 11
year = 2020
hour = 05


minute1 = 45
second1 = 30
minute2 = 46
second2 = 15


data = stix_display_pixels_trans(path_sci_file,path_bkg_file, anytim(time_range), xy_flare = position, energy_range, silent=1)


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


;;;; LOAD DATA


pixel_size = 1.
fov=128.0


vis_aux = dummy_vis
vis_aux.obsvis = abs(dummy_vis.obsvis)


;smcsources[0].srcx = -smcsources[0].srcx
;smcsources[0].srcy = -smcsources[0].srcy
;smcsources[1].srcx = -smcsources[1].srcx
;smcsources[1].srcy = -smcsources[1].srcy

;map_asmc = asmc_bayes_showmap_stix(fov,pixel_size, vis_aux[0].xyoffset, smcsources, phase=phase, cerberus=cerberus)
;loadct, 5
;window, 1
;plot_map, map_asmc ;,title= anytim(time_title,/ecs)

;;plot_map_stix, map_asmc, bla, savefolder + '/reconstruction.png'
;;map2fits, map_asmc, savefolder +'/reconstruction.fits'

;vis_bayes_histograms_cerberus(par_s, wei_s)


;;; CONSTRUCT VISIBILITIES for fit
subc_str_fit = stx_construct_subcollimator()
ddet_idx_fit = where(subc_str.label ne 'cfl' and  subc_str.label ne 'bkg')
syserr = 0.05
ampobs = data.amp_both[ddet_idx_fit]
sigamp = data.amp_both_error[ddet_idx_fit]
sigamp = SQRT(sigamp^2  + syserr^2 * ampobs^2)
dummy_vis_fit = stx_construct_visibility(subc_str[ddet_idx_fit])

dim_vis = size(dummy_vis.u)
if cerberus eq 1 then n_par_fit = 6
if cerberus eq 0 then n_par_fit = 4
;plot_stix_fit, map_asmc.data, data.amp_both[ddet_idx_fit], dummy_vis_fit.u, dummy_vis_fit.v, $
;  sigamp, dim_vis[1]-n_par_fit, dim_vis[1], [fov, fov], [pixel_size, pixel_size], savefolder + '/fit.png'


;;; clustering

;par_s = par_s[*,*,0:1600]
;wei_s = wei_s[0:1600]

plot, par_s[0,0,*]

num = 1000

idx_0 = where((par_s[0,0,*] le -10) or (par_s[0,0,*] ge 10))
dim_0 = size(idx_0)
par_s[0,0, idx_0] = par_s[0,0, num:num+dim_0[1]-1]

idx_0 = where((par_s[0,1,*] le -10) or (par_s[0,1,*] ge 10))
dim_0 = size(idx_0)
par_s[0,1, idx_0] = par_s[0,1,  num:num+dim_0[1]-1]

idx_0 = where((par_s[1,0,*] le -10) or (par_s[1,0,*] ge 10))
dim_0 = size(idx_0)
par_s[1,0, idx_0] = par_s[1,0, num:num+dim_0[1]-1]

idx_0 = where((par_s[1,1,*] le -10) or (par_s[1,1,*] ge 10))
dim_0 = size(idx_0)
par_s[1,1, idx_0] = par_s[1,1, num:num+dim_0[1]-1]
plot, par_s[0,0,*]


;;; clustering a posteriori if we have two possible configurations
; what if we try clustering position and just after flux

tmp_par = par_s[*,[5],*] ; positions and flux (no fwhm)
par = reform([(tmp_par[0,*,*]), (tmp_par[1,*,*])])

iplot, par,  LINESTYLE = 6, SYM_INDEX = 4
h2d = HIST_2D(reform(par[0,*]), reform(par[1,*]),bin1=.02, bin2=.02)
;ct = COLORTABLE(0, /REVERSE)
;IMAGE(h2d);, RGB_TABLE=ct)

weights = CLUST_WTS(par, N_CLUSTERS = 3)
result = CLUSTER(par, weights, N_CLUSTERS = 3)



; first and second configurations (now with all parameters):
par_complete = [reform(par_s[0,*,*]), reform(par_s[1,*,*])]
config_0 = par_complete[*,where(reform(result) eq 0)]
config_1 = par_complete[*,where(reform(result) eq 1)]
config_2 = par_complete[*,where(reform(result) eq 2)]

; weights of the 2 configurations:
wei_0 = wei_s[where(reform(result) eq 0)]
wei_1 = wei_s[where(reform(result) eq 1)]
wei_2 = wei_s[where(reform(result) eq 2)]

; re-shape in standard [N_SOURCES,N_PAR,N_PARTICLES] shape
dim0=size(where(reform(result) eq 0))
dim1=size(where(reform(result) eq 1))
dim2=size(where(reform(result) eq 2))

par_s_clustered = make_array(2,7,5000);
par_s_clustered[0,*,0:dim0[1]-1] = config_0[0:6,*]
par_s_clustered[1,*,0:dim0[1]-1] = config_0[7:13,*]

par_s_clustered[0,*,dim0[1]:dim0[1]+dim1[1]-1] = config_1[0:6,*]
par_s_clustered[1,*,dim0[1]:dim0[1]+dim1[1]-1] = config_1[7:13,*]

par_s_clustered[0,*,dim0[1]+dim1[1]:dim0[1]+dim1[1]+dim2[1]-1] = config_2[0:6,*]
par_s_clustered[1,*,dim0[1]+dim1[1]:dim0[1]+dim1[1]+dim2[1]-1] = config_2[7:13,*]

wei_s_clustered = [wei_0,wei_1, wei_2]
wei_s_clustered = [wei_0,wei_1];, wei_2]

save, config_0,  par_s_clustered, wei_0, filename=savefolder + "/file_for_histo.sav"

save, config_1,  par_s_clustered, wei_1, filename=savefolder+"/file_for_histo.sav"



;;;; plot of the histograms

vis_bayes_histograms_positions(par_s_clustered[0,*,0:dim0[1]-1], wei_0)
vis_bayes_histograms_positions(par_s_clustered[1,*,0:dim0[1]-1], wei_0)

vis_bayes_histograms_cerberus(par_s_clustered[0,*,0:dim0[1]-1], wei_0)
vis_bayes_histograms_cerberus(par_s_clustered[1,*,0:dim0[1]-1], wei_0)

vis_bayes_histograms_cerberus(par_s_clustered[0,*,dim0[1]:dim0[1]+dim1[1]-1], wei_1)
vis_bayes_histograms_cerberus(par_s_clustered[1,*,dim0[1]:dim0[1]+dim1[1]-1], wei_1)

vis_bayes_histograms_cerberus(par_s_clustered[0,*,dim0[1]+dim1[1]:dim0[1]+dim1[1]+dim2[1]-1], wei_2)
vis_bayes_histograms_cerberus(par_s_clustered[1,*,dim0[1]+dim1[1]:dim0[1]+dim1[1]+dim2[1]-1], wei_2)



;;;



; mean values and std values
mean_values_0_0 = mean(config_0[0:6,*], dimension = 2)
std_values_0_0 = stddev(config_0[0:6,*], dimension = 2)
mean_values_0_1 = mean(config_0[7:13,*], dimension = 2)
std_values_0_1 = stddev(config_0[7:13,*], dimension = 2)

mean_values_1_0 = mean(config_1[0:6,*], dimension = 2)
std_values_1_0 = stddev(config_1[0:6,*], dimension = 2)
mean_values_1_1 = mean(config_1[7:13,*], dimension = 2)
std_values_1_1 = stddev(config_1[7:13,*], dimension = 2)

mean_values_2_0 = mean(config_2[0:6,*], dimension = 2)
std_values_2_0 = stddev(config_2[0:6,*], dimension = 2)
mean_values_2_1 = mean(config_2[7:13,*], dimension = 2)
std_values_2_1 = stddev(config_2[7:13,*], dimension = 2)

source_0 = reform(par_s[0,*,*])
source_1 = reform(par_s[1,*,*])

print, 'Original'
print, mean(source_0, dimension = 2)
print, stddev(source_0, dimension = 2)
print, mean(source_1, dimension = 2)
print, stddev(source_1, dimension = 2)

print, 'Config 0'
print, mean_values_0_0
print, std_values_0_0
print, mean_values_0_1
print, std_values_0_1

print, 'Config 1'
print, mean_values_1_0
print, std_values_1_0
print, mean_values_1_1
print, std_values_1_1

print, 'Config 2'
print, mean_values_2_0
print, std_values_2_0
print, mean_values_2_1
print, std_values_2_1

new_smcsources_0 = smcsources
new_smcsources_0[0].srcx = - mean_values_0_0[0]
new_smcsources_0[0].srcy = - mean_values_0_0[1]
new_smcsources_0[0].SRCFWHM = mean_values_0_0[4]
new_smcsources_0[0].SRCFLUX = mean_values_0_0[5]

new_smcsources_0[1].srcx = - mean_values_0_1[0]
new_smcsources_0[1].srcy = - mean_values_0_1[1]
new_smcsources_0[1].SRCFWHM = mean_values_0_1[4]
new_smcsources_0[1].SRCFLUX = mean_values_0_1[5]

new_smcsources_1 = smcsources
new_smcsources_1[0].srcx = mean_values_1_0[0]
new_smcsources_1[0].srcy = mean_values_1_0[1]
new_smcsources_1[0].SRCFWHM = mean_values_1_0[4]
new_smcsources_1[0].SRCFLUX = mean_values_1_0[5]

new_smcsources_1[1].srcx = mean_values_1_1[0]
new_smcsources_1[1].srcy = mean_values_1_1[1]
new_smcsources_1[1].SRCFWHM = mean_values_1_1[4]
new_smcsources_1[1].SRCFLUX = mean_values_1_1[5]

new_smcsources_2 = smcsources
new_smcsources_2[0].srcx = - mean_values_2_0[0]
new_smcsources_2[0].srcy = - mean_values_2_0[1]
new_smcsources_2[0].SRCFWHM = mean_values_2_0[4]
new_smcsources_2[0].SRCFLUX = mean_values_2_0[5]

new_smcsources_2[1].srcx = - mean_values_2_1[0]
new_smcsources_2[1].srcy = - mean_values_2_1[1]
new_smcsources_2[1].SRCFWHM = mean_values_2_1[4]
new_smcsources_2[1].SRCFLUX = mean_values_2_1[5]



;; E6
;new_smcsources = smcsources
;new_smcsources[0].srcx = 7.9
;new_smcsources[0].srcy = -2.9
;new_smcsources[0].SRCFWHM = 17.2
;new_smcsources[0].SRCFLUX = 13.2
;
;new_smcsources[1].srcx = -7.9
;new_smcsources[1].srcy = 2.9
;new_smcsources[1].SRCFWHM = 8.1
;new_smcsources[1].SRCFLUX = 3.2

;; E5
;new_smcsources = smcsources
;new_smcsources[0].srcx = 7.3
;new_smcsources[0].srcy = -2.7
;new_smcsources[0].SRCFWHM = 17.0
;new_smcsources[0].SRCFLUX = 7.3
;
;new_smcsources[1].srcx = -7.3
;new_smcsources[1].srcy = 2.7
;new_smcsources[1].SRCFWHM = 10.7
;new_smcsources[1].SRCFLUX = 3.2

;; E4
;new_smcsources = smcsources
;new_smcsources[0].srcx = 7.3
;new_smcsources[0].srcy = -3.3
;new_smcsources[0].SRCFWHM = 16.9
;new_smcsources[0].SRCFLUX = 5.1
;
;new_smcsources[1].srcx = -7.3
;new_smcsources[1].srcy = 3.3
;new_smcsources[1].SRCFWHM = 8.3
;new_smcsources[1].SRCFLUX = 1.9

;; E3
;new_smcsources = smcsources
;new_smcsources[0].srcx = 7.7
;new_smcsources[0].srcy = -3.3
;new_smcsources[0].SRCFWHM = 14.3
;new_smcsources[0].SRCFLUX = 3.1
;
;new_smcsources[1].srcx = -7.7
;new_smcsources[1].srcy = 3.3
;new_smcsources[1].SRCFWHM = 10.5
;new_smcsources[1].SRCFLUX = 2.1

;; E2
;new_smcsources = smcsources
;new_smcsources[0].srcx = 8.0
;new_smcsources[0].srcy = -3.2
;new_smcsources[0].SRCFWHM = 12.0
;new_smcsources[0].SRCFLUX = 2.3
;
;new_smcsources[1].srcx = -8.0
;new_smcsources[1].srcy = 3.2
;new_smcsources[1].SRCFWHM = 9.1
;new_smcsources[1].SRCFLUX = 1.6

; E1
new_smcsources = smcsources
new_smcsources[0].srcx = 7.9
new_smcsources[0].srcy = -3.4
new_smcsources[0].SRCFWHM = 13.8
new_smcsources[0].SRCFLUX = 2.3

new_smcsources[1].srcx = -7.9
new_smcsources[1].srcy = 3.4
new_smcsources[1].SRCFWHM = 4.8
new_smcsources[1].SRCFLUX = 0.8

map_asmc_config = asmc_bayes_showmap_stix(fov,pixel_size, vis_aux[0].xyoffset, new_smcsources, phase=phase, cerberus=cerberus)
loadct, 5
window, 1
plot_map, map_asmc_config ;,title= anytim(time_title,/ecs)

plot_map_stix, map_asmc_config, bla, savefolder + '/reconstruction_config_final.png'
map2fits, map_asmc_config, savefolder +'/reconstruction_config_final.fits'

plot_stix_fit, map_asmc_config.data, data.amp_both[ddet_idx_fit], dummy_vis_fit.u, dummy_vis_fit.v, $
  sigamp, dim_vis[1]-n_par_fit, dim_vis[1], [fov, fov], [pixel_size, pixel_size], savefolder + '/fit_config_final.png'


end


