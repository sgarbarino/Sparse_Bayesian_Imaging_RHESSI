
;;;;;;;;;;;;;;;; Add the path of the folder containing the code
folder      = ''
add_path, folder + 'PSO/'
add_path, folder + 'Code/'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                                          ;
;                                        READ DATA                                         ;
;                                                                                          ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;; Restore the data - November 18, 2020

data_folder = folder + 'Data June 7 2020/'
cd, data_folder

sci_file     = data_folder + 'solo_L1_stix-sci-xray-l1-1178428688_20200607T213708-20200607T215208_V01.fits'
bkg_file     = data_folder + 'solo_L1_stix-sci-xray-l1-1178451984_20200607T225959-20200607T235900_V01.fits'
time_range   = anytim(['7-Jun-2020 21:41:49', '7-Jun-2020 21:42:49'])
energy_range = [6,10]
xy_flare     = [-1600., -800.]

data = stix_display_pixels_trans(sci_file, bkg_file, time_range, energy_range, xy_flare = [-1600., -800.], silent=1)



;;;;;;;;;;;;;;;;;; Subcollimators used for the reconstruction

;           '3a' '3b' '3c' '4a' '4b' '4c' '5a' '5b' '5c' '6a' '6b' '6c'
subc_3_10 = [7,   29,  1,   25,  5,   23,  6,   30,  2,   15,  27,  31, $
;           '7a' '7b' '7c' '8a' '8b' '8c' '9a' '9b' '9c' '10a' '10b' '10c'
             24,  8,   28,  21,  26,  4,   16,  14,  32,  3,   20,  22]-1


;;;;;;;;;;;;;;;; Visibility amplitudes (sum top and bottom pixels)

ampobs = data.amp_both[subc_3_10]


;;;;;;;;;;;;;;;; Error on the visibility amplitudes. 5% of systematic error is added

sigamp = data.amp_both_error[subc_3_10]
syserr = 0.05
sigamp = SQRT(sigamp^2  + syserr^2 * ampobs^2)


;;;;;;;;;;;;;;;; (u,v)-points

uv = stx_uv_points()
u = uv.u
v = uv.v
u = u[subc_3_10]
v = v[subc_3_10]


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                                          ;
;                              PARTICLE SWARM OPTIMIZATION                                 ;
;                                                                                          ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;; Type of the parametric shape used in the forward fit:
;
; - 'circle' : Gaussian circular source
; - 'ellipse': Gaussian elliptical source
; - 'multi'  : double Gaussian circular source

type = 'ellipse'


;;;;;;;;;;;;;;;; Degrees of freedom

n_amp = n_elements(ampobs)
if type eq 'circle'  then n_free = n_amp - 2.
if type eq 'ellipse' then n_free = n_amp - 4.
if type eq 'multi'   then n_free = n_amp - 6.


;;;;;;;;;;;;;;;; Set to 1 for the computing the parameter uncertainties (confidence strip approach)

uncertainty = 1


;;;;;;;;;;;;;;;; FORWARD FIT WITH PSO

result = stx_amp_pso(type, ampobs, sigamp, u, v, n_free, uncertainty = uncertainty)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                                          ;
;                                    PLOT OF THE MAP                                       ;
;                                                                                          ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


imsize = [128, 128] ; number of pixels
pixel  = [1., 1.]   ; reconstructed map


pso_map = amp_fwdfit_source2map(result.srcstr, pixel = pixel, imsize = imsize)
stx_plot_map_fit, pso_map, ampobs, u, v, sigamp, n_free, imsize, pixel, wwindow=1, title='PSO reconstruction'



end