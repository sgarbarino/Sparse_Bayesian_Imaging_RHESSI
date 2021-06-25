;+
; NAME:
;   hsi_vis_bayes_compute_vis
;
; PURPOSE:
;   Compute the visibilities given the particle
;
; HISTORY:
;   July 2018 Written by S. Lugaro, F. Sciacchitano and A. Sorrentino
;
; CONTACT:
;   sciacchitano [at] dima.unige.it
;   sorrentino [at] dima.unige.it
;
;-

function hsi_vis_bayes_compute_vis, typessrc, paramsrc, nsources

         COMMON uvdata, u,v, pa, mapcenter, alb_apply_index   
         COMMON srcshape, shape    
         jdum = dblarr(2*size(pa,/dimension))
         visxyobs = dblarr(N_elements(jdum)) 
                 
         if (nsources ge 1) then begin
            
            xx = paramsrc[*,0]
            yy = paramsrc[*,1]
            angle = paramsrc[*,2]
            ecc = paramsrc[*,3]
            fw = paramsrc[*,4]
            fl = paramsrc[*,5]
            loopang = paramsrc[*,6]
            
            truesources = replicate({hsi_vis_src_structure},nsources)

            for sorg=0, nsources-1 do begin
            
                truesources[sorg].srctype=shape[typessrc[sorg]]
                truesources[sorg].albedo_ratio = 0
                truesources[sorg].srcx = xx[sorg] 
                truesources[sorg].srcy = yy[sorg] 
                truesources[sorg].srcfwhm = fw[sorg]  
                truesources[sorg].srcflux = fl[sorg]
                truesources[sorg].srcpa = angle[sorg] 
                truesources[sorg].eccen = ecc[sorg] 
                truesources[sorg].loop_angle = loopang[sorg] 
                truesources[sorg].srcheight = 0   
                ;truesources[sorg].albedo_apply = 0   
                
                srcparm   = hsi_vis_fwdfit_structure2array(truesources[sorg], mapcenter)  
                visxyobs  += hsi_vis_fwdfit_func(jdum, srcparm)  
           endfor
         endif

     
     return, visxyobs 

end
