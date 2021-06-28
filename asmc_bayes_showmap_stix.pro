function asmc_bayes_showmap_stix, fov, pixel, mapcenter, srcstr, phase=phase, cerberus=cerberus

mapsize=fov/pixel
xyoffset=mapcenter
srcstrin=srcstr

toofar    = 2.    ; points more than toofar*FWHM will be set to zero
;
; Define the map and its axes.
data    = FLTARR(mapsize,mapsize)
x       = (FINDGEN(mapsize)-mapsize/2.+0.5)*pixel+ xyoffset[0]
y       = (FINDGEN(mapsize)-mapsize/2.+0.5)*pixel+ xyoffset[1]
;
; Expand loops, if any, in input source structure
ok = WHERE(srcstr.srctype NE 'loop', nok)
IF nok GT 0 THEN srcstrin = srcstr[ok]
iloop = WHERE(srcstr.srctype EQ 'loop', nloop)
IF nloop GT 0 THEN BEGIN
  FOR i = 0, nloop-1 DO BEGIN
    strtemp = hsi_vis_fwdfit_makealoop(srcstr[iloop[i]])
    if nloop ne size(srcstr,/dimension) then begin
    IF N_ELEMENTS(srcstrin) GT 0 THEN srcstrin = [[srcstrin, strtemp]] ELSE srcstrin = strtemp
    endif else if nloop eq size(srcstr,/dimension) then begin
    IF N_ELEMENTS(srcstrin) GT nloop THEN srcstrin = [[srcstrin, strtemp]] ELSE srcstrin = strtemp
    endif
  ENDFOR
ENDIF


; Begin loop over source structure.
nsrc    = N_ELEMENTS(srcstrin)
FOR n = 0, nsrc-1 DO BEGIN
  fwhm    = srcstrin[n].srcfwhm
  fwhmeff2  = fwhm^2
  normfactor  = 4 * ALOG(2.) / !PI * srcstrin[n].srcflux / fwhmeff2   
  eccen   = srcstrin[n].eccen
  srcpa   = srcstrin[n].srcpa * !DTOR       ; Degrees to radians, radians E of N
  b       = fwhm * (1.-eccen^2)^0.25        ; Useful if source is elliptical
  
; Loop over y coordinates
  FOR ny = 0, mapsize-1 DO BEGIN
    
    if phase eq 1 then begin
      dy    = srcstrin[n].srcy- y[ny]
      dx    = srcstrin[n].srcx - x 
    endif else begin
      if cerberus eq 0 then begin 
        dy    =  - y[ny]     ; forcing plot in (0,0)
        dx    =  - x
      endif else begin
        dy    = srcstrin[n].srcy- y[ny]
        dx    = srcstrin[n].srcx - x
      endelse 
    endelse 
    
    dr2   = dx^2 + dy^2
    
    IF srcstrin[n].srctype EQ 'ellipse' THEN BEGIN
      pa      = ATAN(dy,dx)               ; radians W of N
      relpa   = pa - srcpa - !PI/2.           ; angle 'pixel to source' and ellipse axis
      fwhmeff2  = b^2 / (1 - (eccen*COS(relpa))^2)  ; nvis-element vector
    ENDIF
    term = 2.77259*dr2/fwhmeff2 < 20
    ok  = WHERE(term LT 20, nok)           ; avoid exponential underflows.
    IF nok EQ 0 THEN CONTINUE
    dflux = EXP(-term[ok])
    data[ok,ny] = data[ok,ny] + dflux * normfactor    ; add to one row for one source component
  ENDFOR
ENDFOR
  

map_particle=make_map(data,$
    id=' ', $ ;earth orbit
    xc=mapcenter[0],yc=mapcenter[1],$
    dx=pixel,dy=pixel,$
    xunits='arcsec', yunits='arcsec')
    
return, map_particle

end
