FUNCTION stix_display_pixels,sci_file,bkg_file,tr_flare,er_flare,stop=stop, silent=silent 

default, silent, 0

;extracts pixel data and displays top and bottom pixel sets 

;input
; sci_file: L1 data with flare 
; bkg_file: L1 data with background (long time average)
; tr_flare: time range of flare to sum over
; er_flare: energy range to sum over
;
;B6 flare on June 7, 2020
;  sci_file='solo_L1_stix-sci-xray-l1-1178428688_20200607T213708-20200607T215208_V01.fits'
;  ;background is taken after the flare (1 hour integrated spectrum for each det and pixel)
;  bkg_file='/Users/samkrucker/Dropbox/sswidl/stix_flares/s20200607_22_B6/solo_L1_stix-sci-xray-l1-1178451984_20200607T225959-20200607T235900_V01.fits'
;  out=stix_display_pixels(sci_file,bkg_file,anytim(['7-Jun-2020 21:38:49','7-Jun-2020 21:45:00']),[6,7])

;read data
;note: deadtime correction assumes 12.5d-6 and ignores double triggers (that is ok for the June 7 flare, but not necessarily for Nov 2020 flares)
  stix_l1_read,sci_file,spec_all=spec_all,dspec_all=dspec_all,rspec_all=rspec_all,drspec_all=drspec_all,time_spec=time_spec,e0=e0,e1=e1,ee=ee,ddee=ddee

  ;data is now in the following variables
  ;SPEC_ALL        FLOAT     = Array[45, 32, 12, 32]  = [time, det, pix, energy]
  ;DSPEC_ALL       FLOAT     = Array[45, 32, 12, 32]
  ;RSPEC_ALL       FLOAT     = Array[45, 32, 12, 32]
  ;DRSPEC_ALL      FLOAT     = Array[45, 32, 12, 32]
  ;TIME_SPEC1      DOUBLE    = Array[45]
  ;E0              DOUBLE    = Array[32]
  ;E1              DOUBLE    = Array[32]
  ;EE              DOUBLE    = Array[32]

;background is taken after the flare
  stix_l1_read,bkg_file,spec_all=spec_all_bkg,dspec_all=dspec_all_bkg,rspec_all=rspec_all_bkg,drspec_all=drspec_all_bkg,time_spec=time_spec_bkg,time_dur=time_dur_bkg,e0=e0,e1=e1,ee=ee,ddee=ddee
  ;same format
  ;SPEC_ALL_BKG   FLOAT     = Array[1, 32, 12, 32] = [time, det , pix, energy]


  ;subcollimator indices
  g10=[3,20,22]-1
  l10=['10a','10b','10c']
  g09=[16,14,32]-1
  l09=['9a','9b','9c']
  g08=[21,26,4]-1
  l08=['8a','8b','8c']
  g07=[24,8,28]-1
  l07=['7a','7b','7c']
  g06=[15,27,31]-1
  l06=['6a','6b','6c']
  g05=[6,30,2]-1
  l05=['5a','5b','5c']
  g04=[25,5,23]-1
  l04=['4a','4b','4c']
  g03=[7,29,1]-1
  l03=['3a','3b','3c']
  g02=[12,19,17]-1
  l02=['2a','2b','2c']
  g01=[11,13,18]-1
  l01=['1a','1b','1c']
  res32=fltarr(32)
  res32(g10)=178.6
  res32(g09)=124.9
  res32(g08)=87.3
  res32(g07)=61.0
  res32(g06)=42.7
  res32(g05)=29.8
  res32(g04)=20.9
  res32(g03)=14.6
  res32(g02)=10.2
  res32(g01)=7.1
  o32=intarr(32)
  o32(g10)=[150,90,30]
  o32(g09)=[170,110,50]
  o32(g08)=[10,130,70]
  o32(g07)=[30,150,90]
  o32(g06)=[50,170,110]
  o32(g05)=[70,10,130]
  o32(g04)=[90,30,150]
  o32(g03)=[110,50,170]
  o32(g02)=[130,70,10]
  o32(g01)=[150,90,30]
  

  g03_10=[g03,g04,g05,g06,g07,g08,g09,g10]
  g01_10=[g01,g02,g03,g04,g05,g06,g07,g08,g09,g10]
  g_coarse=[g06,g07,g08,g09,g10]
  g_fine=[g01,g02,g03,g04,g05]
  g_plot=[g10,g05,g09,g04,g08,g03,g07,g02,g06,g01]
  l_plot=[l10,l05,l09,l04,l08,l03,l07,l02,l06,l01]
  r_plot=res32(g_plot)
  o_plot=o32(g_plot)
  ;[res32(g10),res32(g05),res09,res04,res08,res03,res07,res02,res06,res01]
  ;o_plot=[o10,o05,o09,o04,o08,o03,o07,o02,o06,o01]
  

  ; list of indices for time range
  tr_flare=anytim(tr_flare)
  tlist=where( (time_spec ge tr_flare(0)) AND (time_spec le tr_flare(1)) )
  ;same for energy range
  elist=where( (ee ge er_flare(0)) AND (ee le er_flare(1)) )
  ;title for plot
  range_title=strmid(anytim(tr_flare(0),/vms),0,11)+' '+strmid(anytim(tr_flare(0),/vms),12,8)+'-'+strmid(anytim(tr_flare(1),/vms),12,8)+'UT & '+strtrim(fix(er_flare(0)),2)+'-'+strtrim(fix(er_flare(1)),2)+' keV'


  ;ELUT
  ;correction should be done depending on the spectral shape
  ;here the lazy way: simply correct for actual bin size
  ;f_elut='/Users/samkrucker/Dropbox/sswidl/stix_ql/ELUT/elut_table_20200519.csv'
  f_elut='elut_table_20200519.csv'
  stx_read_elut, gain, offset, str4096, elut_filename = f_elut, scale1024=0, ekev_a = ekev
  ;ekev=[energy edges, pixel, det], only includes science energy bins, not 0 and last
  this_bin_low=reform(ekev(0:29,*,*))
  this_bin_high=reform(ekev(1:30,*,*))
  this_bin_size=this_bin_high-this_bin_low
  ;this_bin_size=[energy,pixel,det]
  ;THIS_BIN_SIZE   DOUBLE    = Array[30, 12, 32]
  ;rspec_all=[time, det, pix, energy]
  ;RSPEC_ALL       FLOAT     = Array[45, 32, 12, 32]
  this_bin_size_switch=fltarr(32,12,30)
  for i=0,29 do for j=0,11 do this_bin_size_switch(*,j,i)=this_bin_size(i,j,*)
  rspec_all_k=rspec_all
  drspec_all_k=drspec_all
  rspec_all_bkg_k=rspec_all_bkg
  for i=0,n_elements(time_spec)-1 do begin
    rspec_all_k(i,*,*,1:30)=rspec_all(i,*,*,1:30)/this_bin_size_switch
    drspec_all_k(i,*,*,1:30)=drspec_all(i,*,*,1:30)/this_bin_size_switch
  endfor
  rspec_all_bkg_k(0,*,*,1:30)=rspec_all_bkg(0,*,*,1:30)/this_bin_size_switch


  ;get rates  averaged over selected time and energy range: output = rate [det,pix]

  if n_elements(elist) eq 1 then begin
    this_r=average(rspec_all(tlist,*,*,elist),1)
    this_dr=sqrt(total(drspec_all(tlist,*,*,elist)^2,1))/(n_elements(tlist)*n_elements(elist))
    this_bkg=average(rspec_all_bkg(*,*,*,elist),1)
  endif else begin
    this_r=average(average(rspec_all(tlist,*,*,elist),4),1)
    this_dr=sqrt(total(total(drspec_all(tlist,*,*,elist)^2,4),1))/(n_elements(tlist)*n_elements(elist))
    this_bkg=average(average(rspec_all_bkg(*,*,*,elist),4),1)
  endelse

  ;plot to display what range was selected
  this_e=3
  if ~silent then begin
  loadct,5
  window,2,xsize=520,ysize=400
  if n_elements(elist) ne 1 then begin
    utplot,time_spec,total(total(total(rspec_all(*,g01_10,*,elist),4),3),2),psym=10,ytitle='STIX count rate [s!U-1!N]',title=range_title
    outplot,time_spec,time_spec*0.+total(rspec_all_bkg(0,g01_10,*,elist)),color=166
    outplot,time_spec(tlist),total(total(total(rspec_all(tlist,g01_10,*,elist),4),3),2),psym=10,color=122
  endif else begin 
    utplot,time_spec,total(total(rspec_all(*,g01_10,*,elist),3),2),psym=10,ytitle='STIX count rate [s!U-1!N]',title=range_title
    outplot,time_spec,time_spec*0.+total(rspec_all_bkg(0,g01_10,*,elist)),color=166
    outplot,time_spec(tlist),total(total(rspec_all(tlist,g01_10,*,elist),3),2),psym=10,color=122
  endelse
  endif
  ;subtract background
  this_rb=this_r-this_bkg
  ;assume background subtraction is perfect
  this_drb=this_dr


  ;display observed moire pattern for each detector
  if ~silent then begin
  window,0,xsize=900,ysize=800
  xmargin=0.08
  ymargin_top=0.12
  ymargin_bot=0.02
  xleer=0.02
  yleer=0.03
  xim=(1-2*xmargin-xleer)/6.
  yim=(1-ymargin_top-ymargin_bot-4*yleer)/5.
  c_top=122
  c_bot=44
  chs=1.0
  for i=0,29 do begin
    this_resolution=i/3
    this_row=i/6
    this_i=i-6*this_row
    if this_i ge 3 then this_space=xleer else this_space=0
    set_viewport,xmargin+this_i*xim+this_space,xmargin+(this_i+1)*xim+this_space,1-ymargin_top-(this_row+1)*yim-(this_row-1)*yleer,1-ymargin_top-this_row*yim-(this_row-1)*yleer
    this_title=l_plot(i)+'('+strtrim(fix(g_plot(i)),2)+');'+strtrim(fix(res32(g_plot(i))),2)+'";'+strtrim(fix(o32(g_plot(i))),2)+'!Uo!N'
    if i ne 24 then begin
      plot,[0.5,1.5,2.5,3.5],this_rb(g_plot(i),0:3),xtitle=' ',ytitle=' ',psym=-1,charsi=chs,yrange=[0,max(this_r(g_plot,0:7))],noe=i,xtickname=replicate(' ',9),ytickname=replicate(' ',9),xticks=8,xminor=1,xticklen=1d-22,title=this_title
      oplot,[0.5,1.5,2.5,3.5],this_rb(g_plot(i),0:3),psym=-1,color=c_top
      errplot,[0.5,1.5,2.5,3.5],(this_rb(g_plot(i),0:3)-this_dr(g_plot(i),0:3)),(this_rb(g_plot(i),0:3)+this_dr(g_plot(i),0:3)),thick=th3,color=c_top
      oplot,[0.5,1.5,2.5,3.5],this_rb(g_plot(i),4:7),psym=-1,color=c_bot
      errplot,[0.5,1.5,2.5,3.5],(this_rb(g_plot(i),4:7)-this_dr(g_plot(i),4:7)),(this_rb(g_plot(i),4:7)+this_dr(g_plot(i),4:7)),thick=th3,color=c_bot
    endif else begin
      plot,[0.5,1.5,2.5,3.5],this_rb(g_plot(i),0:3),xtitle=' ',ytitle='cts/s/keV',psym=-1,charsi=chs,yrange=[0,max(this_r(g_plot,*))],noe=i,xtickname=[' ','A',' ','B',' ','C',' ','D',' '],xticks=8,xminor=1,xticklen=1d-22,title=this_title
      oplot,[0.5,1.5,2.5,3.5],this_rb(g_plot(i),0:3),psym=-1,color=c_top  
      errplot,[0.5,1.5,2.5,3.5],(this_rb(g_plot(i),0:3)-this_dr(g_plot(i),0:3)),(this_rb(g_plot(i),0:3)+this_dr(g_plot(i),0:3)),thick=th3,color=c_top
      oplot,[0.5,1.5,2.5,3.5],this_rb(g_plot(i),4:7),psym=-1,color=c_bot
      errplot,[0.5,1.5,2.5,3.5],(this_rb(g_plot(i),4:7)-this_dr(g_plot(i),4:7)),(this_rb(g_plot(i),4:7)+this_dr(g_plot(i),4:7)),thick=th3,color=c_bot
    endelse
  endfor
  xyouts,0.5,1-ymargin_top/2.5,range_title,/normal,chars=1.6,ali=0.5  
  endif
  
  ;make amplitude
  y_top=this_rb(*,2)-this_rb(*,0)
  x_top=this_rb(*,3)-this_rb(*,1)
  dy_top=sqrt( this_drb(*,2)^2+this_drb(*,0)^2 )
  dx_top=sqrt( this_drb(*,3)^2+this_drb(*,1)^2 )
  r2_top=x_top^2+y_top^2
  y_bot=this_rb(*,6)-this_rb(*,4)
  x_bot=this_rb(*,7)-this_rb(*,5)
  dy_bot=sqrt( this_drb(*,6)^2+this_drb(*,4)^2 )
  dx_bot=sqrt( this_drb(*,7)^2+this_drb(*,5)^2 )

  r2_bot=x_bot^2+y_bot^2
  amp_top=sqrt(r2_top)
  amp_bot=sqrt(r2_bot)
  ;error on amplitudes
  damp_top=sqrt( ((x_top)/amp_top*dx_top)^2+((y_top)/amp_top*dy_top)^2 )
  damp_bot=sqrt( ((x_bot)/amp_bot*dx_bot)^2+((y_bot)/amp_bot*dy_bot)^2 )
  
  ;top and bottom summed
;  x_both=(x_top+x_bot)/2.
;  y_both=(y_top+y_bot)/2.
;  dx_both=(dx_top+dx_bot)/2./sqrt(2.)
;  dy_both=(dy_top+dy_bot)/2./sqrt(2.)
;  r2_both=x_both^2+y_both^2
;  amp_both=sqrt(r2_both)
;  damp_both=sqrt( ((x_both)/amp_both*dx_both)^2+((y_both)/amp_both*dy_both)^2 )
  
  y_both = (this_rb(*,2)+this_rb(*,6)) - (this_rb(*,0)+this_rb(*,4))
  x_both = (this_rb(*,3)+this_rb(*,7)) - (this_rb(*,1)+this_rb(*,5))
  dy_both = sqrt( (this_drb(*,2)+this_drb(*,6))^2 + (this_drb(*,0)+this_drb(*,4))^2 )
  dx_both = sqrt( (this_drb(*,3)+this_drb(*,7))^2 + (this_drb(*,1)+this_drb(*,5))^2 )
  r2_both = x_both^2+y_both^2
  amp_both = sqrt(r2_both)

  ;error on amplitudes
  damp_both = sqrt( ((x_both)/amp_both*dx_both)^2+((y_both)/amp_both*dy_both)^2 )
  
  ;make plot of amplitude vs resolution
  if ~silent then begin
  window,3,xsize=520,ysize=400,xpos=0,ypos=40
  clearplot
  ;shift display for bottom pixel to avoid overlap
  this_ff=1.1

  plot_oo,(1./res32)^2,amp_top,psym=1,xtitle='1/resolution^2',ytitle='amplitudes',yrange=[min([amp_top,amp_bot]),max([amp_top+damp_top,amp_bot+damp_bot])],yst=1,/nodata,title='RED: top; BLUE: bottom',xrange=[2d-5,1d-1],xsty=1
  oplot,(1./res32)^2,amp_top,psym=1,color=c_top,symsize=this_ss,thick=th3
  err_plot,(1./res32)^2,(amp_top-damp_top)>0.001,amp_top+damp_top,color=c_top,thick=th3,width=1d-13
  oplot,(1./res32)^2*this_ff,amp_bot,psym=1,color=c_bot,symsize=this_ss,thick=th3
  err_plot,(1./res32)^2*this_ff,(amp_bot-damp_bot)>0.001,amp_bot+damp_bot,color=c_bot,thick=th3,width=1d-13
  ;color coded according to orientation
  c32=o32/180.*255.
  ;for i=0,29 do plots,(1./res32(g01_10(i)))^2,amp_top(g01_10(i)),psym=1,color=c32(g01_10(i))
  
  ;make amplutide vs resolution plot for summed, color coded by orientation (3 colors for each set)
  window,5,xsize=520,ysize=400,xpos=0,ypos=40
  clearplot
  ;color code orientation in each triplette: smallest orientation: red, center: white, largest: green)
  cc33=o32*0.+255
  cc33(where(o32 lt 60))=6
  cc33(where(o32 ge 120))=4
  ;shift values slightly so that they do not overlap. shift is again grouped by orientation as above
  ss33=o32*0.+1.1
  ss33(where(o32 lt 60))=1.0
  ss33(where(o32 ge 120))=1.2
  loadct2,5
  plot_oo,(1./res32)^2,amp_both,psym=1,xtitle='1/resolution^2',ytitle='amplitudes',yrange=[min([amp_both]),max([amp_both+damp_both])],yst=1,/nodata,title='combined: red=<60!Uo!N; white=60-120!Uo!N; green=>120!Uo!N',xrange=[2d-5,1d-1],xsty=1
  oplot,(1./res32)^2*ss33,amp_both,psym=1,color=6,symsize=this_ss,thick=th3
  err_plot,(1./res32)^2*ss33,(amp_both-damp_both)>0.001,amp_both+damp_both,thick=th3,width=1d-13
  for i=0,31 do plots,(1./res32(i))^2*ss33(i),amp_both(i),psym=1,color=cc33(i),symsize=this_ss,thick=th3
  loadct,5
  endif

  
  ;out={amp_top: amp_top, amp_top_error: damp_top, amp_bot: amp_bot, amp_bot_error: damp_bot, rate_pixel: this_rb, rate_pixel_error: this_drb, tr_flare: tr_flare, er_flare: er_flare}
  out={amp_top: amp_top, $
    amp_top_error: damp_top, $
    amp_bot: amp_bot, $
    amp_bot_error: damp_bot, $
    amp_both: amp_both, $
    amp_both_error: damp_both, $
    rate_pixel: this_rb, $
    rate_pixel_error: this_drb, $
    tr_flare: tr_flare, $
    er_flare: er_flare}
    
  if keyword_set(stop) then stop
  ;format for Gordon
  ;save,p6_cts,p6_cts_e,p6_rate,p6_rate_e,g10,g09,g08,g07,g06,g05,g04,g03,g02,g01,filename=
  ;restore,'/Users/samkrucker/Dropbox/sswidl/stix_flares/s20200607_22_B6/pixel_B6_peak_420s_5-10keV.dat'
  ;IDL> help,p6_cts,p6_cts_e,p6_rate,p6_rate_e
  ;P6_CTS          FLOAT     = Array[32, 12]
  ;P6_CTS_E        FLOAT     = Array[32, 12]
  ;P6_RATE         FLOAT     = Array[32, 12]
  ;P6_RATE_E       FLOAT     = Array[32, 12]
  
  g10=[3,20,22]-1
  g09=[16,14,32]-1
  g08=[21,26,4]-1
  g07=[24,8,28]-1
  g06=[15,27,31]-1
  g05=[6,30,2]-1
  g04=[25,5,23]-1
  g03=[7,29,1]-1
  g02=[12,19,17]-1
  g01=[11,13,18]-1
  p6_cts=out.rate_pixel*(tr_flare(1)-tr_flare(0))
  p6_cts_e=out.rate_pixel_error*(tr_flare(1)-tr_flare(0))
  p6_rate=out.rate_pixel
  p6_rate_e=out.rate_pixel_error
  gh_filename=strmid(anytim(tr_flare(0),/vms),0,11)+'_'+strmid(anytim(tr_flare(0),/vms),12,2)+strmid(anytim(tr_flare(0),/vms),15,2)+'_'+strtrim(fix(er_flare(0)),2)+'-'+strtrim(fix(er_flare(1)),2)+'keV'+'_'+strtrim(fix(tr_flare(1)-tr_flare(0)),2)+'s.dat'
  ;save,p6_cts,p6_cts_e,p6_rate,p6_rate_e,tr_flare,er_flare,g10,g09,g08,g07,g06,g05,g04,g03,g02,g01,filename='pixel_'+gh_filename
  return,out

END

