PRO oplot_int_res,this_g,o32,res32,this_a,this_aa,this_a_range

;interpolate amplitudes between observed angles and plots it (used in stix_display_pixel_trans.pro)



  this_o=min(o32(this_g))+findgen(1001)/1000.*(max(o32(this_g))-min(o32(this_g)))
  this_oo=o32(this_g)
  slist=sort(this_oo)
  this_oo=this_oo(slist)
  this_amp=this_aa(this_g)
  this_amp=this_amp(slist)
  this_int=interpol(this_amp,this_oo,this_o)
  this_c_int=(alog10(this_int)-min(this_a))*(235./this_a_range)+10
  for i=0,999 do plots,sin(this_o(i)/180.*!pi)*res32(this_g(0)),cos(this_o(i)/180.*!pi)*res32(this_g(0)),color=this_c_int(i),psym=-1,symsiz=0.5

  END