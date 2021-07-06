
pro stx_diplay_counts, data, wwindow=wwindow


loadct,5
range_title = ''

this_rb = data.rate_pixel
this_dr = data.rate_pixel_error

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
res10=[7.1,10.2,14.6,20.9,29.8,42.7,61.0,87.3,124.9,178.6]
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
 
window,wwindow,xsize=900,ysize=800
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
     
plot,[0.5,1.5,2.5,3.5],this_rb(g_plot(i),0:3),xtitle=' ',ytitle=' ',psym=-1,charsi=chs,yrange=[0,max(this_rb(g_plot,0:7))],noe=i,xtickname=replicate(' ',9),ytickname=replicate(' ',9),xticks=8,xminor=1,xticklen=1d-22,title=this_title
oplot,[0.5,1.5,2.5,3.5],this_rb(g_plot(i),0:3),psym=-1,color=c_top
errplot,[0.5,1.5,2.5,3.5],(this_rb(g_plot(i),0:3)-this_dr(g_plot(i),0:3)),(this_rb(g_plot(i),0:3)+this_dr(g_plot(i),0:3)),thick=th3,color=c_top
oplot,[0.5,1.5,2.5,3.5],this_rb(g_plot(i),4:7),psym=-1,color=c_bot
errplot,[0.5,1.5,2.5,3.5],(this_rb(g_plot(i),4:7)-this_dr(g_plot(i),4:7)),(this_rb(g_plot(i),4:7)+this_dr(g_plot(i),4:7)),thick=th3,color=c_bot
   
endif else begin

plot,[0.5,1.5,2.5,3.5],this_rb(g_plot(i),0:3),xtitle=' ',ytitle='cts/s/keV',psym=-1,charsi=chs,yrange=[0,max(this_rb(g_plot,*))],noe=i,xtickname=[' ','A',' ','B',' ','C',' ','D',' '],xticks=8,xminor=1,xticklen=1d-22,title=this_title
     
oplot,[0.5,1.5,2.5,3.5],this_rb(g_plot(i),0:3),psym=-1,color=c_top
errplot,[0.5,1.5,2.5,3.5],(this_rb(g_plot(i),0:3)-this_dr(g_plot(i),0:3)),(this_rb(g_plot(i),0:3)+this_dr(g_plot(i),0:3)),thick=th3,color=c_top
oplot,[0.5,1.5,2.5,3.5],this_rb(g_plot(i),4:7),psym=-1,color=c_bot
errplot,[0.5,1.5,2.5,3.5],(this_rb(g_plot(i),4:7)-this_dr(g_plot(i),4:7)),(this_rb(g_plot(i),4:7)+this_dr(g_plot(i),4:7)),thick=th3,color=c_bot
   
endelse
endfor
 
xyouts,0.5,1-ymargin_top/2.5,range_title,/normal,chars=1.6,ali=0.5
 
!p.position = [0, 0, 0, 0]
 
 
 
 end