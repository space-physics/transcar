load esr_cp4
for i=1:length(per)
   lat(i)=per(i).lat;
end
Zpol=interp1(lat_acf/100,r_h,lat);
Lpol=interp1(lat_acf/100,lon_acf/100,lat);
nb_pos=length(per);
n_temps=fix(length(time)/nb_pos);
for i=1:n_temps;
   ind_t0=(i-1)*nb_pos;
   t_pos(i)=time(ind_t0+1);
   for j=1:nb_pos;
      ind_t=ind_t0+j;
      Ne_time(j,i)=exp(interp1(z,log(ne_time(:,ind_t)),Zpol(j)));
      Te_time(j,i)=    interp1(z,    te_time(:,ind_t) ,Zpol(j));
      T1_time(j,i)=    interp1(z,    t1t_time(:,ind_t) ,Zpol(j));
      V1_time(j,i)=    vperpnord_time(ind_t);
   end
end

w_cp4=figure(11);
time_limits=[9.5 10.5];
time_limits=[min(time) max(time)];
color_w_cp4=get(w_cp4,'color');
clf;

pos_fig=[.10 .15 .74 .75];

pos_axe=[pos_fig(1) pos_fig(2) pos_fig(3) pos_fig(4)/4-.01];
pos_labely=[1.05 .5 0];


hNe=axes('Position',pos_axe+3*[0 pos_axe(4)+.01 0 0]);
pcolor(t_pos,lat,log10(Ne_time));
shading interp
set(hNe,'xcolor','blue','ycolor','blue','xlim',time_limits,'clim',[10.5 11.5])
set(hNe,'XTick',[])
tit='N_e [ m^{-3} ]';
ht=pal(hNe(1),tit);
hNe=[hNe(1),ht];
label_y(hNe(1),'center');

hTe=axes('Position',pos_axe+2*[0 pos_axe(4)+.01 0 0]);
pcolor(t_pos,lat,Te_time);
shading interp
set(hTe,'xcolor','blue','ycolor','blue','xlim',time_limits,'clim',[1000 2000])
set(hTe,'XTick',[])
tit='T_e [ K ]';
ht=pal(hTe(1),tit);
hTe=[hTe(1),ht];
label_y(hTe(1),'bottom');

hT1=axes('Position',pos_axe+[0 pos_axe(4)+.01 0 0]);
pcolor(t_pos,lat,T1_time);
shading interp
set(hT1,'xcolor','blue','ycolor','blue','xlim',time_limits,'clim',[500 1500])
set(hT1,'XTick',[])
tit='T_i [ K ]';
ht=pal(hT1(1),tit);
hT1=[hT1(1),ht];
label_y(hT1(1),'bottom');

hV1=axes('Position',pos_axe);
pcolor(t_pos,lat,V1_time);
shading interp
set(hV1,'xcolor','blue','ycolor','blue','xlim',time_limits)
set(hV1,'XTick',[])
tit='V_i [ m.s^{-1} ]';
ht=pal(hV1(1),tit);
hV1=[hV1(1),ht];
label_y(hV1(1),'bottom');
fig=hV1(1);
lab_time;
  
h=axes('Position',pos_fig,'unit','norm','Visible','off');
lab_alti=get(h,'ylabel');
set(lab_alti,'position',[-.07 .5 0],'unit','norm','string','Geographic Latitude','color','blue','vertical','bottom','Visible','on');
