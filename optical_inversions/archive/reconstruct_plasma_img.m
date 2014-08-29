%LOAD THE DATA FROM THE RECONSTRUCTED SIMULATION
datelab='11Feb2002';
if ~exist('ti_time')
    load(['~/sims/',datelab,'/plasma.mat']);
    load(['~/sims/',datelab,'/dat.mat']);
    load(['./NIR_proc/datafiles/',datelab,'/filelist.mat']);
end


%CORRECT TIME VARIABLE FOR PLOTTING PURPOSES
tinds=find(time>=mintime+5/60 & time<=maxtime);
zinds=find(z>=0 & z<=1000);
minz=min(z(zinds)); maxz=max(z(zinds)); zspan=maxz-minz;


%TICK MARKS FOR TIMES AXIS
%xloc=.5:5/60:(maxtime-24);         %20 Nov 2001
xloc=23.5:5/60:23.6667;     %11Feb2002
for k=1:length(xloc)
    toi=xloc(k);
    labhr=int2str(floor(toi));
    labmin=int2str(max(0,round(toi*60-str2num(labhr)*60)));
    if length(labmin)<2
        labmin=['0',labmin];
    end
    labsec=int2str(max(0,round(toi*3600-str2num(labmin)*60-str2num(labhr)*3600)));
    if length(labsec)<2
        labsec=['0',labsec];
    end
    %xlab{k}=[labhr,':',labmin,':',labsec];
    xlab{k}=[labhr,':',labmin];
end



%MEAN ION VELOCITY AND TEMPERATURE
vi_time=(n1_time.*v1_time+n2_time.*v2_time+n3_time.*v3_time+nm_time.*vm_time)./(n1_time+n2_time+n3_time+nm_time);
ti_time=(n1_time.*t1_time+n2_time.*t2_time+n3_time.*t3_time+nm_time.*tm_time)./(n1_time+n2_time+n3_time+nm_time);



%INTERPOLATE VARIABLES FOR PLOTTING PURPOSES
if max(time(tinds))>24
  tplot=time(tinds)-24;
else
  tplot=time(tinds); 
end
zplot=90:5:995;
for k=1:length(tinds);
   tiplot(:,k)=interpolate(ti_time(zinds,tinds(k)),z(zinds),zplot,'lin','lin');
   teplot(:,k)=interpolate(te_time(zinds,tinds(k)),z(zinds),zplot,'lin','lin');
   neplot(:,k)=interpolate(ne_time(zinds,tinds(k)),z(zinds),zplot,'loglin','loglin') ;
   viplot(:,k)=interpolate(vi_time(zinds,tinds(k)),z(zinds),zplot,'lin','lin');
end



%IONOSPHERIC RESPONSE PLOT
FS=9;
datadir=['./inversion_plots/',datelab];

%nelims=[10 11.9];           %default
tilims=[600 2500];
nelims=[10 12];
telims=[600 3500];
vilims=[-100 800];

figure;
set(gcf,'PaperPosition',[0,0,8.5,11]);

subplot(412);  pbaspect([8.5 11/6 1]);
h=imagesc(tplot,zplot,flipud(tiplot),tilims);
c=colorbar;
h=gca;
set(h,'FontSize',FS)
set(h,'Xtick',xloc);
set(h,'XtickLabel',xlab);
yloc0=get(h,'YTick');
dyloc=yloc0-round(minz);
newyloc=round(maxz)-dyloc;
set(h,'YTick',sort(newyloc));
set(h,'YTickLabel',fliplr(yloc0));
xlabel('time (UT hrs)');
ylabel('altitude (km)');
%colormap('gray')
ax=axis;
text(ax(1)+.015,ax(3)+100,'(b)','FontSize',FS+10,'Color',[1 1 1],'FontWeight','bold');
ylabel(c,'T_i (K)')

subplot(413); pbaspect([8.5 11/6 1]);
h=imagesc(tplot,zplot,flipud(teplot),telims);
c=colorbar;
h=gca;
set(h,'FontSize',FS)
set(h,'Xtick',xloc);
set(h,'XtickLabel',xlab);
yloc0=get(h,'YTick');
dyloc=yloc0-round(minz);
newyloc=round(maxz)-dyloc;
set(h,'YTick',sort(newyloc));
set(h,'YTickLabel',fliplr(yloc0));
xlabel('time (UT hrs)');
ylabel('altitude (km)');
%colormap('gray')
ax=axis;
text(ax(1)+.015,ax(3)+100,'(c)','FontSize',FS+10,'Color',[1 1 1],'FontWeight','bold');
ylabel(c,'T_e (K)')

subplot(411); pbaspect([8.5 11/6 1]);
h=imagesc(tplot,zplot,flipud(log10(neplot)),nelims);
c=colorbar;
h=gca;
set(h,'FontSize',FS)
set(h,'Xtick',xloc);
set(h,'XtickLabel',xlab);
yloc0=get(h,'YTick');
dyloc=yloc0-round(minz);
newyloc=round(maxz)-dyloc;
set(h,'YTick',sort(newyloc));
set(h,'YTickLabel',fliplr(yloc0));
xlabel('time (UT hrs)');
ylabel('altitude (km)');
%colormap('gray')
ax=axis;
text(ax(1)+.015,ax(3)+100,'(a)','FontSize',FS+10,'Color',[1 1 1],'FontWeight','bold');
ylabel(c,'log_{10}(n_e) (m^{-3})')

subplot(414); pbaspect([8.5 11/6 1]);
h=imagesc(tplot,zplot,flipud(viplot),vilims);
h=gca;
set(h,'FontSize',FS);
set(h,'Xtick',xloc);
set(h,'XtickLabel',xlab);
yloc0=get(h,'YTick');
dyloc=yloc0-round(minz);
newyloc=round(maxz)-dyloc;
set(h,'YTick',sort(newyloc));
set(h,'YTickLabel',fliplr(yloc0));
c=colorbar('FontSize',FS);
xlabel('time (UT hrs)');
ylabel('altitude (km)');
%colormap('gray')
ax=axis;
text(ax(1)+.015,ax(3)+100,'(d)','FontSize',FS+10,'Color',[1 1 1],'FontWeight','bold');
ylabel(c,'v_i (m/s)')

print('-depsc',[datadir,'reconstructed_plasma_resp_img.eps']);



% %TE/TI FOR BARKER CODE CORRECTIONS
% figure;
% h=imagesc(tplot,zplot,flipud(teplot./tiplot));
% h=gca;
% set(h,'FontSize',16);
% set(h,'Xtick',xloc);
% set(h,'XtickLabel',xlab);
% yloc0=get(h,'YTick');
% dyloc=yloc0-round(minz);
% newyloc=round(maxz)-dyloc;
% set(h,'YTick',sort(newyloc));
% set(h,'YTickLabel',fliplr(yloc0));
% c=colorbar('FontSize',16);
% xlabel('time (UT hrs)');
% ylabel('altitude (km)');
% %colormap('gray')
% ax=axis;
% ylabel(c,'T_e/T_i')
% print('-depsc',[datadir,'Te_Ti.eps']);
% 
% tbark=time(tinds)-24;
% zbark=z(zinds);
% te_ti=te_time(zinds,tinds)./ti_time(zinds,tinds);
% nesim=ne_time(zinds,tinds);
% 
% save barker.mat tbark zbark te_ti nesim;
