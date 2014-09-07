%LOAD THE DATA FROM THE RECONSTRUCTED SIMULATION
datelab='11Feb2002_25s_p1';
datafile='isrdata_11Feb2002_25s_p1.mat';
load(['./datafiles/11Feb2002/',datafile]);



%RENAME VARIABLES SO THEY DON'T GET OVERWRITTEN
isz1=isz;
istime1=istime;
isne1=isne;
iste1=iste;
isti1=isti;
isvi1=isvi;

    

%LOAD OTHER PART OF DATA
datelab='11Feb2002_25s_p2';
datafile='isrdata_11Feb2002_25s_p2.mat';
load(['./datafiles/11Feb2002/',datafile]);


%CORRECT TIME VARIABLE FOR PLOTTING PURPOSES
tinds1=1:length(istime1);
tinds=1:length(istime);
zinds=find(isz>=0 & isz<=1000);
minz=min(isz(zinds)); maxz=max(isz(zinds)); zspan=maxz-minz;



%CREATE A DIRECTORY FOR OUTPUT PLOTS
datelab='11Feb2002_combined';
datadir=['./plots/',datelab];
system(['mkdir ',datadir]);



%INTERPOLATE VARIABLES FOR PLOTTING PURPOSES
tplot1=istime1(tinds1);
zplot=minz:5:maxz;
for k=1:length(tinds1);
   tiplot1(:,k)=interpolate(isti1(zinds,tinds1(k)),isz1(zinds),zplot,'lin','lin');
   teplot1(:,k)=interpolate(iste1(zinds,tinds1(k)),isz1(zinds),zplot,'lin','lin');
   neplot1(:,k)=interpolate(isne1(zinds,tinds1(k)),isz1(zinds),zplot,'lin','lin') ;
   viplot1(:,k)=interpolate(isvi1(zinds,tinds1(k)),isz1(zinds),zplot,'lin','lin');
end
tplot=istime(tinds);
for k=1:length(tinds);
   tiplot(:,k)=interpolate(isti(zinds,tinds(k)),isz(zinds),zplot,'lin','lin');
   teplot(:,k)=interpolate(iste(zinds,tinds(k)),isz(zinds),zplot,'lin','lin');
   neplot(:,k)=interpolate(isne(zinds,tinds(k)),isz(zinds),zplot,'lin','lin') ;
   viplot(:,k)=interpolate(isvi(zinds,tinds(k)),isz(zinds),zplot,'lin','lin');
end



%IONOSPHERIC RESPONSE PLOT
FS=9;

tilims=[600 2500];
nelims=[10 12];
telims=[600 3500];
vilims=[-100 800];

figure;
set(gcf,'PaperPosition',[0,0,8.5,11]);

subplot(412);  pbaspect([8.5 11/6 1]); hold on;
h=imagesc(tplot1,zplot,tiplot1,tilims);
imagesc(tplot,zplot,tiplot,tilims);
c=colorbar;
h=gca;
set(h,'FontSize',FS)
xloc=get(h,'Xtick');
xlab=time2str(xloc,'clock','hm');
set(h,'XtickLabel',xlab);
% yloc0=get(h,'YTick');
% dyloc=yloc0-round(minz);
% newyloc=round(maxz)-dyloc;
% set(h,'YTick',sort(newyloc));
% set(h,'YTickLabel',fliplr(yloc0));
xlabel('time (UT hrs)');
ylabel('altitude (km)');
%colormap('gray')
ax=axis;
text(ax(1)+.015,ax(4)-100,'(b)','FontSize',FS+10,'Color',[1 1 1],'FontWeight','bold');
ylabel(c,'T_i (K)')

subplot(413); pbaspect([8.5 11/6 1]); hold on;
h=imagesc(tplot1,zplot,teplot1,telims);
imagesc(tplot,zplot,teplot,telims)
c=colorbar;
h=gca;
set(h,'FontSize',FS)
set(h,'Xtick',xloc);
set(h,'XtickLabel',xlab);
% yloc0=get(h,'YTick');
% dyloc=yloc0-round(minz);
% newyloc=round(maxz)-dyloc;
% set(h,'YTick',sort(newyloc));
% set(h,'YTickLabel',fliplr(yloc0));
xlabel('time (UT hrs)');
ylabel('altitude (km)');
%colormap('gray')
ax=axis;
text(ax(1)+.015,ax(4)-100,'(c)','FontSize',FS+10,'Color',[1 1 1],'FontWeight','bold');
ylabel(c,'T_e (K)')

subplot(411); pbaspect([8.5 11/6 1]); hold on;
h=imagesc(tplot1,zplot,neplot1,nelims);
imagesc(tplot,zplot,neplot,nelims);
c=colorbar;
h=gca;
set(h,'FontSize',FS)
set(h,'Xtick',xloc);
set(h,'XtickLabel',xlab);
% yloc0=get(h,'YTick');
% dyloc=yloc0-round(minz);
% newyloc=round(maxz)-dyloc;
% set(h,'YTick',sort(newyloc));
% set(h,'YTickLabel',fliplr(yloc0));
xlabel('time (UT hrs)');
ylabel('altitude (km)');
%colormap('gray')
ax=axis;
text(ax(1)+.015,ax(4)-100,'(a)','FontSize',FS+10,'Color',[1 1 1],'FontWeight','bold');
ylabel(c,'log_{10}(n_e) (m^{-3})')

subplot(414); pbaspect([8.5 11/6 1]); hold on;
h=imagesc(tplot1,zplot,viplot1,vilims);
imagesc(tplot,zplot,viplot,vilims);
h=gca;
set(h,'FontSize',FS);
set(h,'Xtick',xloc);
set(h,'XtickLabel',xlab);
% yloc0=get(h,'YTick');
% dyloc=yloc0-round(minz);
% newyloc=round(maxz)-dyloc;
% set(h,'YTick',sort(newyloc));
% set(h,'YTickLabel',fliplr(yloc0));
c=colorbar('FontSize',FS);
xlabel('time (UT hrs)');
ylabel('altitude (km)');
%colormap('gray')
ax=axis;
text(ax(1)+.015,ax(4)-100,'(d)','FontSize',FS+10,'Color',[1 1 1],'FontWeight','bold');
ylabel(c,'v_i (m/s)')

print('-depsc',[datadir,'/isrdata_25s.img'])
