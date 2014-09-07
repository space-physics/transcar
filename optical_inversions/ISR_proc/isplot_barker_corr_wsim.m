%GET THE PROCESSED ISR DATA
if ~exist('time')
    load ./isrdata_zerolag_17Feb2001.mat;
    load ./isrdata_17Feb2001.mat;
    load ~/sims/17Feb2001/plasma.mat;
end



%TICK MARKS FOR TIMES AXIS
iminx=min(find(istimeb>=1+1/60));
xloc=istimeb(iminx):2/60:istimeb(length(istimeb));
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
    xlab{k}=[labhr,':',labmin];
end



%GET RID OF SHITTY DATA POINTS
ibad=find(isneb==0);
igood=find(isneb(:)~=0);
mindat=min(isneb(igood));
isneb(ibad)=1e-6*mindat;



%INTERPOLATE FOR IMAGE SCALING
zindsb=find(iszb>=0 & iszb<=600);
zindsf=find(isz>=0 & isz<=600);
zindss=find(z>=0 & z<=600);
minz=min(iszb(zindsb)); maxz=max(iszb(zindsb)); zspan=maxz-minz;

tplotb=istimeb;
zplot=100:5:600; zplot=zplot(:);
for k=1:length(tplotb);
   neplotb(:,k)=interpolate(isneb(zindsb,k),iszb(zindsb),zplot,'lin','lin');
end

mint=min(istimeb); maxt=max(istimeb);
it1f=min(find(istime>=mint)); it2f=min(find(istime>=maxt));
tplotf=istime(it1f:it2f);
for k=1:length(tplotf)
   neplotf(:,k)=interpolate(isne(zindsf,it1f+k-1),isz(zindsf),zplot,'lin','lin');
end

if(max(time)>=24)
  time=time-24;
end
it1s=min(find(time>=mint)); it2s=min(find(time>=maxt));
tplots=time(it1s:it2s);
for k=1:length(tplots)
   neplots(:,k)=log10(interpolate(ne_time(zindss,it1s+k-1),z(zindss),zplot,'loglin','loglin'));
end



%INTERPOLATE TE/TI ONTO BARKER CODED PLOT GRID
ti_time=(n1_time.*t1_time+n2_time.*t2_time+n3_time.*t3_time+nm_time.*tm_time)./(n1_time+n2_time+n3_time+nm_time);
Te_Ti=te_time./ti_time;
for k=1:length(tplotb)
    its=min(find(time>=tplotb(k)));      %ZOH for finding Te/Ti from simulation
    Te_Ti_grid(:,k)=interpolate(Te_Ti(:,its),z,zplot,'lin','lin')';
    fact(:,k)=(Te_Ti_grid(:,k)+1)/2;                     %Te/Ti correction factor vs. altitude
    
    neplotbcorr(:,k)=neplotb(:,k)+log10(fact(:,k));
end



%PLOT AS A RASTERIZED IMAGE W/ SCALABLE FONTS
figure;
set(gcf,'PaperPosition',[0,0,8.5,11]);
FS=9;
clims=[10 11.9];


subplot(511); pbaspect([8.5 11/6 1])
h=imagesc(tplotb,zplot,flipud(neplotb),clims);
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
xlabel('UT (hrs)');
ylabel('altitude (km)');
%colormap('gray')
axb=axis;
ylabel(c,'log_{10}(raw n_{e}) (m^{-3})')
text(axb(1)+.015,axb(3)+100,'(a)','FontSize',FS+10,'Color',[1 1 1],'FontWeight','bold');

subplot(512); pbaspect([8.5 11/6 1])
h=imagesc(tplotf,zplot,flipud(neplotf),clims);
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
xlabel('UT (hrs)');
ylabel('altitude (km)');
%colormap('gray')
axis(axb);
ylabel(c,'log_{10}(fitted n_{e}) (m^{-3})')
text(axb(1)+.015,axb(3)+100,'(b)','FontSize',FS+10,'Color',[1 1 1],'FontWeight','bold');

subplot(513); pbaspect([8.5 11/6 1])
h=imagesc(tplots-1/240,zplot,flipud(neplots),clims);
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
xlabel('UT (hrs)');
ylabel('altitude (km)');
%colormap('gray')
axis(axb);
ylabel(c,'log_{10}(simulated n_{e}) (m^{-3})')
text(axb(1)+.015,axb(3)+100,'(c)','FontSize',FS+10,'Color',[1 1 1],'FontWeight','bold');

subplot(515); pbaspect([8.5 11/6 1])
h=imagesc(tplotb,zplot,flipud(neplotbcorr),clims);
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
xlabel('UT (hrs)');
ylabel('altitude (km)');
%colormap('gray')
axis(axb);
ylabel(c,'log_{10}(T_e/T_i adj. n_{e}) (m^{-3})')
text(axb(1)+.015,axb(3)+100,'(e)','FontSize',FS+10,'Color',[1 1 1],'FontWeight','bold');

subplot(514); pbaspect([8.5 11/6 1])
h=imagesc(tplotb,zplot,flipud(Te_Ti_grid));
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
xlabel('UT (hrs)');
ylabel('altitude (km)');
%colormap('gray')
axis(axb);
ylabel(c,'T_e/T_i')
text(axb(1)+.015,axb(3)+100,'(d)','FontSize',FS+10,'Color',[1 1 1],'FontWeight','bold');

print -depsc ne_corr_wsim_img.eps;