%GET THE PROCESSED ISR DATA
load ./isrdata_barker.mat;



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
zinds=find(iszb>=0 & iszb<=600);
minz=min(iszb(zinds)); maxz=max(iszb(zinds)); zspan=maxz-minz;

tplot=istimeb;
zplot=100:5:600;
for k=1:length(tplot);
   neplot(:,k)=interpolate(isneb(zinds,k),iszb(zinds),zplot,'loglin','loglin') ;
end



%PLOT AS A RASTERIZED IMAGE W/ SCALABLE FONTS
figure;
set(gcf,'PaperPosition',[0,0,8.5,3]);
FS=12;
h=imagesc(tplot,zplot,flipud(neplot),[10 11.8]);
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
ax=axis;
ylabel(c,'log_{10}(raw n_{e}) (m^{-3})')
print -depsc ne_barker_img.eps;