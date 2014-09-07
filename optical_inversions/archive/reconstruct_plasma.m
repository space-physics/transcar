function reconstruct_plasma(dataDir,writePlots)
%% 
% after running process_transcar() to convert TRANSCAR output to Matlab format,
% reconstruct_plasma() outputs the Ne, Ti, Te, vi plasma parameters
%
addpath('../matlabshared')
%% user parameters
if nargin<2, writePlots = false; end

MinZkm = 0; 
MaxZkm = 1000;
%% LOAD THE DATA FROM THE RECONSTRUCTED SIMULATION
datelab='20Nov2001'; %'11Feb2002';

plasmaFN = [dataDir,'/plasma.mat'];
datFN = [dataDir,'/dat.mat'];
flFN = ['NIR_proc/datafiles/',datelab,'/filelist.mat'];
plotdir=['inversion_plots/',datelab];
%% do it
display(['loading:  ',plasmaFN])
load(plasmaFN,'-regexp',...
     'n\d_time','v\d_time','t\d_time',...
     '[ntv]m_time','[nt]e_time',...
     '^time$','^z$')

display(['loading TRANSCAR input data file:  ',datFN])
load(datFN,'mintime','maxtime')

display(['Loading filelist:  ',flFN])
load(flFN)
%% CORRECT TIME VARIABLE FOR PLOTTING PURPOSES
tinds=find(time>=mintime+5/60 & time<=maxtime);
zinds=find(z>=MinZkm & z<=MaxZkm);
minz=min(z(zinds)); maxz=max(z(zinds)); 
%zspan=maxz-minz;
%% TICK MARKS FOR TIMES AXIS
switch datelab
    case '20Nov2001',xloc=.5:5/60:(maxtime-24);
    case '11Feb2002', xloc=23.5:5/60:23.6667; 
    otherwise, warning(['unknown date ',datelab])
end
%{
for k=1:length(xloc)
    toi=xloc(k);
    labhr=int2str(floor(toi));
    labmin=num2str(max(0,round(toi*60-str2num(labhr)*60)),'%02d'); %#ok<ST2NM>

    labsec=num2str(max(0,round(toi*3600-str2num(labmin)*60-str2num(labhr)*3600)),'%02d'); %#ok<ST2NM>

    %xlab{k}=[labhr,':',labmin,':',labsec];
    xlab{k}=[labhr,':',labmin];
end
%}


%% MEAN ION VELOCITY AND TEMPERATURE
% the n1 etc. variable come from read_tra.m -- a scary script indeed
vi_time=(n1_time.*v1_time + n2_time.*v2_time + n3_time.*v3_time + nm_time.*vm_time) ./ ...
        (n1_time+n2_time+n3_time+nm_time);
ti_time=(n1_time.*t1_time + n2_time.*t2_time + n3_time.*t3_time + nm_time.*tm_time) ./ ...
        (n1_time+n2_time+n3_time+nm_time);
%te_time is loaded from plasma.dat
%ne_time is loaded from plasma.dat
%% INTERPOLATE VARIABLES FOR PLOTTING PURPOSES
if max(time(tinds))>24
  tplot=time(tinds)-24;
else
  tplot=time(tinds); 
end
zplot=90:5:995;
for k=1:length(tinds);
   tiplot(:,k)=interpolateTRANSCAR(ti_time(zinds,tinds(k)),z(zinds),zplot,'lin','lin');
   teplot(:,k)=interpolateTRANSCAR(te_time(zinds,tinds(k)),z(zinds),zplot,'lin','lin');
   neplot(:,k)=interpolateTRANSCAR(ne_time(zinds,tinds(k)),z(zinds),zplot,'loglin','loglin') ;
   viplot(:,k)=interpolateTRANSCAR(vi_time(zinds,tinds(k)),z(zinds),zplot,'lin','lin');
end
%% IONOSPHERIC RESPONSE PLOT
FS=9;

%nelims=[10 11.9];           %default
tilims=[600 2500];
nelims=[10 12];
telims=[600 3500];
vilims=[-100 800];

figure(1)
%set(gcf,'PaperPosition',[0,0,8.5,11]);

h=subtightplot(4,1,2);
pbaspect([8.5 11/6 1]);
imagesc(tplot,zplot,flipud(tiplot),tilims)
c=colorbar;
set(h,'ydir','normal','xtick',[])
%datetick(h,'x','keeplimits','keepticks')
%xlabel('time (UT hrs)');
ylabel('altitude (km)');
%colormap('gray')
ax=axis;
text(ax(1)+.015,ax(3)+100,'(b)','FontSize',FS+10,'Color',[1 1 1],'FontWeight','bold');
ylabel(c,'T_i (K)')

h=subtightplot(4,1,3);
pbaspect([8.5 11/6 1]);
imagesc(tplot,zplot,flipud(teplot),telims)
c=colorbar;
set(h,'ydir','normal','xtick',[])
%datetick(h,'x','keeplimits','keepticks')
%set(h,'FontSize',FS)
%set(h,'Xtick',xloc);
%set(h,'XtickLabel',xlab);
%yloc0=get(h,'YTick');
%dyloc=yloc0-round(minz);
%newyloc=round(maxz)-dyloc;
%set(h,'YTick',sort(newyloc));
%set(h,'YTickLabel',fliplr(yloc0));
%xlabel('time (UT hrs)');
ylabel('altitude (km)');
%colormap('gray')
ax=axis;
text(ax(1)+.015,ax(3)+100,'(c)','FontSize',FS+10,'Color',[1 1 1],'FontWeight','bold');
ylabel(c,'T_e (K)')

h=subtightplot(4,1,1); 
pbaspect([8.5 11/6 1]);
imagesc(tplot,zplot,flipud(log10(neplot)),nelims)
c=colorbar;
set(h,'ydir','normal','xtick',[])
%datetick(h,'x','keeplimits','keepticks')
%set(h,'FontSize',FS)
%set(h,'Xtick',xloc);
%set(h,'XtickLabel',xlab);
%yloc0=get(h,'YTick');
%dyloc=yloc0-round(minz);
%newyloc=round(maxz)-dyloc;
%set(h,'YTick',sort(newyloc));
%set(h,'YTickLabel',fliplr(yloc0));
%xlabel('time (UT hrs)');
ylabel('altitude (km)');
%colormap('gray')
ax=axis;
text(ax(1)+.015,ax(3)+100,'(a)','FontSize',FS+10,'Color',[1 1 1],'FontWeight','bold');
ylabel(c,'log_{10}(n_e) (m^{-3})')

h=subtightplot(4,1,4);
pbaspect([8.5 11/6 1]);
imagesc(tplot,zplot,flipud(viplot),vilims)
set(h,'ydir','normal')
datetick(h,'x','keeplimits','keepticks')
%set(h,'FontSize',FS);
%set(h,'Xtick',xloc);
%set(h,'XtickLabel',xlab);
%yloc0=get(h,'YTick');
%dyloc=yloc0-round(minz);
%newyloc=round(maxz)-dyloc;
%set(h,'YTick',sort(newyloc));
%set(h,'YTickLabel',fliplr(yloc0));
c=colorbar('FontSize',FS);
xlabel('time (UT hrs)');
ylabel('altitude (km)');
%colormap('gray')
ax=axis;
text(ax(1)+.015,ax(3)+100,'(d)','FontSize',FS+10,'Color',[1 1 1],'FontWeight','bold');
ylabel(c,'v_i (m/s)')

if writePlots
    epsFN = [plotdir,'/reconstructed_plasma_resp_img.eps'];
    display(['Writing to ',epsFN])
    print('-depsc',epsFN);
end

%% cleanup
rmpath('../matlabshared')
end %function

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
