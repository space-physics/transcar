function [sumploth,hfigspec ]=sumplot_transcar(datadir,timingfile,plasmafile,emfile,...
                           altlims,onsettime,plotAllSimTime)

doISRinterp = false


sumploth=figure(680); clf(680)
cp = get(sumploth,'pos');
set(sumploth,'pos',[cp(1), cp(2), 900, 740])
%set(gcf,'PaperPosition',[0,0,8.5,11]);

%% LOAD FLUID AND TIMING DATA
%from read_input.m
timingFN=[datadir,'/',timingfile];
display(['sumplot_transcar: Loading ',timingFN])
load(timingFN,'mintime','maxtime','midevent')

plasmaFN = [datadir,'/',plasmafile];
display(['sumplot_transcar: Loading ',plasmaFN])
load(plasmaFN,'-regexp',...
     'n\d_time','v\d_time','t\d_time',...
     '[ntv]m_time','[nt]e_time',...
     '^time$','^z$')
%% FONT AND LINEWIDTH SETTINGS
FS=8; LW=4;
%% CORRECT TIME VARIABLE FOR PLOTTING PURPOSES
%variables saved in read_input.m
%timing(1) dt, (2) starttime, (3) duration, (4) precStart, (5) precEnd
%tspan=maxtime-mintime;
if plotAllSimTime
    tInds = 1:length(time); %take all values
else %restrict to around precipitation time 
    tpad= 0; %MH %1/6*tspan;  %show quiescent ionosphere before precipitation is switched on
    tInds=find(time >= mintime + tpad & ...
               time <= maxtime         );
end

lt = length(tInds);
%% select desireed altitudes
zInds=find(z>=altlims(1) & z<=altlims(2));
%minz=min(z(zinds)); maxz=max(z(zinds)); %zspan=maxz-minz;
minZ = z(zInds(1)); 
maxZ = z(zInds(end));
%% FIGURE OUT HOW TO BEST DISPLAY TIME AXIS
nlab=6;                                                 %nominal number of labels on time axis
xlims=[time(tInds(1)), time(tInds(end))];
dx=(xlims(2)-xlims(1))/nlab;
if dx<1/60                                              %labels need seconds to be unambiguous
    dtfmt='HH:MM:SS';
else                                                    %including seconds is pointless
    dtfmt='HH:MM';
end

simdate = zeros(lt,6); %datevec for datenum
simdate(:,4) = time(tInds);
xax=datenum(simdate);
xOnset= datenum([0 0 0 onsettime 0 0]);
%% MEAN ION VELOCITY AND TEMPERATURE
vi_time=(n1_time.*v1_time+n2_time.*v2_time+n3_time.*v3_time+nm_time.*vm_time)./(n1_time+n2_time+n3_time+nm_time);
ti_time=(n1_time.*t1_time+n2_time.*t2_time+n3_time.*t3_time+nm_time.*tm_time)./(n1_time+n2_time+n3_time+nm_time);
%% INTERPOLATE VARIABLES FOR IMAGESC PLOTTING PURPOSES
% FIXME: maybe this is not necessary? at least, use interp2() for big speedup

%tplot=time(tinds);
%{
nTime
size(ti_time)
size(te_time)
size(ne_time)
size(vi_time)
%}
if doISRinterp
    zPlot=minZ:5:maxZ;
    for k=1:lt
        tiplot(:,k)=interpolateTRANSCAR(ti_time(zInds,tInds(k)),z(zInds),zPlot,'lin','lin');
        teplot(:,k)=interpolateTRANSCAR(te_time(zInds,tInds(k)),z(zInds),zPlot,'lin','lin');
        neplot(:,k)=interpolateTRANSCAR(ne_time(zInds,tInds(k)),z(zInds),zPlot,'loglin','loglin') ;
        viplot(:,k)=interpolateTRANSCAR(vi_time(zInds,tInds(k)),z(zInds),zPlot,'lin','lin');
    end %for
else 
    zPlot = z(zInds);
    tiplot = ti_time(zInds,tInds);
    teplot = te_time(zInds,tInds);
    neplot = ne_time(zInds,tInds);
    viplot = vi_time(zInds,tInds);
end %if doISRinterp
%% IONOSPHERIC RESPONSE AS AN ISR WOULD SEE
FS=9;

h=subtightplot(4,1,2);  %pbaspect([8.5 11/6 1]);
imagesc(xax,zPlot,tiplot)
axis('xy')
c=colorbar;
%set(h,'FontSize',FS)
%datetick(h,'x',dtfmt)
set(h,'xtick',[])
axis('tight')
%xlabel('UT')
ylabel('altitude (km)')
%colormap('gray')
%ax=axis;
%text(ax(1)+.0025,ax(4)-500,'(b)','FontSize',FS+10,'Color',[1 1 1],'FontWeight','bold');
ylabel(c,'T_i (K)')
%plot onset time line
line([xOnset,xOnset],[zPlot(1),zPlot(end)],'color','w','linestyle','--')

%=---------------------

h=subtightplot(4,1,3); %pbaspect([8.5 11/6 1]);
imagesc(xax,zPlot,teplot)
axis('xy')
c=colorbar;
%set(h,'FontSize',FS)
%datetick(h,'x',dtfmt)
set(h,'xtick',[])
axis('tight')
%xlabel('UT')
ylabel('altitude (km)')
%colormap('gray')
%ax=axis;
%text(ax(1)+.0025,ax(4)-500,'(c)',...
           %         'FontSize',FS+10,'Color',[1 1 1],'FontWeight','bold')
ylabel(c,'T_e (K)')
%plot onset time line
line([xOnset,xOnset],[zPlot(1),zPlot(end)],'color','w','linestyle','--')

%=---------------------
h=subtightplot(4,1,1); %pbaspect([8.5 11/6 1]);
maxlogne=max(log10(neplot(:)));
imagesc(xax,zPlot,log10(neplot),[10 maxlogne])
axis('xy')
c=colorbar;
set(h,'FontSize',FS)
%datetick(h,'x',dtfmt)
set(h,'xtick',[])
axis('tight')
%xlabel('UT')
ylabel('altitude (km)')
%colormap('gray')
%ax=axis;
%text(ax(1)+.0025,ax(4)-500,'(a)','FontSize',FS+10,...
   %         'Color',[1 1 1],'FontWeight','bold')
ylabel(c,'log_{10}(n_e) (m^{-3})')
title(h,['Plasma parameters, ',datadir,'  interp=',int2str(doISRinterp)])
%plot onset time line
line([xOnset,xOnset],[zPlot(1),zPlot(end)],'color','w','linestyle','--')

%=---------------------
h=subtightplot(4,1,4); %pbaspect([8.5 11/6 1]);
imagesc(xax,zPlot,viplot)
axis('xy')
c=colorbar;
set(h,'FontSize',FS)
datetick(h,'x',dtfmt)
axis('tight')
xlabel('UT')
ylabel('altitude (km)')
%colormap('gray')
%ax=axis;
%text(ax(1)+.0025,ax(4)-500,'(d)','FontSize',FS+10,...
    %       'Color',[1 1 1],'FontWeight','bold')
ylabel(c,'v_i (m/s)')
%plot onset time line
line([xOnset,xOnset],[zPlot(1),zPlot(end)],'color','w','linestyle','--')

%=---------------------
hfigspec = plot_spectrum(datadir,emfile,midevent,mintime,maxtime,doISRinterp);

end %function
