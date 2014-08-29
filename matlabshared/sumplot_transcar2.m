function sumploth=sumplot_transcar2(datadir,timingfile,plasmafile,emfile,altlims)

%
%function sumploth=sumplot_transcar(datadir,timingfile,plasmafile,emfile,altlims)
%

sumploth=figure;
set(gcf,'PaperPosition',[0,0,8.5,11]);



%LOAD FLUID AND TIMING DATA
load([datadir,'/',timingfile]);
load([datadir,'/',plasmafile]);



%FONT AND LINEWIDTH SETTINGS
FS=8; LW=4;



%CORRECT TIME VARIABLE FOR PLOTTING PURPOSES
tspan=maxtime-mintime;
tpad=1/6*tspan;                                         %show quiescent ionosphere before event is switched on
tinds=find(time>=mintime+tpad & time<=maxtime);
zinds=find(z>=altlims(1) & z<=altlims(2));
minz=min(z(zinds)); maxz=max(z(zinds)); zspan=maxz-minz;



%FIGURE OUT HOW TO BEST DISPLAY TIME AXIS
% nlab=6;                                                 %nominal number of labels on time axis
% xlims=[time(tinds(1)), time(tinds(length(tinds)))];
% dx=(xlims(2)-xlims(1))/nlab;
% if dx<1/60                                              %labels need seconds to be unambiguous
%     fmt='hms';
% else                                                    %including seconds is pointless
%     fmt='hm';
% end
% xloc=xlims(1):dx:xlims(2);
% xlab=time2str(xloc,'clock',fmt);
zs=zeros(size(time(tinds)));
simdate=[zs', zs', zs', time(tinds)', zs', zs'];
xax=datenum(simdate);



%MEAN ION VELOCITY AND TEMPERATURE
vi_time=(n1_time.*v1_time+n2_time.*v2_time+n3_time.*v3_time+nm_time.*vm_time)./(n1_time+n2_time+n3_time+nm_time);
ti_time=(n1_time.*t1_time+n2_time.*t2_time+n3_time.*t3_time+nm_time.*tm_time)./(n1_time+n2_time+n3_time+nm_time);



%INTERPOLATE VARIABLES FOR IMAGESC PLOTTING PURPOSES
tplot=time(tinds);
zplot=minz:5:maxz;
for k=1:length(tinds);
    tiplot(:,k)=interpolateTRANSCAR(ti_time(zinds,tinds(k)),z(zinds),zplot,'lin','lin');
    teplot(:,k)=interpolateTRANSCAR(te_time(zinds,tinds(k)),z(zinds),zplot,'lin','lin');
    neplot(:,k)=interpolateTRANSCAR(ne_time(zinds,tinds(k)),z(zinds),zplot,'loglin','loglin') ;
    viplot(:,k)=interpolateTRANSCAR(vi_time(zinds,tinds(k)),z(zinds),zplot,'lin','lin');
end



%IONOSPHERIC RESPONSE AS AN ISR WOULD SEE
FS=9;

subplot(412);  pbaspect([8.5 11/6 1]);
h=imagesc(xax,zplot,tiplot);
axis xy;
c=colorbar;
h=gca;
set(h,'FontSize',FS)
datetick;
axis tight;
xlabel('UT');
ylabel('altitude (km)');
%colormap('gray')
ax=axis;
text(ax(1)+.0025,ax(4)-500,'(b)','FontSize',FS+10,'Color',[1 1 1],'FontWeight','bold');
ylabel(c,'T_i (K)')

subplot(413); pbaspect([8.5 11/6 1]);
h=imagesc(xax,zplot,teplot);
axis xy;
c=colorbar;
h=gca;
set(h,'FontSize',FS)
datetick;
axis tight;
xlabel('UT');
ylabel('altitude (km)');
%colormap('gray')
ax=axis;
text(ax(1)+.0025,ax(4)-500,'(c)','FontSize',FS+10,'Color',[1 1 1],'FontWeight','bold');
ylabel(c,'T_e (K)')

subplot(411); pbaspect([8.5 11/6 1]);
maxlogne=max(log10(neplot(:)));
h=imagesc(xax,zplot,log10(neplot),[10 maxlogne]);
axis xy;
c=colorbar;
h=gca;
set(h,'FontSize',FS)
datetick;
axis tight;
xlabel('UT');
ylabel('altitude (km)');
%colormap('gray')
ax=axis;
text(ax(1)+.0025,ax(4)-500,'(a)','FontSize',FS+10,'Color',[1 1 1],'FontWeight','bold');
ylabel(c,'log_{10}(n_e) (m^{-3})')

subplot(414); pbaspect([8.5 11/6 1]);
h=imagesc(xax,zplot,viplot);
axis xy;
c=colorbar;
h=gca;
set(h,'FontSize',FS);
datetick;
axis tight;
xlabel('UT');
ylabel('altitude (km)');
%colormap('gray')
ax=axis;
text(ax(1)+.0025,ax(4)-500,'(d)','FontSize',FS+10,'Color',[1 1 1],'FontWeight','bold');
ylabel(c,'v_i (m/s)')



% %LOAD OPTICAL DATA
% load([datadir,'/',emfile]);
% 
% 
% 
% %WHAT TIME INDEX?
% tind=min(find(timeop>=midevent));
% tinds=find(timeop>=mintime & timeop<=maxtime);
% lt=length(tinds);
% bavg=mean(blambda(:,tinds),2);
% 
% %plot sizing parameters
% FS=8;
% LW=2;
% 
% 
% 
% %BUILD LEGEND (THIS IS A BIT OF A HACK TO GET AROUND MATLAB'S STUPID LEGEND
% %CODE)
% subplot(515); pbaspect([8.5 11/6.5 1]);
% plot(557.7,1e0,'g',557.7,1e0,'r',557.7,1e0,'b',557.7,1e0,'y',557.7,1e0,'m',557.7,1e0,'c','LineWidth',LW);
% legend('O meta-stable','O prompt','N_2^+ 1N','N_2^+ Meinel','N_2 2P','N_2 1P','Location','NorthWest','Orientation','Horizontal');
% legend boxoff;
% 
% 
% 
% %BUILD LINES FOR PLOTTING
% linebase=0.5*min(bavg);
% blines=[bavg,linebase*ones(size(bavg))];
% 
% 
% 
% %IDENTIFY EMISSION GROUPS AND PLOT ON SYNTHETIC SPECTROGRAPH
% %metastable oxygen lines
% hold on;
% for k=1:length(lambdams)
%     ilk=find(lambda==lambdams(k));
%     plot([lambda(ilk),lambda(ilk)],log10(blines(ilk,:)),'g','LineWidth',LW)
% end
% 
% %prompt oxygen emissions
% for k=1:length(lambdapa)
%     ilk=find(lambda==lambdapa(k));
%     plot([lambda(ilk),lambda(ilk)],log10(blines(ilk,:)),'r','LineWidth',LW)
% end
% 
% %1NG
% for k=1:length(lambda1NGvec)
%     ilk=find(lambda==lambda1NGvec(k));
%     plot([lambda(ilk),lambda(ilk)],log10(blines(ilk,:)),'b','LineWidth',LW)
% end
% 
% %Meinel
% for k=1:length(lambdaMeinvec)
%     ilk=find(lambda==lambdaMeinvec(k));
%     plot([lambda(ilk),lambda(ilk)],log10(blines(ilk,:)),'y','LineWidth',LW)
% end
% 
% %2PG
% for k=1:length(lambda2PGvec)
%     ilk=find(lambda==lambda2PGvec(k));
%     plot([lambda(ilk),lambda(ilk)],log10(blines(ilk,:)),'m','LineWidth',LW)
% end
% 
% %1PG
% for k=1:length(lambda1PGvec)
%     ilk=find(lambda==lambda1PGvec(k));
%     plot([lambda(ilk),lambda(ilk)],log10(blines(ilk,:)),'c','LineWidth',LW)
% end
% hold off;
% xlabel('wavelength (nm)','FontSize',FS);
% ylabel('log_{10} Brightness (R)','FontSize',FS);
% %title('Simulated auroral emissions','FontSize',FS)
% set(gca,'FontSize',FS);
% ax=axis;
% axis([200,1000,-1,6]);
% ax=axis;
% text(225,3.5,'(e)','FontSize',20,'Color',[0 0 0],'FontWeight','bold');

end
