%LOAD DATA KERNELS AND DATA
%datelab='11Feb2002';
datelab='20Nov2001_short';
load(['./kernel_build/kernel_resampled_',datelab,'.mat']);
load(['./NIR_proc/datafiles/',datelab,'/filelist.mat']);
system(['mkdir inversion_plots/',datelab]);



%RESAMPLED DATA KERNEL
M=bchrresamp;



LW=2;
FS=16;
for k=1:length(filelist)
    load(['./NIR_proc/datafiles/',datelab,'/',filelist{k},'_rayleighs.mat']);

    %INVERSION PARAMS
%    w=100*ones(size(b)); w(2:4)=1; w=w/sum(w);          %give little weight to 1PG, does well for 17Feb2001
%    w=ones(size(b)); w=w/sum(w);                        %works well for 19Nov2001... sort of
    w=100*ones(size(b)); w(2:4)=1; w=w/sum(w);                        %works well for 20Nov2001.    
    beta=2.5;                                           %works quite well for BOT and 20Nov2001
%    beta=2;                                             %works well for 11Feb2002
    ic=ones(nphi,1)./phiNbin(:)*1e20;                   %works for FT and BOT


    %FORGET NOISE CALCULATIONS FOR NOW, CHECK WITH EMBRY-RIDDLE GROUP LATER
    sigman=zeros(size(b));


    %ME INVERSION
    [phimem(:,k),chir]=meberg(b,M,ic,beta,w,sigman);
    phibasis(:,k)=phimem(:,k);
    phimem(:,k)=phimem(:,k).*phiNbin(:);
    
    
    %RELATIVE ERROR COMPUTATION
    relerr_binned(:,k)=abs(M*phibasis(:,k)-b)./b;
    errbinned(k)=sum(relerr_binned(:,k).*w);


    %ENERGY BIN CENTERS
    erplot=er(1:length(er)-1)+0.5*diff(er);


    %PLOT THE RECONSTRUCTION
    loglog(erplot,phimem(:,k),'o-','LineWidth',LW,'MarkerSize',10,'MarkerFaceColor',[1 1 1]);
    axis([70,2e4,3e2,5e7]);
    set(gca,'FontSize',FS);
    ylabel('\phi(E) (cm^{-2} s^{-1} eV^{-1} sr^{-1})');
    xlabel('energy (eV)');
    title(sprintf('Inversion of data from:  %s UT',lablist{k}));
    ax=axis;
%    text(ax(1)+20,ax(4)-5e6,pan{k},'FontSize',FS+8);
    print('-depsc',['inversion_plots/',datelab,'/invdata_',lablist{k},'.eps'])
    clf;
end



%ENERGY-TIME SPECTROGRAM SUMMARY
[h]=E_t_spectrogram(tset,er,phimem);
print(h,'-depsc',['inversion_plots/',datelab,'/spectrogram.eps']);
figure(h);
xloc=get(gca,'XTick');
xlab=get(gca,'XTickLabel');
ax=axis;



%ERROR ANALYSIS PLOT
zs=zeros(size(tset(:)));
dnum=datenum([zs,zs,zs,tset(:),zs,zs]);
figure;
set(gcf,'PaperPosition',[0,0,8.5,3.5]);
plot(dnum,errbinned*100,'LineWidth',2);
set(gca,'FontSize',12);
set(gca,'XTick',xloc);
set(gca,'XTickLabel',xlab);
axis([ax(1:2),0,100]);
xlabel('UT');
ylabel('% relative error');
title('Average relative error for each reconstruction');
print('-depsc',['inversion_plots/',datelab,'/recon_errors.eps']);



%CREATE AN INPUT FILE FOR TRANSCAR
create_precinput(tset,er,phimem,[datelab,'.dat']);
close all;



%OLD CODE

% %TICK MARKS FOR TIMES AXIS
% %xloc=5/60:5/60:4/3;            %20Nov2001
% xloc=23.5:5/60:23.6667;     %11Feb2002
% for k=1:length(xloc)
%     toi=xloc(k);
%     labhr=int2str(floor(toi));
%     labmin=int2str(max(0,round(toi*60-str2num(labhr)*60)));
%     if length(labmin)<2
%         labmin=['0',labmin];
%     end
%     labsec=int2str(max(0,round(toi*3600-str2num(labmin)*60-str2num(labhr)*3600)));
%     if length(labsec)<2
%         labsec=['0',labsec];
%     end
%     %xlab{k}=[labhr,':',labmin,':',labsec];
%     xlab{k}=[labhr,':',labmin];
% end



% fid=fopen([datelab,'.dat'],'w');
% for k=1:length(filelist)
%     data=[er',[phimem(:,k);-1.0]]';
%     fprintf(fid,'%f\n',round(tset(k)*3600));
%     fprintf(fid,'%e %e\n',data);
% end
% lt=length(tset)
% fprintf(fid,'%f\n',round((tset(lt)+(tset(lt)-tset(lt-1)))*3600));
% fprintf(fid,'-1.0 -1.0\n');
% fclose(fid);



% espec=logspace(log10(50),log10(2e4),20);
% lt=length(tset);
% for k=1:lt
%   phiEmem(:,k)=erplot(:).*phimem(:,k);
%   phiEmemplot(:,k)=interpolate(phiEmem(:,k),erplot,espec,'loglin','loglin');
% end
% 
% figure;
% set(gcf,'PaperPosition',[0,0,8.5,3.5]);
% imagesc(tset,log10(espec),flipud(log10(phiEmemplot)));
% h=gca;
% %set(h,'LineStyle','none');
% c=colorbar;
% set(h,'FontSize',12);
% set(h,'XTick',xloc);
% set(h,'XTickLabel',xlab);
% maxE=max(log10(espec));
% minE=min(log10(espec));
% yloc0=get(h,'YTick');
% dyloc=yloc0-round(minE);
% newyloc=round(maxE)-dyloc;
% set(h,'YTick',sort(newyloc));
% set(h,'YTickLabel',fliplr(yloc0));
% title('Energy-time spectrogram of reconstructed e^- precip.');
% ylabel(c,'log_{10} E{\cdot}\phi (cm^{-2} s^{-1} sr^{-1})')
% xlabel('UT');
% ylabel('energy (eV)');
% ax=axis;
% print('-depsc',['inversion_plots/',datelab,'/spectrogram.eps']);