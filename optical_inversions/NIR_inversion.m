function [phimem, errbinned,utprec,energy,diffnumflux] = NIR_inversion(writePrecip,writePlots)
%% user parameters
datelab='20Nov2001'; %'11Feb2002';

if nargin<1, writePrecip=false; end
if nargin<2, writePlots=false; end

addpath('../matlabshared')

outputPlotDir = ['inversion_plots/',datelab];
%% LOAD DATA KERNELS AND DATA

% get phiNbin, er, phiEbin
kernelmat = ['kernel_build/archive/kernel_resampled_',datelab,'.mat'];
display(['loading ',kernelmat])
load(kernelmat,'er','phiNbin','phiEbin','bchrresamp','nphi')
phiNbin = phiNbin(:); %#ok<NODEF> %needs to be column vector

% get filelist 
flistmat = ['NIR_proc/datafiles/',datelab,'/filelist.mat'];
display(['loading ',flistmat])
load(flistmat,'filelist','lablist','tset')
nFiles = length(filelist); %#ok<USENS>

% make output directory if it doesn't already exist
if ~exist(outputPlotDir,'dir')
    display(['creating ',outputPlotDir])
    mkdir(outputPlotDir)
end %if
%% RESAMPLED DATA KERNEL
M=bchrresamp;
LW=2;
%% preallocation
phimem = nan(nphi,nFiles); errbinned = nan(nphi,1);

k=1;
 currNIRfn = ['NIR_proc/datafiles/',datelab,'/',filelist{k},'_rayleighs.mat'];
[phimem(:,k),errbinned(k),erplot] = invert_NIR(M,phiNbin,nphi,er,currNIRfn,datelab);
figure(1)
hLog = loglog(erplot,phimem(:,k),'o-','LineWidth',LW,'MarkerSize',10,'MarkerFaceColor',[1 1 1]);
axis([70,2e4,3e2,5e7])
ylabel('$\phi(E)$ (cm$^{-2}$ s$^{-1}$ eV$^{-1}$ sr$^{-1}$)','interpreter','latex')
xlabel('energy (eV)')



for k=1:nFiles
    currNIRfn = ['NIR_proc/datafiles/',datelab,'/',filelist{k},'_rayleighs.mat'];

    % INVERT
    [phimem(:,k),errbinned(k),erplot] = invert_NIR(M,phiNbin,nphi,er,currNIRfn,datelab);

    %PLOT THE RECONSTRUCTION
    %loglog(erplot,phimem(:,k),'o-','LineWidth',LW,'MarkerSize',10,'MarkerFaceColor',[1 1 1]);
    set(hLog,'xdata',erplot,'ydata',phimem(:,k))
    title(sprintf('Inversion of data from:  %s UT',lablist{k})); %#ok<USENS>
    %ax=axis;   text(ax(1)+20,ax(4)-5e6,pan{k},'FontSize',FS+8);
    drawnow
    if writePlots
        currEPSfn = ['inversion_plots/',datelab,'/invdata_',lablist{k},'.eps'];
        display(['Writing ',currEPSfn])
        print('-depsc',currEPSfn) 
    end
end %for
%% ENERGY-TIME SPECTROGRAM SUMMARY
[h,dnum]=E_t_spectrogram(tset,er,phimem);
print(h,'-depsc',['inversion_plots/',datelab,'/spectrogram.eps']);
figure(h)
%xloc=get(gca,'XTick');
%xlab=get(gca,'XTickLabel');
%ax=axis;
%% ERROR ANALYSIS PLOT
figure(2)
set(gcf,'PaperPosition',[0,0,8.5,3.5]);
plot(dnum,errbinned*100,'LineWidth',2);
%set(gca,'FontSize',12,'XTick',xloc,'XTickLabel',xlab);
%axis([ax(1:2),0,100]);
datetick(gca,'x','HH:MM','keeplimits')
xlabel('UT'); ylabel('% relative error');
title('Average relative error for each reconstruction');
print('-depsc',['inversion_plots/',datelab,'/recon_errors.eps']);
%% CREATE AN INPUT FILE FOR TRANSCAR
if writePrecip
    PrecipDataFN = [datelab,'.dat'];
    display(['Writing ASCII Precipitation TRANSCAR input file to: ',PrecipDataFN])
    create_precinput(tset,er,phimem,PrecipDataFN)
end
%% CREATE A .MAT FILE FOR PRECIPITATION
utprec=tset; 
le=length(er);
energy=er(1:le-1)+1/2*diff(er);
energy=energy(:); 
diffnumflux=phimem; 

if writePrecip
    PrecipFN = [datelab,'_precip.mat'];
    display(['Saving ',PrecipFN])
    save(PrecipFN,'utprec','energy','diffnumflux')
end %if
%% cleanup
rmpath('../matlabshared')
if nargout==0, clear, end
end %function

function [phimem,errbinned,erplot] = invert_NIR(M,phiNbin,nphi,er,currNIRfn,datelab)
    display(['Loading ',currNIRfn])
    load(currNIRfn)    
%% INVERSION PARAMS
switch datelab
    case '20Nov2001', w=100*ones(size(b)); w(2:4)=1; w=w/sum(w); %works well for 20Nov2001.    
    case '17Feb2001', w=100*ones(size(b)); w(2:4)=1; w=w/sum(w); %give little weight to 1PG, does well for 17Feb2001
    case '19Nov2001', w=    ones(size(b));           w=w/sum(w); %works well for 19Nov2001... sort of
    otherwise, warning(['Unknown date ',datelab])
end
   
    betaPhi=2.5;                                           %works quite well for BOT and 20Nov2001
%    betaPhi=2;                                             %works well for 11Feb2002
    ic=ones(nphi,1) ./ phiNbin*1e20;                    %works for FT and BOT
%% FORGET NOISE CALCULATIONS FOR NOW, CHECK WITH EMBRY-RIDDLE GROUP LATER
    sigman=zeros(size(b));
%% Max Entropy INVERSION
    phimem=meberg(b,M,ic,betaPhi,w,sigman);
    phibasis=phimem;
    phimem=phimem .* phiNbin;
%% RELATIVE ERROR COMPUTATION
    relerr_binned=abs( M*phibasis - b ) ./ b;
    errbinned=sum(relerr_binned.*w);
%% ENERGY BIN CENTERS
    erplot=er(1:length(er)-1)+0.5*diff(er);
end %function

%% OLD CODE

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