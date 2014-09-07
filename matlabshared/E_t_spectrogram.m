function [figh,xax]=E_t_spectrogram(tset,er,phimem,onsettime,datadir)

figh=figure(76);
cp = get(figh,'pos');
set(figh,'pos',[cp(1), cp(2), 950, 540])
%set(gcf,'PaperPosition',[0,0,8.5,6.5]);
 
%% ARE WE USING ENERGY BINS?
[le,lt]=size(phimem);
if length(er)~=le
  erplot=er(1:length(er)-1)+0.5*diff(er);             %energy bin centers
else
  erplot=er;
end
%% INTERPOLATE ONTO A REGULAR GRID FOR PLOTTING
espec=logspace(log10(50),log10(1.5e4),20);
for k=1:lt
  phiEmem(:,k)=erplot(:).*phimem(:,k);              %compute differential energy flux
  phiEmemplot(:,k)=interpolateTRANSCAR(phiEmem(:,k),erplot,espec,'loglin','loglin');
end
%% FIGURE OUT HOW TO BEST DISPLAY TIME AXIS
nlab=6;                                                 %nominal number of labels on time axis
xlims=[tset(1), tset(lt)];
dx=(xlims(2)-xlims(1))/nlab;
if dx<1/60                                              %labels need seconds to be unambiguous
    dtfmt='HH:MM:SS';
else                                                    %including seconds is pointless
    dtfmt='HH:MM';
end
simdate = zeros(lt,6); %datevec for datenum
simdate(:,4) = tset;
xax=datenum(simdate);
xOnset= datenum([0 0 0 onsettime 0 0]);
%% COMPUTE TOTAL ENERGY FLUX
ermat=repmat(erplot(:),1,lt);
[totalEflux,dE]=intrap(phiEmem,ermat);
totalEflux=totalEflux(le-1,:)*1.602e-15;           %conversion from Ev/cm2/s/sr to W/m2/sr
totalEflux=totalEflux*2*pi;                        %integration over downward pitch angles
totalEflux=totalEflux*1e3;                          %W to mW
%% IMAGE MAP OF DIFFERENTIAL ENERGY FLUX
%FS=16;
subtightplot(2,1,1)
imagesc(xax,log10(espec),log10(phiEmemplot))
axis xy
%datetick('x',fmt,'keepticks','keeplimits');
datetick(gca,'x',dtfmt,'keeplimits')
%set(gca,'FontSize',12)
%caxis([6.5,8])
%set(gca,'FontSize',FS)
c=colorbar;
%set(c,'FontSize',FS)
ylabel(c,'$\log_{10}{E{\cdot}\phi}$ (cm$^{-2}$ s$^{-1}$ sr$^{-1}$)','interpreter','latex')
%xlabel('UT')
ylabel('log_{10} energy (eV)')
%plot onset time line
line([xOnset,xOnset],[0,100],'color','w','linestyle','--')
title([datadir,' energy, energy flux'])


%TOTAL ENERGY FLUX VS. TIME
subtightplot(2,1,2)
plot(xax,totalEflux,'LineWidth',2)
%datetick('x',fmt,'keepticks','keeplimits');
datetick(gca,'x',dtfmt,'keeplimits')
%set(gca,'FontSize',16)
axis('tight')
%xlabel('UT')
ylabel('energy flux (mW/m^2)')
%plot onset time line
line([xOnset,xOnset],[0,1000],'color','k','linestyle','--','linewidth',1)
end
