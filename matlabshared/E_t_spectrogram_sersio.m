function figh=E_t_spectrogram_sersio(tset,er,phimem)

figh=figure(685);
%set(gcf,'PaperPosition',[0,0,8.5,3.5]);
%% ARE WE USING ENERGY BINS?
[le,lt]=size(phimem);
if length(er)~=le
  erplot=er(1:length(er)-1)+0.5*diff(er);             %energy bin centers
else
  erplot=er;
end
%% INTERPOLATE ONTO A REGULAR GRID FOR PLOTTING
espec=logspace(log10(10),log10(1.1e4),50);
for k=1:lt
  phiEmem(:,k)=erplot(:).*phimem(:,k);              %compute differential energy flux
  phiEmemplot(:,k)=interpolateTRANSCAR(phiEmem(:,k),erplot,espec,'loglin','loglin');
end
%% FIGURE OUT HOW TO BEST DISPLAY TIME AXIS
nlab=6;                                                 %nominal number of labels on time axis
xlims=[tset(1), tset(lt)];
dx=(xlims(2)-xlims(1))/nlab;
if dx<1/60                                              %labels need seconds to be unambiguous
    fmt='hms';
else                                                    %including seconds is pointless
    fmt='hm';
end
xloc=xlims(1):dx:xlims(2);
xlab=time2str(xloc,'clock',fmt);
%% IMAGE MAP OF DIFFERENTIAL ENERGY FLUX
imagesc(tset,log10(espec),flipud(log10(phiEmemplot)),[8 8.75]);
h=gca;
c=colorbar;
set(h,'FontSize',12);
set(h,'XTick',xloc); %FIXME: use datetick instead
set(h,'XTickLabel',xlab);
maxE=max(log10(espec));
minE=min(log10(espec));
yloc0=get(h,'YTick');
dyloc=yloc0-round(minE);
newyloc=round(maxE)-dyloc;
set(h,'YTick',sort(newyloc));
set(h,'YTickLabel',fliplr(yloc0));
title('Energy-time spectrogram of e^- precip.');
ylabel(c,'log_{10} E{\cdot}\phi (cm^{-2} s^{-1} sr^{-1})')
xlabel('UT');
ylabel('log_{10} energy (eV)');

end %function