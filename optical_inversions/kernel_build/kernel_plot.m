%LOAD DATA FROM CHARACTERISTIC MATRIX BUILDER
datelab='17Feb2001';
load(['kernel_',datelab,'.mat']);



%NORMALIZE QUANTITIES TO SHOW ON SAME PLOT
normlist=[1,4,5,7,8];            %band emissions have same characteristics so we dont' need to plot everything
for il=1:length(normlist);
   bchrnorm(il,:)=bchr(normlist(il),:)/max(bchr(normlist(il),:));
end
iz=min(find(z>=750));
phiichrnorm=phiichr(iz,:)/max(phiichr(iz,:));



%CREATE DIRECTORY FOR OUTPUT PLOTS
if ~ (exist(['./plots']) == 7)
    command=['mkdir ./plots'];
    system(command);
end



%PLOT CHARACTERISTIC BRIGHTNESSES
LW=2;
FS=16;
eplot=energies(1:length(energies)-1)+0.5*diff(energies);
figure;
%set(0,'DefaultAxesLineStyleOrder','-|--|^-|-.|*-|o-');
semilogx(eplot,[bchrnorm',phiichrnorm(:)],'LineWidth',LW);
set(gca,'FontSize',FS)%,'FontWeight','bold');
legend('O^+(^2P) 732nm','N_2 1PG','O(3p^5P) 777.4nm','N_2^+ Meinel','O(3p^3P) 844.6nm','\phi_i(750 km)');
legend boxoff;
xlabel('energy (eV)');
ylabel('fraction of max value');
axis([5e1,1e4,-0.1,1.5]);
saveas(gcf,'./plots/charresp.fig','fig');
print('-depsc','./plots/charresp.eps');



%PLOT CHARACTERISTIC PLASMA RESPONSE
i1=1;
i2=max(find(z<1000));

figure;
h=pcolor(log10(eplot),z(i1:i2),phiichr(i1:i2,:));
set(h,'LineStyle','none');
colorbar;
h=gca;
set(h,'FontSize',FS)%,'XTick',[23,23.2,23.4])
ylabel('altitude (km)');
xlabel('log_{10} energy (eV)');
title('\Delta \phi_i (m^{-2} s^{-1}) vs. beam energy and alt.')
colormap('gray')
saveas(gcf,'./plots/chrnvz.fig','fig');
print -depsc ./plots/chrnvz.eps

figure;
h=pcolor(log10(eplot),z(i1:i2),techr(i1:i2,:));
set(h,'LineStyle','none');
colorbar;
h=gca;
set(h,'FontSize',FS)%,'XTick',[23,23.2,23.4])
ylabel('altitude (km)');
xlabel('log_{10} energy (eV)');
title('\Delta T_e (K) vs. beam energy and alt.')
colormap('gray')
saveas(gcf,'./plots/chrtez.fig','fig');
print -depsc ./plots/chrtez.eps

figure;
h=pcolor(log10(eplot),z(i1:i2),nechr(i1:i2,:));
set(h,'LineStyle','none');
colorbar;
h=gca;
set(h,'FontSize',FS)%,'XTick',[23,23.2,23.4])
ylabel('altitude (km)');
xlabel('log_{10} energy (eV)');
title('\Delta n_e (m^{-3}) vs. beam energy and alt.')
colormap('gray')
saveas(gcf,'./plots/chrnez.fig','fig');
print -depsc ./plots/chrnez.eps

figure;
h=pcolor(log10(eplot),z(i1:i2),vichr(i1:i2,:));
set(h,'LineStyle','none');
colorbar;
h=gca;
set(h,'FontSize',FS)%,'XTick',[23,23.2,23.4])
ylabel('altitude (km)');
xlabel('log_{10} energy (eV)');
title('\Delta v_i (m^{-1} s^{-1}) vs. beam energy and alt.')
colormap('gray')
saveas(gcf,'./plots/chrviz.fig','fig');
print -depsc ./plots/chrviz.eps