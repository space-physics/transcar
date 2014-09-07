isrz=100:10:310;
isrneraw=[.5 1.4 2.4 3.8 4 3.9 3.8 3.5 3.3 3.4 3.7 3.8 3.8 3.85 3.85 3.5 3.25 2.6 2.4 2.3 2.2 2.1]*10^11;

load barker.mat;            

tind=min(find(tbark>=1+3/60));
te_tiavg=mean(te_ti(:,tind:tind+5),2);            %use average Te/Ti to make up for the fact that the scanning took 2 minutes
te_tibark=interpolate(te_tiavg,zbark,isrz,'lin','lin');
nesimavg=mean(nesim(:,tind:tind+5),2);

correction_fact=(te_tibark+1)/2;        %multiplicative correction factor
nereal=correction_fact.*isrneraw;       %correct ne

FS=16;
plot([isrneraw(:),nereal(:)],isrz(:),'LineWidth',2);
set(gca,'FontSize',FS)%,'FontWeight','bold');
ax=axis;
axis([ax(1:3),375])
xlabel('n_e (m^{-3})');
ylabel('altitude (km)');
title('Barker-coded data from 19 Nov. 2001');
legend('n_{e} with T_e/T_i=1','n_{e} using estimated T_e/T_i','Location','NorthWest');
print -depsc barker_corr.eps