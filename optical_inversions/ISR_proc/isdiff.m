%CALCULATE CORRECTIONS TO BARKER-CODED DATA
if ~exist('neplotbcorr')
    isplot_barker_corr_wsim;
end


%SMOOTH THE BARKER PROFILES TO MAKE IT EASIER TO FIND MAXIMA
mawindow=15;
for k=1:length(tplotb)
   neplotbsmooth(:,k)=smooth(neplotb(:,k),mawindow);
   neplotbcorrsmooth(:,k)=smooth(neplotbcorr(:,k),mawindow);
end
for k=1:length(tplotf)
   neplotfsmooth(:,k)=smooth(neplotf(:,k),4);
end



%FIND NMF2 AND HMF2 FOR EACH SIMULATION
isearch=find(zplot>=175 & zplot<=400);
for k=1:length(tplotb)
   [nmf2b(k),imaxb]=max(neplotbsmooth(isearch,k));
   hmf2b(k)=zplot(isearch(imaxb));
end

for k=1:length(tplotb)
   [nmf2c(k),imaxc]=max(neplotbcorrsmooth(isearch,k));
   hmf2c(k)=zplot(isearch(imaxc));
end

for k=1:length(tplotf)
   [nmf2f(k),imaxf]=max(neplotfsmooth(isearch,k));
   hmf2f(k)=zplot(isearch(imaxf));
end



%INTERPOLATE FITTED DATA ONTO BARKER TIME GRID
nmf2fplot=interpolate(nmf2f',tplotf',tplotb','lin','lin')
errb=abs(100*(10.^nmf2b-10.^nmf2fplot)./10.^nmf2fplot);
errbcorr=abs(100*(10.^nmf2c-10.^nmf2fplot)./10.^nmf2fplot);
itmean=min(find(tplotb>=1+2/60)):length(tplotb);
avgerrb=mean(errb(itmean)); avgerrbcorr=mean(errbcorr(itmean));



%PLOTS
figure;

subplot(211);
plot(tplotb,10.^nmf2b,tplotf,10.^nmf2f,tplotb,10.^nmf2c,'LineWidth',2);
axtemp=axis;
ax=[min(tplotb),max(tplotb),2.5e11,8e11];
axis(ax);
set(gca,'FontSize',12);
l=legend('raw n_e','fitted n_e','T_e/T_i adj. n_e','Location','NorthWest','Orientation','Horizontal');
set(l,'FontSize',12);
xlabel('UT (hrs)');
ylabel('N_mF_2');
set(gca,'XTick',xloc,'XTickLabel',xlab);

subplot(212);
plot(tplotb,errb,tplotb,errbcorr,'LineWidth',2);
axtemp=axis;
ax=[min(tplotb),max(tplotb),0,70];
axis(ax);
set(gca,'FontSize',12);
l=legend(sprintf('raw n_e.  Avg:  %3.1f %%',avgerrb),sprintf('T_e/T_i adj. n_e.  Avg:  %3.1f %%',avgerrbcorr),'Location','NorthWest','Orientation','Horizontal');
set(l,'FontSize',12);
xlabel('UT (hrs)');
ylabel('% Error in N_mF_2');
set(gca,'XTick',xloc,'XTickLabel',xlab);

print -depsc error_barker.eps