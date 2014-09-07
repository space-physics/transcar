function sumploth=sumplot_conduct(datadir,timingfile,confile,onsettime,plotAllSimTime)

%% loading data 
timingFN = [datadir,'/',timingfile];
display(['sumplot_conduct: loading ',timingFN])
load(timingFN,'maxtime','mintime','maxtime')


confFN = [datadir,'/',confile];
display(['sumplot_conduct: loading ',confFN])
load(confFN,'SigmaP','SigmaH','time','EmagVm','Joulediss')

sumploth=figure(686);
cp = get(sumploth,'pos');
set(sumploth,'pos',[cp(1), cp(2), 950, 520])
%set(gcf,'PaperPosition',[0,0,8.5,11]);
%% CORRECT TIME VARIABLE FOR PLOTTING PURPOSES
if plotAllSimTime
    tInds = 1:length(time); %take all values
else %restrict to around precipitation time 
    tspan=maxtime-mintime;
    tpad=1/6*tspan;                 %show quiescent ionosphere before event is switched on
    tInds=find(time>=mintime+tpad & time<=maxtime);
end

lt = length(tInds);
%% FIGURE OUT HOW TO BEST DISPLAY TIME AXIS
nlab=6;                                                 %nominal number of labels on time axis
xlims=[time(tInds(1)), time(tInds(end))];
dx=(xlims(2)-xlims(1))/nlab;
if dx<1/60                                              %labels need seconds to be unambiguous
    dtfmt='HH:MM:SS';
else                                                    %including seconds is pointless
    dtfmt='HH:MM';
end
% xloc=xlims(1):dx:xlims(2);
% xlab=time2str(xloc,'clock',fmt);
simdate = zeros(lt,6); %datevec for datenum
simdate(:,4) = time(tInds);
xax=datenum(simdate);
xOnset= datenum([0 0 0 onsettime 0 0]);
%% PLOT CONDUCTIVITIES
h=subtightplot(3,1,1);
plot(xax,SigmaP(tInds),...
     xax,SigmaH(tInds),...
                'LineWidth',2)
%set(h,'FontSize',12)
%set(gca,'XTick',xloc,'XTickLabel',xlab);
%ax=axis;
%xloc=time(tInds);
%axis([min(xloc), max(xloc), ax(3:4)])
axis('tight')
datetick(h,'x',dtfmt)
legend('\Sigma_P','\Sigma_H','location','northwest')
%xlabel('UT')
ylabel('conductivity (mhos)')
%plot onset time line
line([xOnset,xOnset],[0,1000],'color','k','linestyle','--')
title([datadir, ' conductivities, E field, Dissipation'])
%% PLOT ELECTRIC FIELD
h = subtightplot(3,1,2);
plot(xax,EmagVm(tInds)*1e3,'LineWidth',2);
%set(h,'FontSize',12)
%set(gca,'XTick',xloc,'XTickLabel',xlab);
%ax=axis;
%axis([min(xloc), max(xloc), ax(3:4)])
axis('tight')
datetick(h,'x',dtfmt)
%xlabel('UT')
ylabel('|E_\perp| (mV/m)')
%plot onset time line
line([xOnset,xOnset],[0,1000],'color','k','linestyle','--')
%% PLOT JOULE DISSIPATION
h = subtightplot(3,1,3);
plot(xax,Joulediss(tInds)*1e3,'LineWidth',2);
%set(h,'FontSize',12);
%set(gca,'XTick',xloc,'XTickLabel',xlab);
%ax=axis;
%axis([min(xloc), max(xloc), ax(3:4)])
axis('tight')
datetick(h,'x',dtfmt)
%xlabel('UT')
ylabel('Joule dissipation (mW/m^2)')

%plot onset time line
line([xOnset,xOnset],[0,1000],'color','k','linestyle','--')

%=---------------------


end %function
