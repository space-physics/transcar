function hf = plot_spectrum(datadir,emfile,midevent,mintime,maxtime,doISRinterp)

%% LOAD OPTICAL DATA
optFN = [datadir,'/',emfile];
display(['sumplot_transcar: Loading ',optFN])
load(optFN,'timeop','blambda','lambda*')
%% WHAT TIME INDEX?
%tind=find(timeop>=midevent, 1 );
tInds= timeop>=mintime & timeop<=maxtime;
%lt=length(tInds);
bavg=mean(blambda(:,tInds),2);
%plot sizing parameters
FS=8;
LW=2;
%%
%BUILD LEGEND (THIS IS A BIT OF A HACK TO GET AROUND MATLAB'S STUPID LEGEND CODE)
%subtightplot(5,1,5); %pbaspect([8.5 11/6.5 1]);
hf = figure(681);
clf(681)
cp = get(681,'pos');
set(681,'pos',[cp(1),cp(2),775,420])
line(557.7,1e0,'color','g')
plot(557.7,1e0,'r',557.7,1e0,'b',557.7,1e0,'y',557.7,1e0,'m',557.7,1e0,'c','LineWidth',LW);
legend('O meta-stable','O prompt','N_2^+ 1N','N_2^+ Meinel','N_2 2P','N_2 1P',...
        'Location','NorthWest','Orientation','Horizontal');
legend('boxoff')
title(['Auroral emission spectrum,',datadir,'  interp=',int2str(doISRinterp)])
%% BUILD LINES FOR PLOTTING
linebase=0.5*min(bavg);
blines=[bavg,linebase*ones(size(bavg))];
%% IDENTIFY EMISSION GROUPS AND PLOT ON SYNTHETIC SPECTROGRAPH
%metastable oxygen lines
hold on;
for k=1:length(lambdams)
    ilk=find(lambda==lambdams(k));
    plot([lambda(ilk),lambda(ilk)],log10(blines(ilk,:)),'g','LineWidth',LW)
end

%prompt oxygen emissions
for k=1:length(lambdapa)
    ilk=find(lambda==lambdapa(k));
    plot([lambda(ilk),lambda(ilk)],log10(blines(ilk,:)),'r','LineWidth',LW)
end

%1NG
for k=1:length(lambda1NGvec)
    ilk=find(lambda==lambda1NGvec(k));
    plot([lambda(ilk),lambda(ilk)],log10(blines(ilk,:)),'b','LineWidth',LW)
end

%Meinel
for k=1:length(lambdaMeinvec)
    ilk=find(lambda==lambdaMeinvec(k));
    plot([lambda(ilk),lambda(ilk)],log10(blines(ilk,:)),'y','LineWidth',LW)
end

%2PG
for k=1:length(lambda2PGvec)
    ilk=find(lambda==lambda2PGvec(k));
    plot([lambda(ilk),lambda(ilk)],log10(blines(ilk,:)),'m','LineWidth',LW)
end

%1PG
for k=1:length(lambda1PGvec)
    ilk=find(lambda==lambda1PGvec(k));
    plot([lambda(ilk),lambda(ilk)],log10(blines(ilk,:)),'c','LineWidth',LW)
end
%hold off
xlabel('wavelength (nm)')%,'FontSize',FS)
ylabel('log_{10} Brightness (R)')%,'FontSize',FS)
%title('Simulated auroral emissions','FontSize',FS)
%set(gca,'FontSize',FS)
axis([200,1000,-1,6])
%text(225,3.5,'(e)','FontSize',20,'Color',[0 0 0],'FontWeight','bold')

end %function