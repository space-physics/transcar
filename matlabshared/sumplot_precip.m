function [fig2h, onsettime]=sumplot_precip(datadir,timingfile,excfile,plotAllSimTime) %#ok<STOUT>

%
%function fig2h=sumplot_precip(datadir,timingfile,excfile)
%

%ENERGY TIME SPECTROGRAM OF INPUT PRECIPITATION
timeFN=[datadir,'/',timingfile];
display(['sumplot_precip: loading ',timeFN])
load(timeFN,'onsettime','maxtime')

excFN = [datadir,'/',excfile];
display(['sumplot_precip: loading ',excFN])
load(excFN,'timeop','e*','fluxdown*')
%% DEFINE Time Of Interest
if plotAllSimTime 
    time = timeop;
    tInds = 1:length(time); %take all values 
else
    tInds=find(timeop>=onsettime & timeop<maxtime);
    time=timeop(tInds);
end
%% HOW MANY FLUX TYPES DO WE HAVE IN THE SIM
it=1; ntypes=0;
while it && ntypes<10
    ntypes=ntypes+1;
    arrname=['e',int2str(ntypes)];
    it=exist(arrname,'var');
end
ntypes=ntypes-1;
%% WHAT ARRAY IS THE TOI IN?
scount=0; it=0;
while scount<=tInds(1) && it<=ntypes
    it=it+1;
    stime=size(eval(['e',int2str(it)]),2);
    scount=scount+stime;
end
arrn=it;

%WHAT INDEX INTO THE ARRAY IS IT?
sprev=0;
if arrn>1
  for k=1:arrn-1
    sprev=sprev+size(eval(['e',int2str(arrn-1)]),2);
  end
else
  sprev=0;
end

energy=eval(['e',int2str(arrn),'(:,1)']);
nflux=eval(['fluxdown',int2str(arrn)])/2/pi;
if arrn>1       %TODO: M.Z.: not entirely sure why code breaks if this is not included...
    nflux=nflux(:,tInds-sprev+1);
else
    nflux=nflux(:,tInds);
end

%% PLOT DIFFERENTIAL ENERGY FLUX VS. ENERGY AND TIME
fig2h=E_t_spectrogram(time,energy,nflux,onsettime,datadir);
end