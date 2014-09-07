function [endtime] = read_input(datadir,datfile,outfile)
%upgraded by Michael Hirsch
%
%function read_input(datadir,datfile,outfile)
%

%OPEN DATCAR FILE AND READ THE TIMING INFORMATION FOR THE PRECIPITATION
inFN = [datadir,'/dir.input/',datfile];
display(['Reading ',inFN])
fid=fopen(inFN,'r');
%{
for iline=1:29
    line=fgetl(fid);
    switch(iline)                       %only bother with timing info for now
        case 4, timestr{1}=line;
        case 6, timestr{2}=line;
        case 7, timestr{3}=line;
        case 28,timestr{4}=line;
        case 29,timestr{5}=line;
        otherwise %do nothing
    end
end
fclose(fid);

%SORT OUT THE TIMES FROM THE CHARACTERS
for k=1:length(timestr)
    is=1; curr=timestr{k};
    while curr(is)>='0' && curr(is)<='9'
      is=is+1;
    end 
    timing(k)=str2num(curr(1:is-1))/3600;        %index 1: dt, 2: starttime, 3:duration, 4:precstart, 5:precend
end
%}

%% upgrade (by MHirsch)
NlinesDATCAR = 29; 
timestr = textscan(fid,'%s %s',NlinesDATCAR,'MultipleDelimsAsOne',true,'Delimiter','\t');
fclose(fid);
timing = str2num(strjoin(timestr{1}([4;6;7;28;29])'))./3600; %#ok<ST2NM>
%timing(1) dt, (2) starttime, (3) duration, (4) precStart, (5) precEnd
%%

%DEFINE INDICES TIMES OF INTEREST FOR PLOTS OF ISR.M AND SPECTROGRAPH.M
minz=90;           %min alt [km?]
maxz=1000;         %max alt [km?]

if timing(4)<timing(2)  %precipitation start < UT starttime simulation
   timing(4:5)=timing(4:5)+24; %shift precipitation forward in time +24 hours
end

duree=timing(5)-timing(4); %how long precipitation lasted [hours]
mintime=timing(4)-.5*duree; %start at precStart - 0.5*duration
maxtime=timing(5)+.5*duree; %stop at precEnd + 0.5*duration
endtime=timing(2)+timing(3); %UT sim start + sim duration 

% TODO: REMOVE THIS LATER.  THIS IS A KLUDGE TO DEAL WITH THE TRANSCAR
% PRECIPITATION NOT CEASING ISSUE.
maxtime=timing(5);
%

if maxtime > endtime
    maxtime=endtime; %#ok<*NASGU>
end
midevent=timing(4)+.5*duree;
onsettime=timing(4); %

writeToFN= [datadir,filesep,outfile];
display(['READ_INPUT.M: Saving timing info to ',writeToFN])

save(writeToFN,'timing','duree','mintime','maxtime','endtime','midevent','onsettime','minz','maxz');

end



%% EXAMPLE DATCAR FILE
% 2					(kiappel, INTEGER)
% 90kmmaxpt.dat                                input file (initial ionospheric conditions)
% 1.					time step (seconds) (REAL)
% 300					number of seconds between fluid code outputs (seconds, REAL)
% 2001048					date of simulation (YYYYDDD)
% 25200.					UT start time for simulation (seconds)
% 68400.					duration of run (seconds)
% 0					jpreci, precipation type (see table below, INTEGER)
% 67.,-50.62    				location of simulation, geodetic latitude,longitude
% 0.					duree de la convection (en secondes) avant la reference (<=0, on ne suit pas les lignes de champ) (tempsconv_1, REAL)
% 0.					duree de la convection (en secondes) apres la reference (<=0, on ne suit pas les lignes de champ) (tempsconv, REAL)
% 0					intervalle de temps (en secondes) separant deux tubes (step,REAL)
% 300					time step between two calls to kinetic code (REAL)
% 0.					transport (en m/s) induit le long de la ligne de champ (vparaB, REAL)
% 126.8					f10.7 index
% 2.					ap index
% 10.					convection electric field (mV/m)
% 1.					O correction factor (REAL)
% 1.					N2 correction factor (REAL)
% 1.					O2 correction factor (REAL)
% 1.					N correction factor (REAL)
% 1.					H correction factor (REAL)
% -1.e-3					topside electron heat flux (mW/m^2, REAL)
% precinput.dat	                            electron precipitation distribution file (energy (eV) in first column, number flux (cm^-2 s^-1 eV^-1 sr^-1) in second column)
% 1					precint, precipation energy interpolation (see table below, INTEGER)
% 1					precext, precipation energy extrapolation (see table below, INTEGER)
% 3600.					precipitation start time (seconds)
% 4800.					precipitation end time (seconds)
% 
% 
% 
% jpreci = 0 if sun only
%        = 1 if electron precipitation only
%        = 2 if sun + electron precipitation
%        = 3 if proton precipitation only
%        = 4 if proton and electron precipitation
%        = 5 if sun + proton precipitation
%        = 6 if proton, electron precipitation + sun
% 
% precint = 0 log-linear interpolation
% 	= 1 sample and hold interpolation
% 
% precext = 0 log-linear extrapolation
% 	= 1 truncation of distribution past energy range of input file