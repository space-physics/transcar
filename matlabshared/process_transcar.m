function onsettime = process_transcar(idlist,writePlots,...
                     writeMat,plotAllSimTime)
%MAIN PROCESSING FILE FOR TRANSCAR DATA, idlist must be defined or else ''
%will be used.  idlist must also contain full path to simulation output data directory.
%% RAW OUTPUT FILES OF TRANSCAR
sim_inputfile='DATCAR';
sim_fluidfile='transcar_output';
sim_excfile='emissions.dat';
%% MAT OUTPUT FILES FOR ANALYSIS CODE (THIS SCRIPT)
timingfile='dat.mat';
plasmafile='plasma.mat';
excfile='excrates.mat';
emfile='aurora.mat';
confile='conductivities.mat';

%IF THERE IS NO LIST OF INPUT SIMULATIONS JUST PROCESS  PWD
if nargin<1
    idlist={pwd};
else
    if ischar(idlist), idlist={idlist}; end
end
if nargin<2, writePlots = false; end
%% READ FILES AND GENERATE MATLAB VARIABLES LOOP OVER LIST OF SIMULATION IDs
for l=1:length(idlist)
    datadir=idlist{l};
    %check for existance, skip if not there
    if ~exist(datadir,'dir')
        display(['Skipping non-existant directory ',datadir])
        continue
    end

    %% PROCESS INPUT FILES FOR THIS SIMULATION
    %read in simulation timing info and save to .mat file
    [endtime] = read_input(datadir,sim_inputfile,timingfile);     
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %read kinetic code output
    read_excrates(datadir,sim_excfile,excfile)
     %kinetic output --> auroral VER
    calc_aurora(datadir,excfile,emfile)           
    
    
        %read fluid code output and save to .mat file
        [ne,vi,Ti,Te,time] = read_fluidmod(datadir,sim_fluidfile,...
                                   plasmafile,writeMat); 
        %process fluid code output into conductivities and Joule dissipation rates
        calc_conductivity(datadir,plasmafile,confile);    
        
        %DID THE SIMULATION COMPLETE PROPERLY?
        dt=mean(diff(time));
        if(time(length(time)) < endtime-2*dt)
         warning('PROCESS_TRANSCAR: There appears to be a problem with these simulation files!')
        end %if
    
    
    %CREATE A DIRECTORY FOR PLOTS
    plotDir=[datadir,filesep,'figures'];
    if ~exist(plotDir,'dir'), 
        display(['creating plot directory ',plotDir])
        mkdir(plotDir);  
    end
%% SUMMARY PLOT OF PLASMA PARAMETERS AND OPTICAL EMISSIONS
    alt_lims=[90, 1000]; %[km]
 
    [fig2h, onsettime] = sumplot_precip(datadir,timingfile,excfile,plotAllSimTime);
    [figh,hfigspec] = sumplot_transcar(datadir,timingfile,plasmafile,emfile,...
                      alt_lims,onsettime,plotAllSimTime); %an 8.5x11 summary plot
    
    fig3h=sumplot_conduct(datadir,timingfile,confile,onsettime,plotAllSimTime);

if writePlots
    plotDir = [datadir,'/figures/'];
    display(['writing EPS plots to ',plotDir])
    print(figh, '-depsc2',[plotDir,'summary_plot.eps']); 
    print(fig2h,'-depsc2',[plotDir,'spectrogram.eps']);
    print(fig3h,'-depsc2',[plotDir,'conductivities.eps']); 
    print(hfigspec,'-depsc2',[plotDir,'OpticalSpecrum.eps']); 
end %if
    
    display(['PROCESS_TRANSCAR: wrote .mat for ID: ',datadir])
    display(['Precipitation onset time: ',num2str(onsettime),' hours.'])
end %for

end %function
