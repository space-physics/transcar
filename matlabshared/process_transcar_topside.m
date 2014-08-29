%MAIN PROCESSING FILE FOR TRANSCAR DATA, idlist must be defined or else ''
%will be used.  idlist must also contain full path to simulation directory.
%close all; clc;



%RAW OUTPUT FILES OF TRANSCAR
sim_inputfile='DATCAR';
sim_fluidfile='transcar_output';
sim_excfile='emissions.dat';



%.MAT OUTPUT FILES FOR ANALYSIS CODE (THIS SCRIPT)
timingfile='dat.mat';
plasmafile='plasma.mat';
excfile='excrates.mat';
emfile='aurora.mat';
confile='conductivities.mat';



%IF THERE IS NO LIST OF INPUT SIMULATIONS JUST ASSUME THAT WE WANT TO
%PROCESS WHAT IS IN PWD
if ~exist('idlist')
    idlist={[pwd,'']};
end



%READ FILES AND GENERATE MATLAB VARIABLES LOOP OVER LIST OF SIMULATION IDs
tic;
for l=1:length(idlist)
    datadir=idlist{l};



    %PROCESS INPUT FILES FOR THIS SIMULATION
    read_input(datadir,sim_inputfile,timingfile);                       %read in simulation timing info and save to .mat file
    read_fluidmod(datadir,sim_fluidfile,plasmafile);                    %read fluid code output and save to .mat file
    read_excrates(datadir,sim_excfile,excfile);                         %read kinetic code output
    calc_aurora(datadir,excfile,emfile);                                %process kinetic output into auroral brightnesses
    calc_conductivity(datadir,plasmafile,confile);                      %process fluid code output into conductivities and Joule dissipation rates


    %CREATE A DIRECTORY FOR PLOTS
    if ~ (exist([datadir,'/figures']) == 7)
        command=['mkdir ',datadir,'/figures'];
        system(command);
    end



    %DID THE SIMULATION COMPLETE PROPERLY?
    load([datadir,'/',timingfile],'endtime');
    load([datadir,'/',plasmafile],'time');
    dt=mean(diff(time));
    if(time(length(time)) < endtime-2*dt)
      fprintf('\nPROCESS_TRANSCAR.M --> There appears to be a problem with these simulation files!');
      break;
    end
    

    
    %SUMMARY PLOT OF PLASMA PARAMETERS AND OPTICAL EMISSIONS
    alt_lims=[90, 2000];
    figh=sumplot_transcar(datadir,timingfile,plasmafile,emfile,alt_lims);                       %an 8.5x11 summary plot
    saveas(figh,[datadir,'/figures/summary_plot.fig'],'fig');                                   %save a matlab figure file
    print(figh,'-depsc',[datadir,'/figures/summary_plot.eps']);                                 %make a postscript file

    fig2h=sumplot_precip(datadir,timingfile,excfile);
    saveas(fig2h,[datadir,'/figures/spectrogram.fig'],'fig');
    print(fig2h,'-depsc',[datadir,'/figures/spectrogram.eps']);
    
    fig3h=sumplot_conduct(datadir,timingfile,confile);
    saveas(fig3h,[datadir,'/figures/conductivities.fig'],'fig');
    print(fig3h,'-depsc',[datadir,'/figures/conductivities.eps']);    
    close all;
    
    fprintf('PROCESS_TRANSCAR.M --> .mat files and plots generated for simulation with ID:  %s\n',datadir);
end
toc
