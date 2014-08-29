function onsettime = process_beams(basedir,simInputFN)
% PROCESS A SET OF run_beams.sh Transcar OUTPUT FILES 
% SAVE THE .eps PLOTS AND .mat FILES
clc
%% setup paths
funcPaths = {'../../matlabshared'};
addpath(funcPaths{:}) % enable common functions
if nargin<1, basedir = '../../AT2'; end
if nargin<2, simInputFN = '../../dir.transcar.server/BT_E1E2prev'; end
writePlots=true;
writeMat=true;
plotAllSimTime=true;
%% read simulation configuration
E1E2prev = csvread(simInputFN); % the reason that E1 and E2 are specified seemingly redundantly is for simplicity of parallel processing from Bash with GNU parallel.
energies(:,1) = E1E2prev(:,1);
idlist = beamDirLister(basedir,energies);
%% run through the list
onsettime = process_transcar(idlist,writePlots,writeMat,plotAllSimTime); 
%display(['Computed onset time: ',num2str(onsettime),' hours.'])
%% SAVE SOME METADATA
if writeMat
    outputFN = [basedir,'/ids.mat']; 
    display(['process_beams: writing to ',outputFN])
    save(outputFN,'idlist','energies','onsettime')
end
%% cleanup
rmpath(funcPaths{:})
end %function