function []=plasmaplots(testname)

% This function returns the plots of the electron heat flux and the electron
% parallel and perpendicular temperatures for a specified simulation

% The argument for this function is the string of the name of the simulation
% that the user is interested in looking at. For example if the directory
% that contains the information that the user is interested in is: models/simulations/test.SERSIO the
% test.SERSIO would be the argument of this function.
% Tyler Parsotan

%cd(['/models/simulations/',testname]) %need to change directories

%load data from simulation
load(['/models/simulations/',testname,'/',plasma.mat])

%% plot data
%plot electron heat flux
figure(660),clf(660)
pcolor(time, z, qe_time);
shading flat
title('Electron Heat Flux');

%plot electron temperature parallel
figure(661),clf(661)
pcolor(time, z, tep_time)
shading flat
title('Electron Parallel Temperature');

%plot electron perpendicular temperature
figure(662),clf(662)
pcolor(time, z, tet_time)
shading flat
title('Electron Perpendicular Temperature');





