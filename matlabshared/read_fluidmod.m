function [ne,vi,Ti,Te,time] = read_fluidmod(...
                        datadir,fluidfile,outfile,writeMat) %#ok<STOUT>
% this function writes plasma.mat based on input file transcar_output
%
% the n1 etc. variable comes from read_tra.m -- a scary script indeed -- M.H.
%
% example:
% [ne,vi,Ti,Te,time] = read_fluidmod('~/transcar/2014-04branch/matt2013local/beam3279.5/','transcar_output')
data_type='tra';
data_path=datadir;
file=fluidfile;
%% FIXME: yuck
run('giveno') 
%% manipulate to get ISR parameters
display('READ_FLUIDMOD: Computing ISR-like parameters')
ne=ne_time; %#ok<*NASGU>
vi=(n1_time.*v1_time + n2_time.*v2_time + n3_time.*v3_time + nm_time.*vm_time)./ne_time;
Tipar=(n1_time.*t1p_time+n2_time.*t2p_time+n3_time.*t3p_time+nm_time.*tmp_time)./ne_time;
Tiperp=(n1_time.*t1t_time+n2_time.*t2t_time+n3_time.*t3t_time+nm_time.*tmt_time)./ne_time;
Ti=1/3*Tipar+2/3*Tiperp;
Te=te_time;
comp_ratio=n1_time./ne_time;
%% write to disk
if writeMat
writeToFN=[datadir,filesep,outfile];
display(['READ_FLUIDMOD: Saving multi-fluid plasma parameters to ',writeToFN])
save(writeToFN,'n*','v*','t*','q*','T*','N*','time','z*',...
               'dipangle','chi_time','comp_ratio','Bmag','E*')
end %if
end %function
