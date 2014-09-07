%LOAD DIRECTORY NAMES FOR THIS SET OF BEAMS
clear,close all;
addpath ~/matlabshared;
basedir='~/simulations/TRANSCAR/beamsims/';
%datelab='01Mar2011v2';
datelab='10Dec2007';
load([basedir,datelab,'/ids.mat']);
le=numel(idlist);
E=1/2*(energies(1:le)+energies(2:le+1));


%FORM CHARACTERISTIC MATRICES COLUMN BY COLUMN (EACH COLUMN COMES FROM A DIFFERENT SIM)
for ie=1:length(idlist)
    fprintf('\nKERNEL_BUILD.M --> loading files for simulation:  %s', idlist{ie});
    
    %load data for this beam
    load([idlist{ie},'/plasma.mat']);
    load([idlist{ie},'/dat.mat']);
    [VER,z,t,lambda,IDs,types,b,diffnumflux_interp,E_interp]=VER_calc(idlist{ie},'excrates.mat');


    %Define wavelength set that we want to use and toi
    if ie==1
        it=min(find(t>onsettime))+0;               %needs to be +1 for plotting purposes, +0 for inversions
        [lz,lt,ll]=size(VER);
        Mp=zeros(lz,le,ll);                             %for wavelength il, the forward model is M=Mp(:,:,il)
        Mb=zeros(ll,le);
    end
    
    
    %characteristic plasma resp.  +1 b/c of the delay of fluid code
    techr(:,ie)=te_time(:,it+1)-te_time(:,it);
    tichr(:,ie)=t1_time(:,it+1)-t1_time(:,it);
    vichr(:,ie)=v1_time(:,it+1)-v1_time(:,it);
    phiichr(:,ie)=n1_time(:,it+1).*v1_time(:,it+1)-n1_time(:,it).*v1_time(:,it);
    nechr(:,ie)=ne_time(:,it+1)-ne_time(:,it);
    nichr(:,ie)=n1_time(:,it+1)-n1_time(:,it);
    
    
    %characteristic optical response
    tmp=squeeze(VER(:,it,:));
    tmp=reshape(tmp,[lz 1 ll]);
    Mp(:,ie,:)=tmp;
    Mb(:,ie)=b(:,it);
    
    %beam intensity
    load([idlist{ie},'/aurora.mat']);
    precchr(:,ie)=fluxdown1(:,it);
%    precchr(:,ie)=fluxdown1(:,1);                   %Be careful here.  The flux variable depends on how the energy binning is set in the simulation
    phiN(ie)=max(precchr(10:numel(precchr(:,ie)),ie))/2/pi;          %The units here are cm-2 s-1 eV-1 sr-1
end
fprintf('\n');



%SAVE RESULTS IN PWD
save(['kernel_',datelab,'.mat'],'Mp','Mb','lambda','energies','z','phiN','datelab','E','IDs','types');
rmpath ~/matlabshared;
