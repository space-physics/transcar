%LOAD DIRECTORY NAMES FOR THIS SET OF BEAMS
datelab='11Feb2002';
load(['~/beamsims/',datelab,'/ids.mat']);
lambdalist=[732; 738.7; 750.5; 775.4; 777.4; 785.3; 808.2; 844.6];
lid={'O^+(^2P)','N_2 1P (5-3)','N_2 1P (4-2)','N_2 1P (2-0)', ...
    'O(3p^5P)','N_2^+ M (2-0)','N_2^+ M (3-1)','O(3p^3P)'};



%FORM CHARACTERISTIC MATRICES COLUMN BY COLUMN (EACH COLUMN COMES FROM A DIFFERENT SIM)
for ie=1:length(idlist)
    fprintf('\nCHARBUILD.M --> loading files for simulation:  %s', idlist{ie});
    
    %load data for this beam
    load([idlist{ie},'/plasma.mat']);
    load([idlist{ie},'/aurora.mat']);
    load([idlist{ie},'/dat.mat']);
    
    %Define wavelength set that we want to use and toi
    if ie==1
        llam=length(lambdalist);
        for il=1:llam
            ilambda(il)=find(lambda==lambdalist(il));
        end
        it=min(find(timeop>onsettime))+0;               %needs to be +1 for plotting purposes, +0 for inversions
        verchr=zeros(llam,size(zop,1),length(idlist));
    end
    
    %characteristic plasma resp.  +1 b/c of the delay of fluid code
    techr(:,ie)=te_time(:,it+1)-te_time(:,it);
    tichr(:,ie)=t1_time(:,it+1)-t1_time(:,it);
    vichr(:,ie)=v1_time(:,it+1)-v1_time(:,it);
    phiichr(:,ie)=n1_time(:,it+1).*v1_time(:,it+1)-n1_time(:,it).*v1_time(:,it);
    nechr(:,ie)=ne_time(:,it+1)-ne_time(:,it);
    nichr(:,ie)=n1_time(:,it+1)-n1_time(:,it);
    
    %characteristic optical response
    for il=1:llam
      verchr(il,:,ie)=permute(plambda(:,it,ilambda(il)),[1 3 2]);
      bchr(il,ie)=blambda(ilambda(il),it);
    end
    
    %beam intensity
    precchr(:,ie)=fluxdown2(:,1);                   %Be careful here.  The flux variable depends on how the energy binning is set in the simulation
    phiN(ie)=max(precchr(10:150,ie))/2/pi;          %The units here are cm-2 s-1 eV-1 sr-1
end
fprintf('\n');



%SAVE RESULTS IN PWD
save(['kernel_',datelab,'.mat'],'*chr*','ilambda','lambdalist','lid','energies','z','zop','phiN','datelab');