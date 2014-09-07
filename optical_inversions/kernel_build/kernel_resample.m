%ADJUSTABLE SCRIPT PARAMETERS
datelab='11Feb2002';
phiE=7.0114e10;                 %This is set by the main beams.sh script that runs the beam simulations in TRANSCAR
mode='const_energy';            %hold total energy flux constant for each beam
%mode='even'                    %hold differential number flux constant for each beam
nbins=8;



%LOAD FULL CHARACTERISTIC MATRIX FOR A SPECIFIED DATE
load(['kernel_',datelab,'.mat']);



%DIVIDE ENERGIES RANGES INTO A SPECIFIED NUMBER OF BINS
le=length(energies)-1;
nperbin=floor(le/nbins);    %number of bins in the original dist. that go into each resampled bin
nused=nperbin*nbins;        %how many do we use?
ntrunc=le-nused;            %how many energy bins are we hacking off of the end of the dist.
phiEbin=nperbin*phiE;       %total energy flux per (new) bin



%RESAMPLE THE CHARACTERISTIC MATRICES IN ENERGY ACCORDING TO NBINS DEFINED
%ABOVE
for k=1:nbins
    %resample energy
    er(k)=energies((k-1)*nperbin+1);
    emin=er(k); emax=energies(k*nperbin+1);
    debin=emax-emin;

    %perform the appropriate scaling and adding to form new matrix
    bchrresamp(:,k)=zeros(size(bchr(:,1)));
    for l=1:nperbin
        ncol=(k-1)*nperbin+l;
        if(k==1)
            debin=emax-emin;
            phiNbin(k)=phiEbin/(0.5*(emin+emax)*debin);
        else
            if(strcmp(mode,'even'))
                phiNbin(k)=phiNbin(1);
            else
                phiNbin(k)=phiEbin/(0.5*(emin+emax)*debin);
            end
        end
        ratio=phiNbin(k)/phiN(ncol);

        bchrresamp(:,k)=bchrresamp(:,k)+ratio*bchr(:,ncol);
    end
end
er(nbins+1)=emax;               %Add in the last energy point so we can difference later



%DEFINE SIZE OF SYSTEM
[nlambda,nphi]=size(bchrresamp);



%SAVE RESULTS IN AN APPROPRIATELY NAMED .MAT FILE
save(['kernel_resampled_',datelab,'.mat'],'bchrresamp','er','phiNbin','phiEbin','nphi','nlambda');