function VERkernel_plot()

%for Hanna Dahlgren's paper (?)

datelab='01Mar2011v2';
%datelab='10Dec2007';
load(['kernel_',datelab,'.mat']);

%my standard NIR set from [Zettergren et al, 2008]
lambdalist=[732; 738.7; 750.5; 775.4; 777.4; 785.3; 808.2; 844.6];  

%sampling of different bright N2 bands
%lambdalist=[844.6,427.8,785.3,405.8,775.4]     

ll=numel(lambdalist);

for il=1:ll
    iarr=find(lambda==lambdalist(il));
    ID=IDs{iarr};
    
    figure(555+il)
    M=Mp(:,:,iarr);
    pcolor(log10(E),z,M)
    ax=axis;
    axis([ax(1:2) 90 400])
    shading flat
    colorbar
    xlabel('energy [eV]')
    ylabel('altitude [km]')
    title(ID)
end %for

end %function
