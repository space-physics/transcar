%LOAD DATA FROM CHARACTERISTIC MATRIX BUILDER
load(['kernel_',datelab,'.mat']);



%PICK OUT ALTITUDES/ENERGIES WHERE CHEMICAL EQUILIBRIUM IS GOOD ASSUMPTION
iz1=3;
iz2=min(find(z>=215));
ie1=min(find(energies>=1e3));
ie2=length(energies)-1;



%FORM THE KERNEL
Mne=nechr(iz1:iz2,ie1:ie2);
zne=z(iz1:iz2);
ene=energies(ie1:ie2);
phiNne=phiN(ie1:ie2);


%SAVE KERNEL TO A SEPARATE FILE
save(['kernel_',date,'_ne.mat'],'Mne','energies','zne','ene','phiNne','datelab');