addpath ~/matlabshared;

E0=[100,300,500,1000,3000,10000,30000];
PhiE_mWm2=[1,5.29,5.29,10,10,10,10];
PhiE=PhiE_mWm2*1e-3/1.6e-19*1e-4;        %conversion to eV/cm2/s
t=[11850,14400]/3600;

for k=1:length(E0)
    %MAXWELLIAN DIST
    E=logspace(1,4.5,25); E=E(:);
    phiN(:,k)=1/4/pi/2/E0(k)^3*PhiE(k)*E.*exp(-E/E0(k));
    filename=['precinput_',num2str(E0(k))];
    create_precinput(t,E,[phiN(:,k),phiN(:,k)],filename);
end

loglog(E,phiN)
axis([30,3e4,1e2,1e8])

rmpath ~/matlabshared;