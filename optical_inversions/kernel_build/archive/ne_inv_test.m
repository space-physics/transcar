%LOAD DATA FROM CHARACTERISTIC MATRIX BUILDER
load(['kernel_',datelab,'_ne.mat']);



%CREATE DENSITY RESPONSE TO A MAXWELLIAN
le=length(ene)-1;
diffe=diff(ene);
E=ene;%(1:le)%+diffe(1:le);
E0=5000;
phiT=7.0114e12;
maxw=1/2*phiT*E/E0^3.*exp(-E/E0)/2/pi;          %cm^-2 s^-1 eV^-1 sr^-1
precipbasis=maxw(:)./phiNne(:);

nedat=Mne*precipbasis;



%CORRUPT WITH NOISE
sigman=0.1*nedat;
nedatnoisy=nedat+sigman.*randn(size(nedat));
%semilogx(nedat,zne,nedatnoisy,zne);



%ME ALGORITHM PARAMETERS
w=1/length(nedatnoisy)*ones(size(nedatnoisy));
beta=1e-7;
ic=ones(size(ene(:)))*1e2;



%PERFORM ME INVERSION
[phi,chir]=mebergne(nedatnoisy,Mne,ic,beta,w,sigman);



phiest=phi.*phiNne(:);
loglog(ene,maxw,'-',ene,phiest,'*-')