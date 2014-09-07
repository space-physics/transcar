function calc_conductivity(datadir,plasmafile,confile)

inFN = [datadir,'/',plasmafile];

display(['CALC_CONDUCTIVITY: Loading densities and temperatures of ionospheric species from: ',inFN]);
load(inFN,'Enord_time','Eest_time','Bmag','n*','N*',...
          'dipangle','time','z','Tn_time','t*')

%% CONVERSION TO MKSA UNITS
%display('CALC_CONDUCTIVITY: Setting up physical constants for conductivity calculations');
EmagVm=sqrt(Enord_time.^2+Eest_time.^2)/1e3;
BmagT=Bmag/1e4;
%% CONVERSION TO CM^-3
n1cm=n1_time*1e-6;
n2cm=n2_time*1e-6;
n3cm=n3_time*1e-6;
n4cm=n4_time*1e-6;
n5cm=n5_time*1e-6;
n6cm=n6_time*1e-6;
necm=ne_time*1e-6;
nocm=No_time*1e-6;
nn2cm=Nn2_time*1e-6;
no2cm=No2_time*1e-6;
%% SETUP UP MASS AND CHARGE OF EACH SPECIES
amu=1.6605e-27;
elcharge=1.6022e-19;
m=[16, 1, 14, 28, 30, 32, 5.4858e-4]*amu;           %column dimension is ion species
q=[1, 1, 1, 1, 1, 1, -1]*elcharge;
omega=q*BmagT./m;


%3RD DIMENSION IS ION SPECIES 4TH IS NEUTRAL WE ARE COLLIDING WITH.
%COLLISIONS WITH NO ARE BEING IGNORED HERE WHICH SHOULD BE OKAY I THINK
%ALL CROSS SECTIONS ARE FROM (SCHUNK AND NAGY, 2000)
%display('CALC_CONDUCTIVITY: Computing ion-neutral collision frequencies for all 7 species');
lz=length(z);
lt=length(time);
ls=7;
nuin=zeros(lz,lt,6,3);

Tr1=(t1_time+Tn_time)/2;
nuin(:,:,1,1)=3.67e-11*nocm.*sqrt(Tr1).*(1-0.064*log10(Tr1)).^2;
nuin(:,:,1,2)=6.82e-10*nn2cm;
nuin(:,:,1,3)=6.64e-10*no2cm;

nuin(:,:,2,1)=6.61e-11*nocm.*sqrt(t2_time).*(1-0.047*log10(t2_time)).^2;
nuin(:,:,2,2)=33.6e-10*nn2cm;
nuin(:,:,2,3)=32e-10*no2cm;

nuin(:,:,3,1)=4.42e-10*nocm;
nuin(:,:,3,2)=7.47E-10*nn2cm;
nuin(:,:,3,3)=7.25e-10*no2cm;

Tr4=(tm_time+Tn_time)/2;
nuin(:,:,4,1)=2.58e-10*nocm;
nuin(:,:,4,2)=5.14e-11*nn2cm.*sqrt(Tr4).*(1-0.069*log10(Tr4)).^2;
nuin(:,:,4,3)=4.49e-10*no2cm;

nuin(:,:,5,1)=2.44e-10*nocm;
nuin(:,:,5,2)=4.34e-10*nn2cm;
nuin(:,:,5,3)=4.27e-10*no2cm;

Tr6=(tm_time+Tn_time)/2;
nuin(:,:,6,1)=2.31e-10*nocm;
nuin(:,:,6,2)=4.13e-10*nn2cm;
nuin(:,:,6,3)=2.59e-11*no2cm.*sqrt(Tr6).*(1-0.073*log10(Tr6)).^2;

nuin(:,:,7,1)=8.9e-11*nocm.*(1+5.7e-4*te_time).*sqrt(te_time);
nuin(:,:,7,2)=2.33e-11*nn2cm.*(1-1.21e-4*te_time).*te_time;
nuin(:,:,7,3)=1.82e-10*no2cm.*(1+3.6e-2*sqrt(te_time)).*sqrt(te_time);



%CONDUCTIVITY CALCULATIONS (AFTER SCHUNK AND NAGY, 2000)
display('CALC_CONDUCTIVITY: Computing Pederson and Hall conductivities');
nui=sum(nuin,4);
sigmaP=zeros(lz,lt);
sigmaH=zeros(lz,lt);
for k=1:ls
  if k<7
    nk=eval(['n',int2str(k),'_time']);
  else
    nk=ne_time;
  end
  nuk=nui(:,:,k);
  
  sigmak=nk*q(k)^2/m(k)./nuk;
  sigmaP=sigmaP+sigmak.*nuk.^2./(nuk.^2+omega(k)^2);
  sigmaH=sigmaH-sigmak.*nuk*omega(k)./(nuk.^2+omega(k)^2);
end
ibad= sigmaH<0;
sigmaH(ibad)=0;

display(['calc_conductivity: set ',int2str(sum(ibad(:))),...
          ' sigmaH values to zero.'])
%% HEIGHT INTEGRATED CONDUCTIVITIES
display('CALC_CONDUCTIVITY: Computing height integrated conductivities')
r=z/cosd(dipangle)*1e3;         %km to m conversion
rmat=repmat(r,1,lt);
[SigmaP,dr]=intrap(sigmaP,rmat);
SigmaP=SigmaP(lz-1,:);
[SigmaH,dr]=intrap(sigmaH,rmat);
SigmaH=SigmaH(lz-1,:);
%% JOULE DISSIPATION
display('CALC_CONDUCTIVITY: Computing Joule dissipation')
Kperpmag=[SigmaP.*EmagVm; SigmaH.*EmagVm];
Joulediss=Kperpmag(1,:).*EmagVm; %#ok<NASGU>
EmagVmmat=repmat(EmagVm,lz,1);
JPmag=sigmaP.*EmagVmmat;
JHmag=sigmaH.*EmagVmmat; %#ok<NASGU>
joulediss=JPmag.*EmagVmmat; %#ok<NASGU>

condFN = [datadir,'/',confile];
display(['CALC_CONDUCTIVITY: Saving conductivities and energy dissipation calculations to: ',condFN])
save(condFN,'EmagVm','BmagT','Kperpmag','joule*','Joule*','JPmag','JHmag','z',...
            'time','Sigma*','sigma*','nuin','m','q','nui','omega') 

end
