function [Tipar,Tiperp,vinperp]=perp_transport(time,z,Emag,miamu,datadate,gloc,activ)



%GEOPHYSICAL PARAMS.
glat=gloc(1); glon=gloc(2);
f107a=activ(1); f107=activ(2); ap=activ(3);
day=datadate(1); month=datadate(2); year=datadate(3);
doy=round(30.5*month+day);
yearshort=mod(year,100);
iyd=yearshort*1000+doy;



%PHYSICAL CONSTANTS
kb=1.38e-23;
amu=1.66e-27;
elchrg=1.6e-19;
B0=0.488e-4;
mi=amu*miamu;
omegai=elchrg*B0/mi;            %gyro-frequency



%MSIS NEUTRAL ATMOSPHERE THIS PARTICULAR TIME AND ALTITUDE INTERVAL
sec=round(time*3600);
neudat=msis_matlab(z,glat,glon,iyd,sec,f107a,f107,ap);
nO=neudat(:,3); nN2=neudat(:,4); nO2=neudat(:,5);
Tn=neudat(:,11);

%COLLISION FREQUENCIES
if miamu==30
    nuiO=2.44e-10*nO*1e-6;
    nuiN2=4.34e-10*nN2*1e-6;
    nuiO2=4.27e-10*nO2*1e-6;
else
    Tr1=(Tn+Tn)/2;          %assume Ti=Tn
    nuiO=3.67e-11*nO*1e-6.*sqrt(Tr1).*(1-0.064*log10(Tr1)).^2;
    nuiN2=6.82e-10*nN2*1e-6;
    nuiO2=6.64e-10*nO2*1e-6;
end

nui=nuiO+nuiN2+nuiO2;
mnnum=16*amu*nuiO/(mi+16*amu)+28*amu*nuiN2/(mi+28*amu)+32*amu*nuiO2/(mi+32*amu);
mndenom=nuiO/(mi+16*amu)+nuiN2/(mi+28*amu)+nuiO2/(mi+32*amu);
mn=mnnum./mndenom;
alphai=nui/omegai;



%CALCULATE ION TEMPERATURE PARTITIONAING COEFFS. VS ALT
if miamu==30
    for k=1:length(z)
    betapar(k)=0.49;
    betaperp(k)=0.73;
    end
else
    %these coeffs. are from linear fits to data in table 2 of [Winkler et al, 1992]
    parmperp=[0.1880; 0.6900];
    parmpar=[-0.3480; 0.5200];
    for k=1:length(z)
        pctO=nO(k)./(nN2(k)+nO(k));
        betapar(k)=parmpar(1)*pctO+parmpar(2);
        betaperp(k)=parmperp(1)*pctO+parmperp(2);
    end
end
betapar=betapar(:);
betaperp=betaperp(:);



%ION TEMPERATURE CALCULATION
Tipar=Tn+mn/3/kb*1.5.*betapar./(1+alphai.^2)*Emag^2/B0^2;
Tiperp=Tn+mn/3/kb*1.5.*betaperp./(1+alphai.^2)*Emag^2/B0^2;



%PERP VELOCITY
vinperp=sqrt(1./(1+alphai.^2)*Emag^2/B0^2);


end
