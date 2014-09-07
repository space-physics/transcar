function [UT,Epmaglsq,dEpmaglsq]=estimate_Emag(istime,isz,isti,isdti,datadate,gloc,activ)



%GEOPHYSICAL PARAMETERS
glat=gloc(1); glon=gloc(2);
f107a=activ(1); f107=activ(2); ap=activ(3);
day=datadate(1); month=datadate(2); year=datadate(3);
doy=round(30.5*month+day);
yearshort=mod(year,100);                                %I don't know why this is an msis input because msis doesn't use it to calculate parameters
iyd=yearshort*1000+doy;


%IDENTIFY TIMES AND RANGES OF INTEREST
UT=istime;
%izdat=find(isz>=130 & isz<=175);        %this is based on comments by [St. Maurice et al, 1999, Section 2.1] in AG
%izdat=find(isz>=130 & isz<=157);        %this is based on comments by [St. Maurice et al, 1999, Section 2.1] and by M. McCready
izdat=find(isz>=130 & isz<=150);
lz=length(izdat);
lt=length(UT);



%WINDOW TEMPERATURE MEASUREMENTS IN ALTITUDE
Tipardat=isti(izdat,:);
dTipar=isdti(izdat,:);
zdat=isz(izdat);



%PHYSICAL CONSTANTS
kb=1.38e-23;
amu=1.66e-27;
elchrg=1.6e-19;
B0=0.488e-4;
miamu=30;
betapar=0.52;                   %partitioning coeff based on St. Maurice's comments
mi=amu*miamu;
omegai=elchrg*B0/mi;            %gyro-frequency



Epmaglsq=zeros(1,lt);
for k=1:lt
    %MSIS NEUTRAL ATMOSPHERE THIS PARTICULAR TIME AND ALTITUDE INTERVAL
    sec=round(UT(k)*3600);
    neudat=msis_matlab(zdat,glat,glon,iyd,sec,f107a,f107,ap);
    nO=neudat(:,3); nN2=neudat(:,4); nO2=neudat(:,5);
    Tn=neudat(:,11);

    %COLLISION FREQUENCIES
    nuiO=2.44e-10*nO*1e-6;
    nuiN2=4.34e-10*nN2*1e-6;
    nuiO2=4.27e-10*nO2*1e-6;
    nui=nuiO+nuiN2+nuiO2;
    mnnum=16*amu*nuiO/(mi+16*amu)+28*amu*nuiN2/(mi+28*amu)+32*amu*nuiO2/(mi+32*amu);
    mndenom=nuiO/(mi+16*amu)+nuiN2/(mi+28*amu)+nuiO2/(mi+32*amu);
    mn=mnnum./mndenom;
    alphai=nui/omegai;

    %ION-NEUTRAL TEMPERATURE DIFFERENCE
    Ti_Tn=max(0,Tipardat(:,k)-Tn);

    %FILTER OUT BAD DATA AND PERFORMS THE LEAST SQUARES FIT
    ifs=find(Ti_Tn>0);
    if length(ifs) == 0
        Epmaglsq(k)=0;
        dEpmaglsq(k)=0;
    else
        y=Ti_Tn(ifs);                                                   %data
        dy=dTipar(ifs,k);                                               %uncertainty in temperature estimates
        A=mn(ifs)/3/kb*1.5*betapar./(1+alphai(ifs).^2)/B0^2;            %kernel
        W=diag(1./dy.^2);                                               %compute weights for fit
%        W=eye(length(ifs));
        Aginv=inv(A'*W*A)*A'*W;                                         %weighted least squares inverse of A
        x=Aginv*y;                                                      %weighted least squares solution for E^2
        dx2=Aginv*diag(dy.^2)*Aginv';                                   %Covariance matrix of x given uncorrelated y
        dx=sqrt(dx2);
        Epmaglsq(k)=sqrt(x);                                            %E estimate
        dEpmaglsq(k)=Epmaglsq(k)*0.5*dx/x;                              %Propogate uncertainty through square root operation (see Wikipedia article on Propagation of uncertainty)
    end
    if isnan(Epmaglsq(k)) | isnan(dEpmaglsq(k))
        Epmaglsq(k)=0;
        dEpmaglsq(k)=0;
    end
end

end