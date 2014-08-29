%This file generate a ISR spectrum quasi-based on the Sheffield approach.
%           Created by M. Diaz
%           Edited by  D. Ismail       (drift velocity added)
%           Edited by  M. Zettergren   (multiple species added)
%
%[Spec w]=IS_spectrum(f,Te,Ti,ni,Vi,miamu,Z,dt,nw)
%
% Output --> variable name : description
% Spec: Spectrum vector 
% w: Radial vector in rad*Hz
% 
% Input --> variable name : description (type, eg. scalar or vector)
% f: Radar frequency in MHz
% Te: Electron temperature in Kelvin        (scalar)
% Ti: Ion temperatures in Kelvin            (vector, one element for each ion species)
% ni: ion densities in 1/m^3                (vector, one element for each ion species)
% Vi: Plasma drify velocity in m/s          (scalar, all species are assumed to drift at same speed)
% miamu: ion masses in atomic mass units    (vector, one element for each ion species) 
% Z: Ionization number                      (vector, one element for each ion species)
% dt: Time sampling step in sec.  When mi=16 a good dt is 5e-5.
% nw: Length desired of the output vectors.  When fmi=16 a good nw is 1000.


function [Spec w]=IS_spectrum(f,Te,Ti,ni,Vi,miamu,Z,dt,nw)
Nu=100;
me=9.1094e-31;
amu=1.6726e-27;
kb=1.3807e-23;
q=1.6022e-19;
eps0=8.8542e-12;
c=3e8;
lamb=c/(f*1e6);             %Radar wavelength 
k=4*pi/lamb;                %Received wavenumber (2*kr, kr: radar wavenumber)

ne=sum(ni);
Dfact=sqrt(eps0*kb*Te/q^2);
%lambD=sqrt(eps0*kb*Te/ne/q^2);   %Debye length
%alpha=1/(k*lambD);


a=sqrt(2*kb*Te/me);                 %Electron thermal speed 
mi=miamu*amu;                       %Ion masses in kg
ns=length(mi);                      %# of ion species
for is=1:ns                         %Thermal speeds for each species
  b(is)=sqrt(2*kb*Ti(is)/mi(is));                 
end

for iw=1:nw
    omega=2*pi*(iw-ceil((nw+1)/2))/(dt*nw);    %Radial frequency

    xe=(omega/k-Vi)/a;                      %Electron terms in the spectrum
    alpha=1/k/(Dfact/sqrt(ne));
    Gexe=(alpha^2)*(1-2*xe*exp(-(xe)^2)*intexp2(xe,Nu)-i*sqrt(pi)*xe*exp(-(xe)^2));
    fe0=ne*(1/sqrt(pi*a^2))*exp(-1*((omega/k-Vi)^2)/(a^2));                                %Maxwellian distribution for electrons
    
    Gixi=0;
    for is=1:ns
        xi=(omega/k-Vi)/b(is);           %Ion terms in the spectrum
        alpha=1/k/(Dfact/sqrt(ni(is)));
        Gixi=Gixi+Z(is)*Te/Ti(is)*alpha^2* ...
            (1-2*xi*exp(-(xi)^2)*intexp2(xi,Nu)- ...
            i*sqrt(pi)*xi*exp(-(xi)^2));  
        fi0(is)=ni(is)*(1/sqrt(pi*b(is)^2))*exp(-1*((omega/k-Vi)^2)/(b(is)^2));            %Maxwellain distribution for ions
    end
    eps=1+Gexe+Gixi;                     %Plasma permitivity

    aux=(2*pi/k)*(abs(1-(Gexe/eps))^2)*fe0+2*pi/k*abs(Gexe/eps)^2*sum(Z.*fi0);
    %This "if" is used to change the "NaN" for "zeros" when the spectrum
    %goes to small numbers
    if isnan(aux)==1
        Spec(iw)=0.0;
    else
        Spec(iw)=aux;
    end
    w(iw)=omega; %Radial frequency vector
end