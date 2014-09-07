%This file generate a ISR spectrum based on the Sheffield approach. 

% Outputs
% Spec: Spectrum vector 
% w: Radial vector in rad*Hz
% 
% Inputs
% f: Radar frequency in MHz
% Te: Electron temperature in Kelvin
% Ti: Ion temperature in Kelvin
% n: Electron (ion) density in 1/m^3
% fmi: Factor of the ion mass - composition (fmi=1 for Hydrogen and fmi=16 for Oxigen) 
% Z: Ionization number (e.g. H+ or O+ has Z=1, H++ or O++ has Z=2)
% dt: Time sampling step in sec
% Ni: Length desired of the output vectors

function [Spec w]=isr(f,Te,Ti,n,fmi,Z,dt,Ni)
Nu=10000;
me=9.1094e-31;
mi=1.6726e-27*fmi;
kb=1.3807e-23;
q=1.6022e-19;
eps0=8.8542e-12;
c=3e8;
lamb=c/(f*1e6); %Radar wavelength 
k=4*pi/lamb; %Received wavenumber (2*kr, kr: radar wavenumber)

lambD=(eps0*kb*Te)/(n*q^2); %Debye length
a=sqrt(2*kb*Te/me); %Electron thermal speed 
b=sqrt(2*kb*Ti/mi); %Ion thermal speed
alpha=1/(k*lambD);

for j=1:Ni
    omega=2*pi*(j-ceil((Ni+1)/2))/(dt*Ni); %Radial frequency
    xe=omega/(k*a);
    xi=omega/(k*b);
    %Integral() is a function to calculate the integral between 0 and x of
    %the function exp(p^2)
    Gexi=(alpha^2)*(1-2*xi*exp(-(xi)^2)*Integral(xi,Nu)-i*sqrt(pi)*xi*exp(-(xi)^2));
    Gexe=(alpha^2)*(1-2*xe*exp(-(xe)^2)*Integral(xe,Nu)-i*sqrt(pi)*xe*exp(-(xe)^2));
    Gixe=(Z*Te/Ti)*Gexe;
    Gixi=(Z*Te/Ti)*Gexi;
    eps=1+Gexe+Gixi; %Plasma permitivity
    
    fe0=(1/sqrt(pi*a^2))*exp(-1*((omega/k)^2)/(a^2)); %Maxwellian distribution for electrons
    fi0=(1/sqrt(pi*b^2))*exp(-1*((omega/k)^2)/(b^2)); %Maxwellain distribution for ions

    aux=(2*pi/k)*(abs(1-(Gexe/eps))^2)*fe0+(2*pi*Z/k)*(abs(Gexe/eps)^2)*fi0;
    %This "if" is used to change the "NaN" for "zeros" when the spectrum goes to small  
    if isnan(aux)==1
        Spec(j)=0.0;
    else
        Spec(j)=aux;
    end
    w(j)=omega; %Radial frequency vector
end

