function C=comp_zc(z,zc)

%H=60;
H=45;
C=2./(1+(1+8*exp(-3*(z-zc)./H)).^0.5);      %after Oliver (1975)

end