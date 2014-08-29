function cmp=composition_sondrestrom(z,datadate,gloc,activ)


%GET THE NEUTRAL TEMPERATURE FROM MSIS
glat=gloc(1); glon=gloc(2);
f107a=activ(1); f107=activ(2); ap=activ(3);
day=datadate(1); month=datadate(2); year=datadate(3);
doy=round(30.5*month+day);
yearshort=mod(year,100);
iyd=yearshort*1000+doy;

sec=3*3600;
neudat=msis_matlab([120;1000],glat,glon,iyd,sec,f107a,f107,ap);
T_120=neudat(1,11);
T_exos=neudat(2,11);



%COMPUTE COMPOSITION AT SONDRESTROM AS PER OLIVER, 1975
H_inf=60.;
S=.01;

x=3.*(z-180.)*T_exos/H_inf./(T_exos-(T_exos-T_120)*exp(-S*(z-120.)));
cmp=2./(1.+sqrt(1.+8.*exp(-x)));
izero=find(x<-85);
cmp(izero)=0;
cmp=reshape(cmp,size(z));
% cmp=zeros(size(z));
% for k=1:length(z)
%     x  = 3.*(z(k)-180.) * T_exos / H_inf ...
%         /( T_exos - (T_exos-T_120)*exp(-S*(z(k)-120.)) );
% 
%     if x < -85.
%         cmp(k)=0.
%     else
%         cmp(k)=2. / ( 1. + sqrt( 1. + 8.*exp(-x) ) );
%     end
% end

end