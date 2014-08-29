function [Up,Sp,Vp,Mp,Mpinv]=tsvd(M,tlev)
  [U,S,V]=svd(M);
  
  Up=U(:,1:tlev);
  Sp=S(1:tlev,1:tlev);
  Vp=V';
  Vp=Vp(1:tlev,:)';
  Mp=Up*Sp*Vp';
  Mpinv=Vp*inv(Sp)*Up';
end