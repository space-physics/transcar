function C=inv_waldteufel_Fi(Fi)

sFi=size(Fi);

ibad=find(Fi>1);
Fi(ibad)=1;
minwald=waldteufel_Fi(1);
ibad=find(Fi<minwald);
Fi(ibad)=minwald;

a=0.785*ones(sFi);
b=-1*Fi-1.2535;
c=2.1058-2.1*Fi;

for k=1:max(sFi)
    soln(1)=(-b(k)+sqrt(b(k)^2-4*a(k)*c(k)))/2/a(k);
    soln(2)=(-b(k)-sqrt(b(k)^2-4*a(k)*c(k)))/2/a(k);

    if isreal(soln(1)) & soln(1)<=1 & soln(1)>=0
        C(k)=soln(1);
    else
        C(k)=soln(2);
    end
end
C=reshape(C,sFi);

end