function [phi,chir]=meberg(b,M,phi0,beta,w,sigman)

%Max number of iterations
maxit=7500;

%Kernel Size
[m,n]=size(M);

%Maximum Entropy iterations
nit=1; chir=1e25; ex=0; phi=phi0;
while nit<=maxit && abs(chir-1)>1e-1 && ~ex
  %MART
  phimat=repmat(phi',m,1);
  t=min(1./(M.*phimat),[],2);
  c=beta*(1-1./b.*(M*phi)).*t;
  phi=phi./(1-phi.*(M'*diag(w)*c));
  
  %Check the convergence
  if mod(nit,100)==0

    %Let the user know what is happening
    if ~mod(nit,1000)
        fprintf('MEBERG.M: iteration %d \r',nit);
    end %if
      
    %Compute solution error.  Inefficient should be done only once per
    %function call.
    chirprev=chir;
    e=M*phi-b;
    if(norm(sigman)>0)
      chir=e'*diag(1./sigman.^2)*e; 
    else
      chir=e'*e; 
    end
      
    %Change beta if we are getting close to err. min
    dchi=chirprev-chir;
    if dchi<=0
        ex=1;
    elseif dchi<1
        beta=dchi;
    end
  end
  nit=nit+1;
end
return;

end