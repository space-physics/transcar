%Integral() is a function to calculate the integral between 0 and x of the function exp(p^2)
function sum = Integral(x,Nu)
aux=0;
for i=1:Nu
    aux=aux+(x/Nu)*(0.5*(exp(((i-1)*x/Nu)^2)+exp((i*x/Nu)^2)));
end
sum=aux;
