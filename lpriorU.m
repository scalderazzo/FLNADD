function [lprior] = lpriorU(par,a,b)
if (par<a) || (par>b) 
lprior=-Inf;
else
lprior=log(1/(b-a));
end

