function [lprior] = lpriorN(par,a,b)
lprior=log(normpdf(par,a,b));
end
