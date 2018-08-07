function [lprior] = lpriorGa(par,a,b)
for kkk=1:length(par)
    if (par(kkk)>0)
    lprior(kkk)=log(gampdf(par(kkk),a(kkk),b(kkk)));
    else
    lprior(kkk)=-Inf;
    end
end
end

