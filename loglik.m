function [llik] = loglik(dt,m,maxd,data,timechange,low,up,pars)
    
    if all(pars>low) && all(pars<up) 
     [~,~,meanout,varout,~,~] = FilterDelay(dt,m,data,maxd,timechange,pars);
      varvec=reshape(varout,length(data),1);
            if isreal(varvec)==1 && sum(sum(isnan(varvec)))==0 && sum(varvec<0)==0
                llik=sum(log(normpdf(data',meanout',sqrt(varvec))));
            else
                llik=NaN;
            end    
    else
         llik=-Inf;
    end
    
end