function [xout,errorout,outlogL] = MCMC_Cry_DA(it,deltat,cry,th0,a,b,up,low,maxd,timechange,m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%initialisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%vector for errors
error=zeros(1,29);
%output samples
out=th0;

%evaluate priors in log parametrisation
lpro(1)=lpriorN(-th0(1),a(1),b(1));
lpro(2)=lpriorN((th0(2)+th0(1))/exp(th0(3)),a(2),b(2));
lpro(3)=lpriorN(th0(3),a(3),b(3));
lpro(4)=lpriorN(th0(4),a(4),b(4));
lpro(5)=lpriorU(th0(5),a(5),b(5));
lpro(6)=lpriorU(th0(5)/th0(6),a(6),b(6));
lpro(7:14)=lpriorN(th0(7:14),a(7:14),b(7:14));

%initialise algortihm samples
xol=th0;

%transform back original parametrisation
rep=[1/exp(th0(1)),exp((th0(2)+th0(1))/exp(th0(3))),exp(th0(3)),exp(th0(4)),th0(5),th0(5)/th0(6),exp(th0(7:14))];

%evaluate loglikelihood in original parametrisation
logLo=loglik(deltat,m,maxd,cry,timechange,low,up,rep);

%sample initial values until loglikelihood is real, not missing and finite
while isnan(logLo)==1 || logLo==-Inf || isreal(logLo)==0
    th0=zeros(1,14);
    th0([1:4,7:14])=normrnd(a([1:4,7:14]),b([1:4,7:14]));
    th0(5:6)=unifrnd(a(5:6),b(5:6));
    lpro(1)=lpriorN(-th0(1),a(1),b(1));
    lpro(2)=lpriorN((th0(2)+th0(1))/exp(th0(3)),a(2),b(2));
    lpro(3)=lpriorN(th0(3),a(3),b(3));
    lpro(4)=lpriorN(th0(4),a(4),b(4));
    lpro(5)=lpriorU(th0(5),a(5),b(5));
    lpro(6)=lpriorU(th0(5)/th0(6),a(6),b(6));
    lpro(7:14)=lpriorN(th0(7:14),a(7:14),b(7:14));
    
    xol=th0;
    out=th0;
    rep=[1/exp(th0(1)),exp((th0(2)+th0(1))/exp(th0(3))),exp(th0(3)),exp(th0(4)),th0(5),th0(5)/th0(6),exp(th0(7:14))];
    logLo=loglik(deltat,m,maxd,cry,timechange,low,up,rep);
end

%initialise stored loglikelihood values
logLik=zeros(1,it);
logLik(1)=logLo;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%pilot run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first 10^3 iterations to obtain an estimate of the
%variance/covariance matrix

for i=1:1000
    
    for j=1:6
        %initialise
        xs=xol;
        lprs=lpro;
        %propose
        xs(j)=xol(j)+normrnd(0,1);
        %evaluate priors
        lprs(1)=lpriorN(-xs(1),a(1),b(1));
        lprs(2)=lpriorN(((xs(2)+xs(1))/exp(xs(3))),a(2),b(2));
        lprs(3:4)=lpriorN(xs(3:4),a(3:4),b(3:4));
        lprs(5)=lpriorU(xs(5),a(5),b(5));
        lprs(6)=lpriorU(xs(5)/xs(6),a(6),b(6));
        %reparametrise
        rep=[1/exp(xs(1)),exp((xs(2)+xs(1))/exp(xs(3))),exp(xs(3)),exp(xs(4)),xs(5),xs(5)/xs(6),exp(xs(7:14))];
        %evaluate loglikelihood
        logLs=loglik(deltat,m,maxd,cry,timechange,low,up,rep);
        %set loglikelihood to -Inf if NA and store errors
        if isnan(logLs)==1
            error=[error;i,xol,xs];
            logLs=-Inf;
        end
        %acceptance ratio (Metropolis-Hastings) - added Jacobian of
        %transformation
        al=min(1,exp(sum(lprs)+logLs-xs(3)+log(xs(5)/xs(6)^2)-(sum(lpro)+logLo-xol(3)+log(xol(5)/xol(6)^2))));
        %accept/reject and update
        u=rand(1);
        if (u<al)
            xol=xs;
            lpro=lprs;
            logLo=logLs;
        end
    end
    

    xs=xol;
    lprs=lpro;
    sig2=(eye(6))/6;
    xs(7:12)=mvnrnd(xol(7:12),sig2);
    lprs(7:12)=lpriorN(xs(7:12),a(7:12),b(7:12));
    rep=[1/exp(xs(1)),exp((xs(2)+xs(1))/exp(xs(3))),exp(xs(3)),exp(xs(4)),xs(5),xs(5)/xs(6),exp(xs(7:14))];
    logLs=loglik(deltat,m,maxd,cry,timechange,low,up,rep);
    if isnan(logLs)==1
        error=[error;i,xol,xs];
        logLs=-Inf;
    end
    al=min(1,exp(sum(lprs)+logLs-xs(3)+log(xs(5)/xs(6)^2)-(sum(lpro)+logLo-xol(3)+log(xol(5)/xol(6)^2))));
    u=rand(1);
    if (u<al)
        xol=xs;
        lpro=lprs;
        logLo=logLs;
    end
    
    xs=xol;
    lprs=lpro;
    sig3=(eye(2))/2;
    xs(13:14)=mvnrnd(xol(13:14),sig3);
    lprs(13:14)=lpriorN(xs(13:14),a(13:14),b(13:14));
    rep=[1/exp(xs(1)),exp((xs(2)+xs(1))/exp(xs(3))),exp(xs(3)),exp(xs(4)),xs(5),xs(5)/xs(6),exp(xs(7:14))];
    logLs=loglik(deltat,m,maxd,cry,timechange,low,up,rep);
    if isnan(logLs)==1
        error=[error;i,xol,xs];
        logLs=-Inf;
    end
    al=min(1,exp(sum(lprs)+logLs-xs(3)+log(xs(5)/xs(6)^2)-(sum(lpro)+logLo-xol(3)+log(xol(5)/xol(6)^2))));
    u=rand(1);
    if (u<al)
        xol=xs;
        lpro=lprs;
        logLo=logLs;
    end
    
    out(i+1,:)=xol;
    logLik(i+1)=logLo;
    
end

%compute mean and variance covariance matrices for each block

muN1=mean(out(1:i,1:6));
sigC1=cov(out(1:i,1:6));
muN2=mean(out(1:i,7:12));
sigC2=cov(out(1:i,7:12));
muN3=mean(out(1:i,13:14));
sigC3=cov(out(1:i,13:14));

%sample according to variance/covariance matrix of accepted values

for i=1001:30000
    
    %initialise
   xs=xol;
    lprs=lpro;
    %update varance/covariance matrix
    muO1=muN1;
    muN1=muO1+(out(i,1:6)-muO1)./i;
    sigC1=sigC1.*((i-2)/(i-1))+(1/(i-1)).*((out(i,1:6)-muN1)'*(out(i,1:6)-muO1));
    %sample from mixture (Roberts and Rosenthal, 2009)
    if (u<0.95)
        sigma1=(sigC1./6);
        xs(1:6)=mvnrnd(xol(1:6),sigma1);
    else
        sigma2=(((0.1)^2.*eye(6))./6);
        xs(1:6)=mvnrnd(xol(1:6),sigma2);
    end
    lprs(1)=lpriorN(-xs(1),a(1),b(1));
    lprs(2)=lpriorN(((xs(2)+xs(1))/xs(3)),a(2),b(2));
    lprs(3:4)=lpriorN(xs(3:4),a(3:4),b(3:4));
    lprs(5)=lpriorU(xs(5),a(5),b(5));
    lprs(6)=lpriorU(xs(5)/xs(6),a(6),b(6));
    rep=[1/exp(xs(1)),exp((xs(2)+xs(1))/exp(xs(3))),exp(xs(3)),exp(xs(4)),xs(5),xs(5)/xs(6),exp(xs(7:14))];
    logLs=loglik(deltat,m,maxd,cry,timechange,low,up,rep);
    if isnan(logLs)==1
        error=[error;i,xol,xs];
        logLs=-Inf;
    end
    al=min(1,exp(sum(lprs)+logLs-xs(3)+log(xs(5)/xs(6)^2)-(sum(lpro)+logLo-xol(3)+log(xol(5)/xol(6)^2))));
    u=rand(1);
    if (u<al)
        xol=xs;
        lpro=lprs;
        logLo=logLs;
    end
    

     xs=xol;
    lprs=lpro;
    muO2=muN2;
    muN2=muO2+(out(i,7:12)-muO2)./i;
    sigC2=sigC2.*((i-2)/(i-1))+(1/(i-1)).*((out(i,7:12)-muN2)'*(out(i,7:12)-muO2));
    if (u<0.95)
        sigma3=((2.38^2).*sigC2./6);
        xs(7:12)=mvnrnd(xol(7:12),sigma3);
    else
        sigma4=(((0.1)^2.*eye(6))./6);
        xs(7:12)=mvnrnd(xol(7:12),sigma4);
    end
    lprs(7:12)=lpriorN(xs(7:12),a(7:12),b(7:12));
    rep=[1/exp(xs(1)),exp((xs(2)+xs(1))/exp(xs(3))),exp(xs(3)),exp(xs(4)),xs(5),xs(5)/xs(6),exp(xs(7:14))];
    logLs=loglik(deltat,m,maxd,cry,timechange,low,up,rep);
    if isnan(logLs)==1
        error=[error;i,xol,xs];
        logLs=-Inf;
    end
    al=min(1,exp(sum(lprs)+logLs-xs(3)+log(xs(5)/xs(6)^2)-(sum(lpro)+logLo-xol(3)+log(xol(5)/xol(6)^2))));
    u=rand(1);
    if (u<al)
        xol=xs;
        lpro=lprs;
        logLo=logLs;
    end
    
    
   xs=xol;
    lprs=lpro;
    muO3=muN3;
    muN3=muO3+(out(i,13:14)-muO3)./i;
    sigC3=sigC3.*((i-2)/(i-1))+(1/(i-1)).*((out(i,13:14)-muN3)'*(out(i,13:14)-muO3));
    if (u<0.95)
        sigma5=((2.38^2).*sigC3./2);
        xs(13:14)=mvnrnd(xol(13:14),sigma5);
    else
        sigma6=(((0.1)^2.*eye(2))./2);
        xs(13:14)=mvnrnd(xol(13:14),sigma6);
    end
    lprs(13:14)=lpriorN(xs(13:14),a(13:14),b(13:14));
    rep=[1/exp(xs(1)),exp((xs(2)+xs(1))/exp(xs(3))),exp(xs(3)),exp(xs(4)),xs(5),xs(5)/xs(6),exp(xs(7:14))];
    logLs=loglik(deltat,m,maxd,cry,timechange,low,up,rep);
    if isnan(logLs)==1
        error=[error;i,xol,xs];
        logLs=-Inf;
    end
    al=min(1,exp(sum(lprs)+logLs-xs(3)+log(xs(5)/xs(6)^2)-(sum(lpro)+logLo-xol(3)+log(xol(5)/xol(6)^2))));
    u=rand(1);
    if (u<al)
        xol=xs;
        lpro=lprs;
        logLo=logLs;
    end
    
    out(i+1,:)=xol;
    logLik(i+1)=logLo;
    pause(0.01)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Delayed Acceptance algorithm run (Christen and Fox, 2005; Golightly et
%%al., 2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

muN1=mean(out(1:i,1:6));
sigC1=cov(out(1:i,1:6));
muN2=mean(out(1:i,7:12));
sigC2=cov(out(1:i,7:12));
muN3=mean(out(1:i,13:14));
sigC3=cov(out(1:i,13:14));

%likelihood evaluated for dt=0.5h
logLoA=logLo;
%likelihood evaluated for dt=0.1h
xs=xol;
rep01=[(1/exp(xs(1)))/5,(exp((xs(2)+xs(1))/exp(xs(3))))/5,exp(xs(3)),exp(xs(4)),xs(5),xs(5)/xs(6),(exp(xs(7:9)))/5,exp(xs(10)),(exp(xs(11:13)))/5,exp(xs(14))];
logLoF=loglik(deltat/5,m,maxd,cry,timechange,low,up,rep01);

for i=30001:it
    
    xs=xol;
    lprs=lpro;
    muO1=muN1;
    muN1=muO1+(out(i,1:6)-muO1)./i;
    sigC1=sigC1.*((i-2)/(i-1))+(1/(i-1)).*((out(i,1:6)-muN1)'*(out(i,1:6)-muO1));
    if (u<0.95)
        sigma1=(sigC1./6);
        xs(1:6)=mvnrnd(xol(1:6),sigma1);
    else
        sigma2=(((0.1)^2.*eye(6))./6);
        xs(1:6)=mvnrnd(xol(1:6),sigma2);
    end
    
    lprs(1)=lpriorN(-xs(1),a(1),b(1));
    lprs(2)=lpriorN(((xs(2)+xs(1))/exp(xs(3))),a(2),b(2));
    lprs(3:4)=lpriorN(xs(3:4),a(3:4),b(3:4));
    lprs(5)=lpriorU(xs(5),a(5),b(5));
    lprs(6)=lpriorU(xs(5)/xs(6),a(6),b(6));
    rep=[1/exp(xs(1)),exp((xs(2)+xs(1))/exp(xs(3))),exp(xs(3)),exp(xs(4)),xs(5),xs(5)/xs(6),exp(xs(7:14))];
    logLsA=loglik(deltat,m,maxd,cry,timechange,low,up,rep);
    if isnan(logLsA)==1
        error=[error;i,xol,xs];
        logLsA=-Inf;
    end
    
    ratioA=exp(sum(lprs)+logLsA-xs(3)+log(xs(5)/xs(6)^2)-(sum(lpro)+logLoA-xol(3)+log(xol(5)/xol(6)^2)));
    al=min(1,ratioA);
    
    u=rand(1);
    if (u<al)
        rep01=[(1/exp(xs(1)))/5,(exp((xs(2)+xs(1))/exp(xs(3))))/5,exp(xs(3)),exp(xs(4)),xs(5),xs(5)/xs(6),(exp(xs(7:9)))/5,exp(xs(10)),(exp(xs(11:13)))/5,exp(xs(14))];
        logLsF=loglik(deltat/5,m,maxd,cry,timechange,low,up,rep01);
        if isnan(logLsF)==1
            error=[error;i,xol,xs];
            logLsF=-Inf;
        end
        al=min(1,(exp(sum(lprs)+logLsF-xs(3)+log(xs(5)/xs(6)^2)-(sum(lpro)+logLoF-xol(3)+log(xol(5)/xol(6)^2))))./ratioA);
        u=rand(1);
        if (u<al)
            xol=xs;
            lpro=lprs;
            logLoA=logLsA;
            logLoF=logLsF;
        end
    end
    %%here
    
    xs=xol;
    lprs=lpro;
    muO2=muN2;
    muN2=muO2+(out(i,7:12)-muO2)./i;
    sigC2=sigC2.*((i-2)/(i-1))+(1/(i-1)).*((out(i,7:12)-muN2)'*(out(i,7:12)-muO2));
    if (u<0.95)
        sigma3=((2.38^2).*sigC2./6);
        xs(7:12)=mvnrnd(xol(7:12),sigma3);
    else
        sigma4=(((0.1)^2.*eye(6))./6);
        xs(7:12)=mvnrnd(xol(7:12),sigma4);
    end
    lprs(7:12)=lpriorN(xs(7:12),a(7:12),b(7:12));
    rep=[1/exp(xs(1)),exp((xs(2)+xs(1))/exp(xs(3))),exp(xs(3)),exp(xs(4)),xs(5),xs(5)/xs(6),exp(xs(7:14))];
    logLsA=loglik(deltat,m,maxd,cry,timechange,low,up,rep);
    if isnan(logLsA)==1
        error=[error;i,xol,xs];
        logLsA=-Inf;
    end
    
    ratioA=exp(sum(lprs)+logLsA-xs(3)+log(xs(5)/xs(6)^2)-(sum(lpro)+logLoA-xol(3)+log(xol(5)/xol(6)^2)));
    al=min(1,ratioA);
    u=rand(1);
    if (u<al)
        rep01=[(1/exp(xs(1)))/5,(exp((xs(2)+xs(1))/exp(xs(3))))/5,exp(xs(3)),exp(xs(4)),xs(5),xs(5)/xs(6),(exp(xs(7:9)))/5,exp(xs(10)),(exp(xs(11:13)))/5,exp(xs(14))];
        logLsF=loglik(deltat/5,m,maxd,cry,timechange,low,up,rep01);
        if isnan(logLsF)==1
            error=[error;i,xol,xs];
            logLsF=-Inf;
        end
        al=min(1,(exp(sum(lprs)+logLsF-xs(3)+log(xs(5)/xs(6)^2)-(sum(lpro)+logLoF-xol(3)+log(xol(5)/xol(6)^2))))./ratioA);
        u=rand(1);
        if (u<al)
            xol=xs;
            lpro=lprs;
            logLoA=logLsA;
            logLoF=logLsF;
        end
    end
    
    
    xs=xol;
    lprs=lpro;
    muO3=muN3;
    muN3=muO3+(out(i,13:14)-muO3)./i;
    sigC3=sigC3.*((i-2)/(i-1))+(1/(i-1)).*((out(i,13:14)-muN3)'*(out(i,13:14)-muO3));
    if (u<0.95)
        sigma5=((2.38^2).*sigC3./2);
        xs(13:14)=mvnrnd(xol(13:14),sigma5);
    else
        sigma6=(((0.1)^2.*eye(2))./2);
        xs(13:14)=mvnrnd(xol(13:14),sigma6);
    end
    lprs(13:14)=lpriorN(xs(13:14),a(13:14),b(13:14));
    rep=[1/exp(xs(1)),exp((xs(2)+xs(1))/exp(xs(3))),exp(xs(3)),exp(xs(4)),xs(5),xs(5)/xs(6),exp(xs(7:14))];
    logLsA=loglik(deltat,m,maxd,cry,timechange,low,up,rep);
    if isnan(logLsA)==1
        error=[error;i,xol,xs];
        logLsA=-Inf;
    end
    
    ratioA=exp(sum(lprs)+logLsA-xs(3)+log(xs(5)/xs(6)^2)-(sum(lpro)+logLoA-xol(3)+log(xol(5)/xol(6)^2)));
    al=min(1,ratioA);
    u=rand(1);
    if (u<al)
        rep01=[(1/exp(xs(1)))/5,(exp((xs(2)+xs(1))/exp(xs(3))))/5,exp(xs(3)),exp(xs(4)),xs(5),xs(5)/xs(6),(exp(xs(7:9)))/5,exp(xs(10)),(exp(xs(11:13)))/5,exp(xs(14))];
        logLsF=loglik(deltat/5,m,maxd,cry,timechange,low,up,rep01);
        if isnan(logLsF)==1
            error=[error;i,xol,xs];
            logLsF=-Inf;
        end
        al=min(1,(exp(sum(lprs)+logLsF-xs(3)+log(xs(5)/xs(6)^2)-(sum(lpro)+logLoF-xol(3)+log(xol(5)/xol(6)^2))))./ratioA);
        u=rand(1);
        if (u<al)
            xol=xs;
            lpro=lprs;
            logLoA=logLsA;
            logLoF=logLsF;
        end
    end
    
    out(i+1,:)=xol;
    logLik(i+1)=logLoF;
    
    if (rem(i,1000)==0)
        save('CryDataMCMC');
    end
    
    pause(0.01)
    
end

xout=out;
errorout=error;
outlogL=logLik;
end
