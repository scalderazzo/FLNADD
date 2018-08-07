function [ths,thp,meanout,varout,Pxxs,Pxxp] = FilterDelay(dt,m,data,maxd,timechange,pars)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%initialisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%data
y=data;
n=length(y);

%parameters : transcription function
Rmax=pars(1);   %maximum transcription rate
Kpc=pars(2);    %thhreshold
Hn=pars(3);     %Hill coefficient
mu=pars(4);     %degradation rate
meanD=pars(5);  %mean delay
sdD=pars(6);    %standard deviation delay

%parameters : initial condition
ni1=pars(7);    %transcription rate before switch
ni2=pars(8);    %transcription rate after switch
beta=pars(9);   %scale coefficient
mu0=pars(10);   %degradation rate
rho0=pars(11);  %initial mean
P0=pars(12);    %initial variance

%parameters : measurement
kappa=pars(13);  %scale
sigma=pars(14);  %measurement error standard deviation

%algorithm parameters
%m                     %number of observations to predict before update
samp=0.5;              %sampling time-interval
integrint=samp/dt;     %number of unobserved states to integrate over (measurement equation)
delaydt=maxd/dt;       %number of unobserved states to integrate over (delay)

%delay evaluation grid:
timegrid=0:dt:(maxd-dt);
timegridH=(timegrid+timegrid+dt)./2;
evalp=maxd-timegridH;

%delay weights
aP=meanD^2/sdD.^2;
bP=sdD.^2/meanD;
pdP = makedist('Gamma','a',aP,'b',bP);
gammatrP = truncate(pdP,0,maxd);        %truncated gamma
weightsM=pdf(gammatrP,evalp);           %delay weights (unnormalised)
weightsMN=weightsM./(sum(weightsM));    %delay weights (normalised)


%initialisation
meanout=zeros(m,(n-maxd/samp)/m);                   %predicted mean  (observed)
varout=zeros(m,m,(n-maxd/samp)/m);                  %predicted variance (observed)

PxxI=zeros(integrint*n,integrint*n);      %current variance/covariance matrix (unobserved)
Pxxs=zeros(integrint*n,integrint*n);      %updated variance/covariance matrix (unobserved)
Pxxp=zeros(1,integrint*n);                     %predicted variance (unobserved)
ths=zeros(1,integrint*n);                      %updated mean  (unobserved)
thp=zeros(1,integrint*n);                      %predicted mean  (unobserved)


%auxiliary matrices and vectors for vectorised computations
A=1-mu*dt;
A1=1-mu0*dt;

FF=zeros(maxd/samp,maxd/samp*integrint);
vecf=zeros(1,integrint);
vecf(1:end)=ones(1,integrint);
ww=1;
for jj=1:maxd/samp
    FF(jj,ww:ww+integrint-1)=vecf;
    ww=ww+integrint;
end


%set counters of unobserved states
vv=0;
M=0;

%set counter of observations
ycount=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%filtering of initial condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%inital unobserved mean (until frst observation)
thI=zeros(1,integrint);
thI(1)=rho0;
if integrint>1
    for jj=2:(integrint)
        thI(jj)=A1*thI(jj-1)+dt*ni1;
    end
end

%inital unobserved variance/covariance (until frst observation)
PxxI(1,1)=P0;
if integrint>1
    gg=2;
    for jj=2:integrint
        PxxI(gg,gg)=A1*PxxI(gg-1,gg-1)*A1'+dt*beta*(mu0*thI(jj-1)+ni1);
        PxxI(1:jj-1,gg)=PxxI(1:jj-1,gg-1)*A1';
        PxxI(gg,1:jj-1)=PxxI(1:jj-1,gg)';
        gg=gg+1;
    end
end

%update unobserved states counters
vv=vv+integrint*m;
M=M+integrint*m;

%initial covariance between x and y
Pxy=PxxI(vv-integrint+1:vv,vv-integrint+1:vv)*(FF(end-m+1:end,end-integrint+1:end))';

%initial variance of y (no measurement error)
Pyy=(FF(end-m+1:end,end-integrint+1:end))*PxxI(vv-integrint+1:vv,vv-integrint+1:vv)*(FF(end-m+1:end,end-integrint+1:end))';

%predicted mean and variance ofunobserved states
thp(vv-(m*integrint)+1:vv)=thI(vv-(m*integrint)+1:vv);
Pxxp(vv-(m*integrint)+1:vv)=diag(PxxI(vv-(m*integrint)+1:vv,vv-(m*integrint)+1:vv));

%prediciton of current data mean and variance
meanout(:,1)=FF(1:m,1:m*integrint)*thI(vv-(m*integrint)+1:vv)';
varout(:,:,1)=Pyy+sigma.^2.*eye(m);

%updated data counter
ycount=ycount+m;

%%update of unobserved states
ths(vv-integrint+1:vv)=thI(vv-integrint+1:vv)'+Pxy/(Pyy+sigma.^2.*eye(m))*(y(ycount+1-m:ycount)'-FF(1:m,1:m*integrint)*thI(vv-integrint+1:vv)');
Pxxs(vv-integrint+1:vv,vv-integrint+1:vv)=PxxI(vv-integrint+1:vv,vv-integrint+1:vv)-Pxy/(Pyy+sigma.^2.*eye(m))*Pxy';

%updated current vector/matrices
thI=ths;
PxxI(vv-integrint+1:vv,vv-integrint+1:vv)=Pxxs(vv-integrint+1:vv,vv-integrint+1:vv);

%repeat until first switch time

for s=2:((maxd/samp)/m)
    
    if s<=(timechange/samp)/m
        ni=ni1;
    else
        ni=ni2;
    end
    
    for jj=M:(M+(integrint*m)-1)
        thI(jj+1)=A1*thI(jj)+dt*ni;
        PxxI(jj+1,jj+1)=A1*PxxI(jj,jj)*A1'+dt*beta*(mu0*thI(jj)+ni);
        PxxI(1:jj+1,jj+2)=PxxI(1:jj+1,jj+1)*A1';
        PxxI(jj+2,1:jj+1)=PxxI(1:jj+1,jj+2)';
        vv=vv+1;
    end
    
    M=M+integrint*m;
    
    Pxy=PxxI(vv-integrint+1:vv,vv-integrint+1:vv)*(FF(end-m+1:end,end-integrint+1:end))';
    Pyy=(FF(end-m+1:end,end-integrint+1:end))*PxxI(vv-integrint+1:vv,vv-integrint+1:vv)*(FF(end-m+1:end,end-integrint+1:end))';
    
    thp(vv-(m*integrint)+1:vv)=thI(vv-(m*integrint)+1:vv)';
    Pxxp(vv-(m*integrint)+1:vv)=diag(PxxI(vv-(m*integrint)+1:vv,vv-(m*integrint)+1:vv));
    
    meanout(:,s)=FF(1:m,1:m*integrint)*thI(vv-(m*integrint)+1:vv)';
    varout(:,:,s)=Pyy+sigma.^2.*eye(m);
    
    ycount=ycount+m;
    ths(vv-integrint+1:vv)=thI(vv-integrint+1:vv)'+Pxy/(Pyy+sigma.^2.*eye(m))*(y(ycount+1-m:ycount)'-FF(1:m,1:m*integrint)*thI(vv-(m*integrint)+1:vv)');
    Pxxs(vv-integrint+1:vv,vv-integrint+1:vv)=PxxI(vv-integrint+1:vv,vv-integrint+1:vv)-Pxy/(Pyy+sigma.^2.*eye(m))*Pxy';
    
    thI=ths;
    PxxI(vv-integrint+1:vv,vv-integrint+1:vv)=Pxxs(vv-integrint+1:vv,vv-integrint+1:vv);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%filtering after the initial condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s=(((maxd/samp)/m)+1):(n/m)
    
    thI=ths;
    PxxI(vv-delaydt+1:vv,vv-delaydt+1:vv)=Pxxs(vv-delaydt+1:vv,vv-delaydt+1:vv);
    
    %loop for the unobserved states until next observation to be predicted
    for j=M:(M+(integrint*m)-1)
        
        %compute weighted input from 0 to tau_m
        thInp=thI((j-delaydt+1):j);
        PFI=thInp*weightsMN';
        
        %compute transcription function and its derivative
        trfun=Rmax/(1+(PFI/Kpc)^Hn);
        trfunDer=-(Rmax*Hn*(PFI/Kpc)^(Hn-1))/(Kpc*(1+(PFI/Kpc)^Hn)^2);
        
        %predict unobserved states mean
        thI(j+1)=A*thI(j)+dt*trfun;
        
        %predict unobserved states variance
        HH=(1-2*mu*dt)*PxxI(vv,vv);
        covmV=dt.*PxxI(vv,vv-delaydt+(1:delaydt))*(trfunDer.*weightsMN');
        covCC=PxxI(vv-delaydt+(1:delaydt),vv-delaydt+(1:delaydt))*(trfunDer.*weightsMN');
        DD=dt*(mu*thI(j)+trfun);
        PxxI(vv+1,vv+1)=HH+2*sum(covmV)+kappa.*DD;
        
        %predict unobserved states variance (back until time of max delay)    
        PxxI(vv-delaydt+(1:delaydt),vv+1)=PxxI(vv-delaydt+(1:delaydt),vv)*A'+dt.*covCC;
        PxxI(vv+1,vv-delaydt+(1:delaydt))=PxxI(vv-delaydt+(1:delaydt),vv+1)';
        
        %updated unobserved states counter
        vv=vv+1;
    end
    
    M=M+integrint*m;

    Pxy=PxxI(vv-delaydt+1:vv,vv-delaydt+1:vv)*(FF(end-m+1:end,:))';
    Pyy=(FF(end-m+1:end,:))*PxxI(vv-delaydt+1:vv,vv-delaydt+1:vv)*(FF(end-m+1:end,:))';
    
    thp(vv-(m*integrint)+1:vv)=thI(vv-(m*integrint)+1:vv)';
    Pxxp(vv-(m*integrint)+1:vv)=diag(PxxI(vv-(m*integrint)+1:vv,vv-(m*integrint)+1:vv));
   
    meanout(:,s)=FF(1:m,1:m*integrint)*thI(vv-(m*integrint)+1:vv)';
    varout(:,:,s)=Pyy+sigma.^2.*eye(m);
    
    ycount=ycount+m;
    ths(vv-delaydt+1:vv)=thI(vv-delaydt+1:vv)'+Pxy/(Pyy+sigma.^2.*eye(m))*(y(ycount+1-m:ycount)'-FF(1:m,1:m*integrint)*thI(vv-(m*integrint)+1:vv)');
    Pxxs(vv-delaydt+1:vv,vv-delaydt+1:vv)=PxxI(vv-delaydt+1:vv,vv-delaydt+1:vv)-Pxy/(Pyy+sigma.^2.*eye(m))*Pxy';
    
end

end