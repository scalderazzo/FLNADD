function [xmat,trfun]= SSA(T,dt,c0,Rmax,Kpc,Hn,mu,meanD,sdD,maxd,Pre,Post)

%T is the total simulation time length
%dt is the sampling interval
%c0 is the initial number of molecules of each species

timegrid=0:dt:(maxd-dt);
timegridH=(timegrid+timegrid+dt)./2;

it=T/dt;
i=1;
target=dt;
timespan=length(timegrid);
trfun=zeros(1,timespan);
xmat=c0;
tt=0;
S=(Post-Pre)'; %stochiometry

aP=meanD^2/sdD^2;
bP=sdD^2/meanD;
pdP = makedist('Gamma','a',aP,'b',bP);
gammatrP = truncate(pdP,0,maxd);
evalp=maxd-timegridH;   
weightsM=pdf(gammatrP,evalp); 
weightsN=weightsM./(sum(weightsM));

x=xmat(:,end);

while i<=it 
    h=[Rmax/(1+((xmat(timespan+i-length(timegrid):timespan+i-1)*weightsN')./Kpc).^Hn),mu*x];
    tt=tt+exprnd(1/sum(h));
    j=gendist(h./sum(h),1,1);
    x=x+(S(:,j));
    while tt>=target && i<=it
        xmat(:,timespan+i)=x;
        trfun(timespan+i)=h(1);
        i=i+1;
        target=target+dt; 
    end
end
end
