function [] = CrySimulationCoverage(Hn)
rng(Hn)

%parameters : transcription function
Rmax=90;   %maximum transcription rate
Kpc=150;    %thhreshold
Hn=Hn;     %Hill coefficient
mu=0.25;     %degradation rate
meanD=9.2;  %mean delay
sdD=sqrt(15);    %standard deviation delay

%parameters : initial condition
ni1=80;    %transcription rate before switch
ni2=30;    %transcription rate after switch
beta=6;   %scale coefficient
muIn=0.2;   %degradation rate

%maximum delay (in hours)
maxd=30;

%observed Cry1-luc
load('Cry1DataPaper.mat')
Cry1DN=Cry1NormDetrended';
 
n=233;                  %data sample size
Dt=0.5;                 %time-interval for observations (in hours)
ot=n*Dt;                %data totoal observation time

%%simulated data

%%Matrices to obatain stochiometry matrix
Pre=[0;
 1];
Post=[1;
 0];

%%obtain frequent data for initial condition from observed normalised
%%Cry1-luc

nSim=n*2;               %sample size of simulated data (double, to forget initial condition)
dtS=0.01;               %time-interval for SSA

in=Cry1DN./mean(Cry1DN);
t1=Dt:Dt:ot;
tt=dtS:dtS:ot;
[yy,~,~]=fit(t1',in,'smoothingspline'); 
insm=feval(yy,tt);

scaling=0.006; %scaling factor to molecule numbers from previous simulation
insmR=insm./scaling;   %rescaled initial condition

%simulate with SSA

for tt=1:10
c0=insmR(1:(maxd/dtS))';
[xmat,trfun]= SSA((nSim*Dt)-maxd,dtS,c0,Rmax,Kpc,Hn,mu,meanD,sdD,maxd,Pre,Post);
gillsimDCS1T(tt,:)=xmat;
gillsimcumIT=cumsum(gillsimDCS1T(tt,:));
gillsimcumS1T(tt,:)=[gillsimcumIT(1).*(Dt/dtS),diff(gillsimcumIT((Dt/dtS):(Dt/dtS):n/dtS))];
size=1/mean(gillsimcumS1T(tt,:));
gillsimdataS1T(tt,:)=size.*gillsimcumS1T(tt,:)+normrnd(0,0.01,1,nSim);
gillsimdataS1HT(tt,:)=size.*gillsimcumS1T(tt,:)+normrnd(0,0.05,1,nSim);
end


%initialise matrices for coverage probabilities
covPo=zeros(10,5,2);
covSo=zeros(10,5,2);

%discretisation time-intervals dt
discr=[0.5,0.25,0.1,0.05,0.01];

%observed data time-vector
t=Dt:Dt:ot;

%%


for discrI=1
    
    
    dt=discr(discrI);
    covP=zeros(ot/dt,10,2);
    covS=zeros(ot/dt,10,2);
    
    for uu=1:10
        
        P0=gillsimDCS1T(uu,ot/dtS+1);
        mu0=gillsimDCS1T(uu,ot/dtS+1);
        
        for sn=2
            
            if (sn==1)
                sigmae=0.01;
                data=gillsimdataS1T(uu,(n+1):end);
            else
                sigmae=0.05;
                data=gillsimdataS1HT(uu,(n+1):end);
            end
            
            
            t1=1:n;
            [yy,~,~]=fit(t1',data','smoothingspline','SmoothingParam',0.3);
            Gilldatasm=feval(yy,t1);
            timechange=(t(find(Gilldatasm(1:39)==max(Gilldatasm(1:39)))));
            
            size=(1/mean(gillsimcumS1T(uu,:)))*(dt/dtS);
            par=[size.*Rmax,size.*Kpc,Hn,mu,meanD,sdD,size*ni1,size*ni2,size.*beta,muIn,size.*mu0,(size.^2).*P0,size,sigmae];
            
            [ths,thp,meanout,varout,Pxxs,Pxxp] = FilterDelay(dt,1,data,maxd,timechange,par);
            
            
            meanoutot=[];
            varoutot=[];
            for rr=1:n
                meanoutot=[meanoutot,meanout(:,rr)'];
                varoutot=[varoutot,diag(varout(:,:,rr))'];
            end
            
            
            varunobp=Pxxp;
            y1p=thp+1.96.*sqrt(varunobp);
            y2p=thp-1.96.*sqrt(varunobp);

            varunob=diag(Pxxs)';
            y1s=ths+1.96.*sqrt(varunob);
            y2s=ths-1.96.*sqrt(varunob);

            unobs=[gillsimDCS1T(uu,(ot/dtS+1):end)];
            unP=size.*unobs(1:(dt/0.01):end);
            
            for yy=((maxd/dt)+1):(ot/dt)
                covP(yy,uu,sn)=(unP(yy)>y2p(yy) & unP(yy)<y1p(yy));
                covS(yy,uu,sn)=(unP(yy)>y2s(yy) & unP(yy)<y1s(yy));
            end
            
            covPo(uu,discrI,sn)=mean(covP(((maxd/dt)+1):end,uu,sn));
            covSo(uu,discrI,sn)=mean(covS(((maxd/dt)+1):end,uu,sn));
            
        end
        
    end
   
    
end

save(['SimCry' num2str(5)])
end