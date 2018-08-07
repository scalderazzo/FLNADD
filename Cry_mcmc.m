function [] = Cry_mcmc

load('Cry1Data.mat');
Cry1DN=Cry1NormDetrended';

maxNumCompThreads(1);

rng(11); %location ID in SCN

t=0.5:0.5:116.5;
t1=1:233;
tt=1:233;
[yy,gof,out]=fit(t1',Cry1DN,'smoothingspline','SmoothingParam',0.3); 
Gilldatasm=feval(yy,tt);
timechange=(t(find(Gilldatasm(1:39)==max(Gilldatasm(1:39)))));

a=zeros(1,14);
b=10.*ones(1,14);

b([9,12:14])=20; 

a(14)=-5;
b(14)=1;

a(3)=log(1.5);
b(3)=5;

a([4,10])=log(0.58);
b([4,10])=0.5; 

b(5)=23;
b(6)=20;

up=Inf.*ones(1,14);
low=-Inf.*ones(1,14);
up(5)=23;
up(6)=20;
low(5)=0;
low(6)=0;

th0=zeros(1,14);
th0([1:4,7:14])=normrnd(a([1:4,7:14]),b([1:4,7:14]));
th0(5:6)=unifrnd(a(5:6),b(5:6));
maxd=30;
m=1;
deltat=0.5;
it=350000;
cry=Cry1DN';

[xout,errorout,outlogL] = MCMC_Cry_DA(it,deltat,cry,th0,a,b,up,low,maxd,timechange,m);

end