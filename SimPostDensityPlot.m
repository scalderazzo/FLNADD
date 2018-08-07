fsz = 8;   
lwt = 0.6;    
msz = 8; 
alw= 0.75;
hl=4.7;

load('SimCry5.mat')
Hn=5;

uu=1;
t=0.5:0.5:116.5;
P0=[gillsimDCS1T(uu,11650)];
mu0=[gillsimDCS1T(uu,11650)];
size=(1/mean(gillsimcumS1T(uu,:)))*50;
sigmae=0.01;
tv=[size.*90,size.*150,Hn,0.25,9.2,sqrt(15),30,size*200,size*50,10,size.*mu0,size.*P0,size,sigmae];
tvL=[log(tv(1:4)),tv(5:6),log(tv(8:end))];
parnames={'log(\kappa R_{max})','log(\kappa K_{pc})','log(n)','log(\mu)','E[\tau]','SD[\tau]','log(\kappa \nu_1)','log(\kappa \nu_2)','log(\kappa \beta)','log(\mu_{0})','log(\kappa E[X(0)])','log(\kappa V[X(0)])','log(\kappa)','log(\sigma_{\epsilon})'};
xl=[-2,-1,0,-2.2,6,0,-7,-10];
xu=[2,0.3,4,0,18,13,-3,-3];
for rr=[1:4]
subplot(4,2,rr)
hold off
for ii=[1:3,5,6,8:10]
load(['crysimMCMC' num2str(ii)]);
P0=[gillsimDCS1T(ii,11650)];
mu0=[gillsimDCS1T(ii,11650)];
size=(1/mean(gillsimcumS1T(ii,:)))*50;
sigmae=0.01;
tv=[size.*90,size.*150,Hn,0.25,9.2,sqrt(15),30,size*200,size*50,10,0.2,size.*mu0,size.*P0,size,sigmae];
tvL=[log(tv(1:4)),tv(5:6),log(tv(8:end))];
xs=out(100000:100:end,:);
rep=[1./exp(xs(:,1)),exp((xs(:,2)+xs(:,1))./exp(xs(:,3))),exp(xs(:,3)),exp(xs(:,4)),xs(:,5),xs(:,5)./xs(:,6),exp(xs(:,7)),exp(xs(:,8)),exp(xs(:,9)),exp(xs(:,10)),exp(xs(:,11)),exp(xs(:,12)),exp(xs(:,13)),exp(xs(:,14))];
rep=[log(rep(:,1:4)),rep(:,5:6),log(rep(:,7:end))];
[xs,ks]=ksdensity(rep(:,rr));
plot(ks,xs,'k','LineWidth',lwt,'MarkerSize',msz)
xlim([xl(rr),xu(rr)])
if rr==1
ylim([0,2.5])
end
if rr==2
ylim([0,13])
end
y=get(gca,'ylim');
l=line([tvL(rr) tvL(rr)],y);
set(l,'color','c','LineWidth',lwt)
title(parnames(rr))
set(gca,'fontSize',fsz,'fontName','CMU Serif','LineWidth', alw) 
set(0,'DefaultAxesTitleFontWeight','normal')
hold on
end
x=get(gca,'xlim');
veclen=length(out(100000:end,rr));
randnumbp=normrnd(a(rr),b(rr),1,veclen);
[val,xi]=ksdensity(randnumbp);
plot(xi,val,'--','LineWidth',1,'MarkerSize',msz);
xlim(x)
set(gca,'fontSize',fsz,'fontName','CMU Serif','LineWidth', alw) 
set(0,'DefaultAxesTitleFontWeight','normal')
end
for rr=[5:6]
subplot(4,2,rr)
hold off
for ii=[1:3,5,6,8:10]
load(['crysimMCMC' num2str(ii)]);
P0=[gillsimDCS1T(ii,11650)];
mu0=[gillsimDCS1T(ii,11650)];
size=(1/mean(gillsimcumS1T(ii,:)))*50;
sigmae=0.01;
tv=[size.*90,size.*150,Hn,0.25,9.2,sqrt(15),30,size*200,size*50,10,0.2,size.*mu0,size.*P0,size,sigmae];
tvL=[log(tv(1:4)),tv(5:6),log(tv(8:end))];
xs=out(100000:100:end,:);
rep=[1./exp(xs(:,1)),exp((xs(:,2)+xs(:,1))./exp(xs(:,3))),exp(xs(:,3)),exp(xs(:,4)),xs(:,5),xs(:,5)./xs(:,6),exp(xs(:,7)),exp(xs(:,8)),exp(xs(:,9)),exp(xs(:,10)),exp(xs(:,11)),exp(xs(:,12)),exp(xs(:,13)),exp(xs(:,14))];
rep=[log(rep(:,1:4)),rep(:,5:6),log(rep(:,7:end))];
[xs,ks]=ksdensity(rep(:,rr));
plot(ks,xs,'k','LineWidth',lwt,'MarkerSize',msz)
xlim([xl(rr),xu(rr)])
y=get(gca,'ylim');
l=line([tvL(rr) tvL(rr)],y);
set(l,'color','c','LineWidth',lwt)
title(parnames(rr))
set(gca,'fontSize',fsz,'fontName','CMU Serif','LineWidth', alw) 
set(0,'DefaultAxesTitleFontWeight','normal')
hold on
end
axis manual
x=get(gca,'xlim');
l=line(x,[1/(b(rr)-a(rr)) 1/(b(rr)-a(rr))]);
set(l,'LineStyle','--','LineWidth',1)
set(gca,'fontSize',fsz,'fontName','CMU Serif','LineWidth', alw) 
set(0,'DefaultAxesTitleFontWeight','normal')
xlim(x)
end
for rr=[13:14]
subplot(4,2,rr-6)
hold off
for ii=[1:3,5,6,8:10]
load(['crysimMCMC' num2str(ii)]);
P0=[gillsimDCS1T(ii,11650)];
mu0=[gillsimDCS1T(ii,11650)];
size=(1/mean(gillsimcumS1T(ii,:)))*50;
sigmae=0.01;
tv=[size.*90,size.*150,Hn,0.25,9.2,sqrt(15),30,size*200,size*50,10,0.2,size.*mu0,size.*P0,size,sigmae];
tvL=[log(tv(1:4)),tv(5:6),log(tv(8:end))];
xs=out(100000:100:end,:);
rep=[1./exp(xs(:,1)),exp((xs(:,2)+xs(:,1))./exp(xs(:,3))),exp(xs(:,3)),exp(xs(:,4)),xs(:,5),xs(:,5)./xs(:,6),exp(xs(:,7)),exp(xs(:,8)),exp(xs(:,9)),exp(xs(:,10)),exp(xs(:,11)),exp(xs(:,12)),exp(xs(:,13)),exp(xs(:,14))];
rep=[log(rep(:,1:4)),rep(:,5:6),log(rep(:,7:end))];
[xs,ks]=ksdensity(rep(:,rr));
plot(ks,xs,'k','LineWidth',lwt,'MarkerSize',msz)
if rr==14
ylim([0,1.5])
end
xlim([xl(rr-6),xu(rr-6)])
y=get(gca,'ylim');
l=line([tvL(rr) tvL(rr)],y);
set(l,'color','c','LineWidth',lwt)
title(parnames(rr))
hold on
end
x=get(gca,'xlim');
veclen=length(out(100000:end,rr));
randnumbp=normrnd(a(rr),b(rr),1,veclen);
[val,xi]=ksdensity(randnumbp);
plot(xi,val,'--','LineWidth',1,'MarkerSize',msz);
set(gca,'fontSize',fsz,'fontName','CMU Serif','LineWidth', alw) 
set(0,'DefaultAxesTitleFontWeight','normal')
xlim(x)
end


width = 3.47;    
height = hl;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); 
h=gcf;
set(findall(h,'type','text'),'fontSize',fsz,'fontName','CMU Serif')
set(h, 'Color', 'w');

export_fig postsim100.pdf -m2

