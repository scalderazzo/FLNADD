fsz = 8;    
lwt = 1;  
msz = 8; 
alw= 0.75;
hl=4.7;

load('CryDataMCMC.mat')

names={'log(\kappa R_{max})','log(\kappa K_{pc})','log(n)','log(\mu)','E[\tau]','SD[\tau]','log(\kappa \nu_1)','log(\kappa \nu_2)','log(\kappa \beta)','log(\mu_{0})','log(\kappa E[X(0)])','log(\kappa V[X(0)])','log(\kappa)','log(\sigma_{\epsilon})'};
xs=out;
rep=[1./exp(xs(:,1)),exp((xs(:,2)+xs(:,1))./exp(xs(:,3))),exp(xs(:,3)),exp(xs(:,4)),xs(:,5),xs(:,5)./xs(:,6),exp(xs(:,7)),exp(xs(:,8)),exp(xs(:,9)),exp(xs(:,10)),exp(xs(:,11)),exp(xs(:,12)),exp(xs(:,13)),exp(xs(:,14))];
rep=[log(rep(:,1:4)),rep(:,5:6),log(rep(:,7:end))];

for rr=[1:4]
subplot(4,2,rr)
hold off
[xs,ks]=ksdensity(rep(100000:100:end,rr));
plot(ks,xs,'k','LineWidth',lwt,'MarkerSize',msz)
title(names(rr))
set(gca,'fontSize',fsz,'fontName','CMU Serif','LineWidth', alw) 
set(0,'DefaultAxesTitleFontWeight','normal')
hold on

x=get(gca,'xlim');
veclen=length(out(100000:100:end,rr));
randnumbp=normrnd(a(rr),b(rr),1,veclen.*1000);
[val,xi]=ksdensity(randnumbp);
plot(xi,val,'--','LineWidth',lwt,'MarkerSize',msz);
set(gca,'fontSize',fsz,'fontName','CMU Serif','LineWidth', alw) 
set(0,'DefaultAxesTitleFontWeight','normal')
xlim(x)
end

for rr=[5:6]
subplot(4,2,rr)
hold off
[xs,ks]=ksdensity(rep(100000:100:end,rr));
plot(ks,xs,'k','LineWidth',lwt,'MarkerSize',msz)
title(names(rr))
set(gca,'fontSize',fsz,'fontName','CMU Serif','LineWidth', alw) 
set(0,'DefaultAxesTitleFontWeight','normal')
hold on
axis manual
x=get(gca,'xlim');
l=line(x,[1/(b(rr)-a(rr)) 1/(b(rr)-a(rr))]);
set(l,'LineStyle','--','LineWidth',lwt,'MarkerSize',msz)
set(gca,'fontSize',fsz,'fontName','CMU Serif','LineWidth', alw) 
set(0,'DefaultAxesTitleFontWeight','normal')
xlim(x)
end
for rr=[13:14]
subplot(4,2,rr-6)
hold off
[xs,ks]=ksdensity(rep(100000:100:end,rr));
plot(ks,xs,'k','LineWidth',lwt,'MarkerSize',msz)
title(names(rr))
set(gca,'fontSize',fsz,'fontName','CMU Serif','LineWidth', alw) 
set(0,'DefaultAxesTitleFontWeight','normal')
hold on

x=get(gca,'xlim');
veclen=length(out(100000:100:end,rr));
randnumbp=normrnd(a(rr),b(rr),1,veclen.*1000);
[val,xi]=ksdensity(randnumbp);
plot(xi,val,'--','LineWidth',lwt,'MarkerSize',msz);
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

export_fig postCry.pdf -m4
