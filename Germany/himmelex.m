%%
 
% The system is written in file <himmelode.html |himmelode.m|> and the
% sum of squares function in <himmelss.html |himmelss.m|>.
clc
clear
clear model data parama options
num=xlsread('Germany data.xlsx');
data.ydata=[num(1:26,1),num(46:71,6:7)];

Iq0 = 586;  E0 = 20; A0 = 100; I0 = 70; 
R10 = 18;  Eq0 = 100; R20=20; Aq0 = 195;


%%
% Initial concentrations are saved in |data| to be used in sum of
% squares function.
Pop = 8.2*10^7    % Region actual total population
L = 100;            % Number of iterations
          

Iqvalue = zeros(L,1001); R1value=zeros(L,1001);
value=zeros(L,5);

for i=1:L
%%
% Refine the first guess for the parameters with |fmincon| and
% calculate residual variance as an estimate of the model error variance.
eta0 = 0.9; theta0 = 0.9; rho0 = 0.9; beta0 = 0.5;
phi0  = 1.9068*10^(-5); epsilon0 = 0.5;  alpha0 = 0.22;
gammaA0 = 0.0001; gamma0 = 0.73; gammaq0 = 0.035;
Na = 0.0025; mu = 0.23; zeta = 0.035;

k00 = [eta0,theta0,rho0,beta0,phi0,epsilon0,alpha0,gammaA0,gamma0,gammaq0,Na,mu,E0,A0,I0,Eq0,R20,Aq0,zeta]';

lb=[0.5,0.1,0.1,0.1,10^(-6),0.1,0.1,0,0,0,0,0,0,0,0,0,0,0,0]';
ub=[0.95,0.95,0.95,0.9,10^(-3),0.9,1,1,1,1,1,1,10^4,10^4,10^4,10^4,10^4,10^4,1]';
[k0,ss0] = fmincon(@(k0) himmelss(k00,data),k00,[],[],[],[],lb,ub);
mse = ss0/length(data.ydata);

%%
params = {
    {'k1', k0(1), 0,0.95}
    {'k2', k0(2), 0,0.95}
    {'k3', k0(3), 0,1}
    {'k4', k0(4), 0}
    {'k5', k0(5), 0}
    {'k6', k0(6), 0}
    {'k7', k0(7), 0}
    {'k8', k0(8), 0}
    {'k9', k0(9), 0}
    {'k10', k0(10), 0.03}
    {'k11', k0(11), 0}
    {'k12', k0(12), 0}
    {'k13', k0(13), 0}
    {'k14', k0(14), 0}
    {'k15', k0(15), 0}
    {'k16', k0(16), 0}
    {'k17', k0(17), 0}
    {'k18', k0(18), 0}
    {'k19', k0(19), 0}
    };
model.ssfun = @himmelss;
model.sigma2 = mse;

options.nsimu = 2000;
options.updatesigma = 1;

%%
[results,chain,s2chain] = mcmcrun(model,data,params,options);

%%
% figure(1); clf
% mcmcplot(chain,[],results,'chainpanel')
% subplot(5,4,19)
% mcmcplot(sqrt(s2chain),[],[],'dens',2)
% title('error std')

%%
% Function |chainstats| lists some statistics, including the
% estimated Monte Carlo error of the estimates.
chainstats(chain,results)

%%
kk=mean(chain);
data.y0 = [num(46,6)-kk(18);kk(13);kk(14);kk(15);R10;kk(16);kk(17);kk(18)];
[t,y] = ode45(@himmelode,linspace(0,100,1001),data.y0,[],mean(chain));

% figure(2);%SEIR
% plot(10*num(1:26,1),num(46:71,6),'o','Color',[0.92 0.5 0.44])
% hold on
% plot(10*num(1:26,1),num(46:71,7),'o','Color',[0.3 0.8 0.9])

IIq=y(:,1)+y(:,8);
% hold on
% plot(IIq,'-','Color',[0.8 0.4 0.1])%Iq
% hold on
% plot(y(:,2),'y-');%E
% hold on
% plot(y(:,3),'m-')%A
% hold on
% plot(y(:,4),'b-')%I
% hold on
% plot(y(:,5),'g-')%R1
% hold on
% plot(y(:,6),'c-')%Eq
% hold on
% plot(y(:,7),'c-')%R2
% 
% legend({'Iqtrue','Rtrue','Iq','E','A','I','R1','Eq','R2'},'Location','best')
% ylim([0,1.5*10^5])
% title('Germany')
% xlabel('Days from 7 March,2020')
% ylabel('Population')

lamda=(kk(6)/kk(7)+kk(1)/(kk(2)+kk(9))+...
    kk(4)*(1-kk(1))/(kk(12)+kk(8)))*(1-kk(3))*kk(5)*kk(11)*Pop; %Basic reproductive number

[aa,index]=max(IIq);
PA=(y(index,3)+y(index,8))/(y(index,1)+y(index,3)+y(index,4)+y(index,8));
% PB=y(index,8)/(y(index,3)+y(index,8))
[Iqq,index1]=max(IIq)
VV=[lamda,PA,Iqq,index1,1/kk(7)];
value(i,:)=VV;

Iqvalue(i,:)=IIq'; R1value(i,:)=y(:,5)'; 

end

Iqmm=prctile(Iqvalue,[2.5 97.5],1); R1mm=prctile(R1value,[2.5 97.5],1);

figure(3);clf
plot(t',Iqmm,'-','Color',[0.96 0.96 0.96],'HandleVisibility','off')
hold on
plot(t',R1mm,'-','color',[0.93 0.89 0.87],'HandleVisibility','off')

ylim([0,1.5*10^5])
set(gca,'XTickLabel',[0:10:100])
title('Germany')
xlabel('Days from 7 March,2020')
ylabel('Population')

hold on
f=zeros(length(t),2);
j=1;
for j=1:1001
    f(j,1)=max([Iqmm(1,j),Iqmm(2,j)]);
    f(j,2)=min([Iqmm(1,j),Iqmm(2,j)]);
    
end

fy=cat(2,f(:,1)',flipdim(f(:,2),1)')';
fx=cat(2,t',flipdim(t,1)')';
H_F1=patch(fx,fy,[0.93 0.91 0.91])
set(H_F1,{'EdgeColor'},{'[0.93 0.91 0.91]'})


hold on
ff=zeros(length(t),2);
jj=1;
for jj=1:1001
    ff(jj,1)=max([R1mm(1,jj),R1mm(2,jj)]);
    ff(jj,2)=min([R1mm(1,jj),R1mm(2,jj)]);
    
end

ffy=cat(2,ff(:,1)',flipdim(ff(:,2),1)')';
ffx=cat(2,t',flipdim(t,1)')';
H_F2=patch(ffx,ffy,[0.93 0.83 0.71],'FaceAlpha',0.5)
set(H_F2,{'EdgeColor'},{'[0.93 0.83 0.71]'})
hold on
plot(num(1:26,1),num(46:71,6),'o','Color',[0.92 0.5 0.44])
hold on
plot(num(1:26,1),num(46:71,7),'o','Color',[0.3 0.8 0.9])
% legend({'Iq predicted value','R1 predicted value','Iq true data','R1 true data'},'Location','best')
set(0,'defaultfigurecolor','w')

%%
valuemm=prctile(value,[2.5 97.5],1);
BRN=mean(value);
Med=median(value);
