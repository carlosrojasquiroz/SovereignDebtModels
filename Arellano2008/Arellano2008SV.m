%---------------------------------------------------
% Arellano's model with stochastic volatility
%---------------------------------------------------
% (c) Carlos Rojas Quiroz
%---------------------------------------------------
% This file requires the following scripts:
% TauchenSV.m
% Utility.m
% This version: 21.03.21
%---------------------------------------------------
tic
clear all;
close all;
clc;
%---------------------------------------------------
% Calibration
%---------------------------------------------------
beta=0.953;
sigma=2.0;
r=0.017;
theta=0.282;
% Output process
muY=0;
rhoY=0.945;
% Stochastic process
rhoS=0.85;
muS=(1-rhoS)*(-3.7);
sigS=0.05;
% Tauchen parameters
mS=3; 
mY=3;
%---------------------------------------------------
% Grids for B and Y
%---------------------------------------------------
Nb=251;
Ny=21;
Ns=11;
B=linspace(-0.4,0.4,Nb);
[Y,S, P]=TauchenSV(muY,muS,rhoY,rhoS,sigS,Ny,Ns,mS,mY);
Y=exp(Y');
S=exp(S');
Yhat=min(0.969*mean(Y),Y);
%---------------------------------------------------
% Value functions
%---------------------------------------------------
V=zeros(Nb,Ny,Ns);
Vc=zeros(Nb,Ny,Ns);
Vd=zeros(Ny,Ns,1);
q=ones(Nb,Ny,Ns).*(1/(1+r));
policy=zeros(Nb,Ny,Ns);
def=zeros(Nb,Ny,Ns);
%---------------------------------------------------
% Parameters for the iteration
%---------------------------------------------------
iterV=0;
distV=10;
tol=0.00001;

while distV>tol
    Vcomp=V;
%---------------------------------------------------
% Default Value
%---------------------------------------------------
zero_ind=find(B>=0,1);
for k=1:Ns
Vd(:,k,1)=Utility(sigma,Yhat)'+beta*P(:,:,k)*(theta*V(zero_ind,:,k)'+(1-theta)*Vd(:,k,1));
end
%---------------------------------------------------
% Repayment and Total Values
%---------------------------------------------------
for k=1:Ns
    for j=1:Ny
        for i=1:Nb            
            c=max(Y(1,j)-B'.*q(:,j,k)+B(1,i),0);
            m=Utility(sigma,c)+beta*(P(j,:,k)*V(:,:,k)')';
            [Vc(i,j,k),policy(i,j,k)]=max(m(:));
            V(i,j,k)=(Vd(j,k,1)>Vc(i,j,k))*Vd(j,k,1)+(1-(Vd(j,k,1)>Vc(i,j,k)))*Vc(i,j,k);
            def(i,j,k)=(Vd(j,k,1)>Vc(i,j,k));
        end
    end
end
%---------------------------------------------------
% Prices
%---------------------------------------------------
for k=1:Ns
q(:,:,k)=((1-def(:,:,k))*P(:,:,k)')./(1+r);
end
distV=max(max(abs(V-Vcomp)./((abs(V)+abs(Vcomp))./2+0.001)));
iterV=iterV+1;
end
toc
%---------------------------------------------------
% Replicating figures from Arellano (2008)
%---------------------------------------------------
% Values for B
highB=0;
lowB=-0.30;
% Indexes 
iYhigh=15;
iYlow=9;
iYmid=13;
iShigh=11;
iSlow=1;
iBhigh=find(B>=highB,1);
iBlow=find(B>=lowB,1);
%---------------------------------------------------
% 1) Bond price schedule
figure(1)
plot(B(1,[iBlow:iBhigh]),q([iBlow:iBhigh],iYlow,iSlow),'Color',[48, 71, 94]/255,'LineWidth',2.5);
xlabel('B´','FontName','Serif','FontSize', 18)
title('Low income','FontName','Serif','FontSize', 18)
hold on;
plot(B(1,[iBlow:iBhigh]),q([iBlow:iBhigh],iYlow,iShigh),'Color',[234, 84, 85]/255,'LineWidth',2.5);
saveas(gcf,'qLow','epsc');
figure(2)
plot(B(1,[iBlow:iBhigh]),q([iBlow:iBhigh],iYmid,iSlow),'Color',[48, 71, 94]/255,'LineWidth',2.5);
xlabel('B´','FontName','Serif','FontSize', 18)
title('Middle income','FontName','Serif','FontSize', 18)
hold on;
plot(B(1,[iBlow:iBhigh]),q([iBlow:iBhigh],iYmid,iShigh),'Color',[234, 84, 85]/255,'LineWidth',2.5);
saveas(gcf,'qMid','epsc');
figure(3)
plot(B(1,[iBlow:iBhigh]),q([iBlow:iBhigh],iYhigh,iSlow),'Color',[48, 71, 94]/255,'LineWidth',2.5);
xlabel('B´','FontName','Serif','FontSize', 18)
title('High income','FontName','Serif','FontSize', 18)
hold on;
plot(B(1,[iBlow:iBhigh]),q([iBlow:iBhigh],iYhigh,iShigh),'Color',[234, 84, 85]/255,'LineWidth',2.5);
legend('\sigma_{Low}','\sigma_{High}','FontName','Serif','FontSize', 20,'Location','Best')
saveas(gcf,'qHigh','epsc');
%---------------------------------------------------
% 2) Savings function
% Values for B
highB=0.15;
lowB=-0.30;
iBhigh=find(B>=highB,1);
iBlow=find(B>=lowB,1);
figure(4)
plot(B(1,[iBlow:iBhigh]),B(1,policy([iBlow:iBhigh],iYlow,iSlow)),'Color',[48, 71, 94]/255,'LineWidth',2.5);
hold on;
plot(B(1,[iBlow:iBhigh]),B(1,policy([iBlow:iBhigh],iYlow,iShigh)),'Color',[234, 84, 85]/255,'LineWidth',2.5);
xlabel('B','FontName','Serif','FontSize', 18)
title('Low income','FontName','Serif','FontSize', 18)
saveas(gcf,'SavsLow','epsc');
figure(5)
plot(B(1,[iBlow:iBhigh]),B(1,policy([iBlow:iBhigh],iYmid,iSlow)),'Color',[48, 71, 94]/255,'LineWidth',2.5);
hold on;
plot(B(1,[iBlow:iBhigh]),B(1,policy([iBlow:iBhigh],iYmid,iShigh)),'Color',[234, 84, 85]/255,'LineWidth',2.5);
xlabel('B','FontName','Serif','FontSize', 18)
title('Middle income','FontName','Serif','FontSize', 18)
saveas(gcf,'SavsMid','epsc');
figure(6)
plot(B(1,[iBlow:iBhigh]),B(1,policy([iBlow:iBhigh],iYhigh,iSlow)),'Color',[48, 71, 94]/255,'LineWidth',2.5);
hold on;
plot(B(1,[iBlow:iBhigh]),B(1,policy([iBlow:iBhigh],iYhigh,iShigh)),'Color',[234, 84, 85]/255,'LineWidth',2.5);
xlabel('B','FontName','Serif','FontSize', 18)
title('High income','FontName','Serif','FontSize', 18)
legend('\sigma_{Low}','\sigma_{High}','FontName','Serif','FontSize', 20,'Location','Best')
saveas(gcf,'SavsHigh','epsc');
%---------------------------------------------------
% Values for B
highB=0.4;
lowB=-0.4;
% Indexes
iBhigh=find(B>=highB,1);
iBlow=find(B>=lowB,1);
% Interest rates (annual): NaN if default
grid=iBlow:iBhigh;
for i=1:length(grid)
    if def(grid(i),iYlow,iSlow)==1
        Rlow(i,1)=NaN;
    else
        Rlow(i,1)=((q(policy(grid(i),iYlow,iSlow),iYlow,iSlow).^(-1)).^4)-1;
    end
    if def(grid(i),iYlow,iShigh)==1
        Rhigh(i,1)=NaN;
    else
        Rhigh(i,1)=((q(policy(grid(i),iYlow,iShigh),iYlow,iShigh).^(-1)).^4)-1;
    end
end
% 4a) Equilibrium interest rate
figure(7)
plot(B(1,[iBlow:iBhigh]),Rhigh,'Color',[234, 84, 85]/255,'LineWidth',2.5);
hold on;
plot(B(1,[iBlow:iBhigh]),Rlow,'Color',[48, 71, 94]/255,'LineWidth',2.5);
xlabel('B','FontName','Serif','FontSize', 18)
title('Low income','FontName','Serif','FontSize', 18)
saveas(gcf,'EqInRLow','epsc');
% Interest rates (annual): NaN if default
grid=iBlow:iBhigh;
for i=1:length(grid)
    if def(grid(i),iYmid,iSlow)==1
        Rlow(i,1)=NaN;
    else
        Rlow(i,1)=((q(policy(grid(i),iYmid,iSlow),iYmid,iSlow).^(-1)).^4)-1;
    end
    if def(grid(i),iYmid,iShigh)==1
        Rhigh(i,1)=NaN;
    else
        Rhigh(i,1)=((q(policy(grid(i),iYmid,iShigh),iYmid,iShigh).^(-1)).^4)-1;
    end
end
% 4b) Equilibrium interest rate
figure(8)
plot(B(1,[iBlow:iBhigh]),Rhigh,'Color',[234, 84, 85]/255,'LineWidth',2.5);
hold on;
plot(B(1,[iBlow:iBhigh]),Rlow,'Color',[48, 71, 94]/255,'LineWidth',2.5);
xlabel('B','FontName','Serif','FontSize', 18)
title('Middle income','FontName','Serif','FontSize', 18)
saveas(gcf,'EqInRMid','epsc');
% Interest rates (annual): NaN if default
grid=iBlow:iBhigh;
for i=1:length(grid)
    if def(grid(i),iYhigh,iSlow)==1
        Rlow(i,1)=NaN;
    else
        Rlow(i,1)=((q(policy(grid(i),iYhigh,iSlow),iYhigh,iSlow).^(-1)).^4)-1;
    end
    if def(grid(i),iYhigh,iShigh)==1
        Rhigh(i,1)=NaN;
    else
        Rhigh(i,1)=((q(policy(grid(i),iYhigh,iShigh),iYhigh,iShigh).^(-1)).^4)-1;
    end
end
% 4c) Equilibrium interest rate
figure(9)
plot(B(1,[iBlow:iBhigh]),Rlow,'Color',[48, 71, 94]/255,'LineWidth',2.5);
hold on;
plot(B(1,[iBlow:iBhigh]),Rhigh,'Color',[234, 84, 85]/255,'LineWidth',2.5);
xlabel('B','FontName','Serif','FontSize', 18)
title('High income','FontName','Serif','FontSize', 18)
legend('\sigma_{Low}','\sigma_{High}','FontName','Serif','FontSize', 20,'Location','Best')
saveas(gcf,'EqInRHigh','epsc');
