%---------------------------------------
% Replication code for Arellano (2008)
%---------------------------------------
% (c) Carlos Rojas Quiroz
%---------------------------------------
% This file requires the following scripts:
% Tauchen.m (from Jan Hannes Lang' website)
% Utility.m
%---------------------------------------
tic
clear all;
close all;
clc;
%---------------------------------------
% Calibration
%---------------------------------------
beta=0.953;
sigma=2.0;
r=0.017;
rho=0.945;
mue=0;
eta=0.025;
theta=0.282;
%---------------------------------------
% Grids for B and Y
%---------------------------------------
Nb=251;
Ny=21;
Nt=200;
B=linspace(-0.4,0.4,Nb);
[Y,P]=Tauchen(mue,rho,eta,Ny);
Y=exp(Y');
Yhat=min(0.969*mean(Y),Y);
%---------------------------------------
% Value functions
%---------------------------------------
V=zeros(Nb,Ny);
Vc=zeros(Nb,Ny);
Vd=zeros(Ny,1);
q=ones(Nb,Ny).*(1/(1+r));
policy=zeros(Nb,Ny);
def=zeros(Nb,Ny);
%---------------------------------------
% Parameters for the iteration
%---------------------------------------
iterV=0;
distV=10;
tol=0.00001;
%---------------------------------------
t=1;
% Default value
Vd=Utility(sigma,Yhat)';
% Repayment and Total Values
for j=1:Ny
    for i=1:Nb            
        c=max(Y(1,j)+B(1,i),0);
        m=Utility(sigma,c);
        [Vc(i,j),policy(i,j)]=max(m(:));
        V(i,j)=(Vd(j,1)>Vc(i,j))*Vd(j,1)+(1-(Vd(j,1)>Vc(i,j)))*Vc(i,j);
        def(i,j)=(Vd(j,1)>Vc(i,j));
    end
end
%---------------------------------------
for t=2:Nt
    Vcomp=V;
%---------------------------------------
% Default Value
%---------------------------------------
zero_ind=find(B>=0,1);
Vd=Utility(sigma,Yhat)'+beta*P*(theta*V(zero_ind,:)'+(1-theta)*Vd);
%---------------------------------------
% Repayment and Total Values
%---------------------------------------
for j=1:Ny
    for i=1:Nb            
        c=max(Y(1,j)-B'.*q(:,j)+B(1,i),0);
        m=Utility(sigma,c)+beta*(P(j,:)*V')';
        [Vc(i,j),policy(i,j)]=max(m(:));
        V(i,j)=(Vd(j,1)>Vc(i,j))*Vd(j,1)+(1-(Vd(j,1)>Vc(i,j)))*Vc(i,j);
        def(i,j)=(Vd(j,1)>Vc(i,j));
    end
end
%---------------------------------------
% Prices
%---------------------------------------
q=((1-def)*P')./(1+r);
distV=max(max(abs(V-Vcomp)./((abs(V)+abs(Vcomp))./2+0.001)));
end
toc

%---------------------------------------
% Replicating figures from Arellano (2008)
%---------------------------------------
% Values for Y
highY=1.05*mean(Y);
lowY=0.95*mean(Y);
% Values for B
highB=0;
lowB=-0.35;
% Indexes 
iYhigh=find(Y>=highY,1)-1;
iYlow=find(Y>=lowY,1)-1;
iBhigh=find(B>=highB,1);
iBlow=find(B>=lowB,1);
% 1) Bond price schedule
figure(1)
plot(B(1,[iBlow:iBhigh]),q([iBlow:iBhigh],iYlow),'Color',[48, 71, 94]/255,'LineWidth',2.5);
hold on;
plot(B(1,[iBlow:iBhigh]),q([iBlow:iBhigh],iYhigh),'Color',[234, 84, 85]/255,'LineWidth',2.5);
legend('y_{Low}','y_{High}','FontName','Serif','FontSize', 12,'Location','Best')
xlabel('B´','FontName','Serif','FontSize', 12)
title('Bond price schedule','FontName','Serif','FontSize', 14)
saveas(gcf,'BondPrice02','epsc');
%---------------------------------------
% Values for B
highB=0.15;
lowB=-0.30;
% Indexes
iBhigh=find(B>=highB,1);
iBlow=find(B>=lowB,1);
% 2) Value function
figure(2)
plot(B(1,[iBlow:iBhigh]),V([iBlow:iBhigh],iYlow),'Color',[48, 71, 94]/255,'LineWidth',2.5);
hold on;
plot(B(1,[iBlow:iBhigh]),V([iBlow:iBhigh],iYhigh),'Color',[234, 84, 85]/255,'LineWidth',2.5);
legend('y_{Low}','y_{High}','FontName','Serif','FontSize', 12,'Location','Best')
xlabel('B','FontName','Serif','FontSize', 12)
title('Value function','FontName','Serif','FontSize', 14)
saveas(gcf,'ValFn02','epsc');
%---------------------------------------
% 3) Savings function
figure(3)
plot(B(1,[iBlow:iBhigh]),B(1,policy([iBlow:iBhigh],iYlow)),'Color',[48, 71, 94]/255,'LineWidth',2.5);
hold on;
plot(B(1,[iBlow:iBhigh]),B(1,policy([iBlow:iBhigh],iYhigh)),'Color',[234, 84, 85]/255,'LineWidth',2.5);
legend('y_{Low}','y_{High}','FontName','Serif','FontSize', 12,'Location','Best')
xlabel('B','FontName','Serif','FontSize', 12)
title('Savings function','FontName','Serif','FontSize', 14)
saveas(gcf,'SavsFn02','epsc');
%---------------------------------------
% Values for B
highB=0.03;
lowB=-0.08;
% Indexes
iBhigh=find(B>=highB,1);
iBlow=find(B>=lowB,1);
% Interest rates (annual): NaN if default
grid=iBlow:iBhigh;
for i=1:length(grid)
    if def(grid(i),iYlow)==1
        Rlow(i,1)=NaN;
    else
        Rlow(i,1)=((q(policy(grid(i),iYlow),iYlow).^(-1)).^4)-1;
    end
    if def(grid(i),iYhigh)==1
        Rhigh(i,1)=NaN;
    else
        Rhigh(i,1)=((q(policy(grid(i),iYhigh),iYhigh).^(-1)).^4)-1;
    end
end
% 4) Equilibrium interest rate
figure(4)
plot(B(1,[iBlow:iBhigh]),Rlow,'Color',[48, 71, 94]/255,'LineWidth',2.5);
hold on;
plot(B(1,[iBlow:iBhigh]),Rhigh,'Color',[234, 84, 85]/255,'LineWidth',2.5);
legend('y_{Low}','y_{High}','FontName','Serif','FontSize', 12,'Location','Best')
xlabel('B','FontName','Serif','FontSize', 12)
title('Equilibrium interest rate','FontName','Serif','FontSize', 14)
saveas(gcf,'EqInR02','epsc');
