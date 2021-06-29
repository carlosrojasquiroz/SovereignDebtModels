%---------------------------------------------------
% Replication code for Cuadra, Sanchez and Sapriza (2010)
%---------------------------------------------------
% (c) Carlos Rojas Quiroz
%---------------------------------------------------
% This file requires the following scripts:
% Tauchen.m (from Jan Hannes Lang' website)
% Utility.m
% This version: 29.06.21
%---------------------------------------------------
tic
clear all;
close all;
clc;
%---------------------------------------------------
% Calibration
%---------------------------------------------------
beta=0.97;
sigma=2.0;
psi=0.455;
theta=0.10;
r=0.01;
rho=0.85;
phi=0.99;
pi=0.30;
mue=0;
eta=0.006;
%---------------------------------------------------
% Grids for B and Y
%---------------------------------------------------
Nb=251;
Ny=25;
B=linspace(-0.15,0.0,Nb);
[A,P]=Tauchen(mue,rho,eta,Ny);
A=exp(A');
Ahat=min(phi*mean(A),A);
%---------------------------------------------------
% Value functions
%---------------------------------------------------
V=zeros(Nb,Ny);
Vc=zeros(Nb,Ny);
Vd=zeros(Ny,1);
q=ones(Nb,Ny).*(1/(1+r));
policy=zeros(Nb,Ny);
def=zeros(Nb,Ny);
%---------------------------------------------------
% Parameters for the iteration
%---------------------------------------------------
iterV=0;
distV=10;
distQ=10;
distVd=10;
tol=0.00001;
TmaxD=TmaxVal(Ahat,0,0,q,sigma,pi,psi)';
TmaxR=TmaxVal(A,B,B,q,sigma,pi,psi);
%---------------------------------------------------
while distV>tol && distQ>tol && distVd>tol
    Vcomp=V;
    Qcomp=q;
    Vdcomp=Vd;
%---------------------------------------------------
% Default Value
%---------------------------------------------------
for j=1:Ny
    ld(1,j)=(Ahat(1,j)/(1+TmaxD(1,j))).^(1/psi);
    cd(1,j)=Ahat(1,j).*ld(1,j)./(1+TmaxD(1,j));
    gd(1,j)=max(TmaxD(1,j).*cd(1,j),0);
end
zero_ind=find(B>=0,1);
Vd=Utility2(sigma,pi,psi,cd,ld,gd)'+beta*P*(theta*V(zero_ind,:)'+(1-theta)*Vd);
%---------------------------------------------------
% Repayment and Total Values
%---------------------------------------------------
for j=1:Ny
    for i=1:Nb
        for h=1:Nb
            lr(j,i,h)=(A(1,j)/(1+TmaxR(j,i,h))).^(1/psi);
            cr(j,i,h)=A(1,j).*lr(j,i,h)./(1+TmaxR(j,i,h));
            gr(j,i,h)=max(TmaxR(j,i,h).*cr(j,i,h)-B(1,h)*q(h,j)+B(1,i),0);
            m(j,i,h)=Utility2(sigma,pi,psi,cr(j,i,h),lr(j,i,h),gr(j,i,h))+beta*(P(j,:)*V(h,:)');
        end
    end
end
for j=1:Ny
    for i=1:Nb
        [Vc(i,j),policy(i,j)]=max(m(j,i,:));
        V(i,j)=(Vd(j,1)>Vc(i,j))*Vd(j,1)+(1-(Vd(j,1)>Vc(i,j)))*Vc(i,j);
        def(i,j)=(Vd(j,1)>Vc(i,j));    
    end
end
%---------------------------------------------------
% Prices
%---------------------------------------------------
q=((1-def)*P')./(1+r);
distV=max(max(abs(V-Vcomp)./((abs(V)+abs(Vcomp))./2+0.001)));
distQ=max(max(abs(q-Qcomp)./((abs(q)+abs(Qcomp))./2+0.001)));
distVd=max(max(abs(Vd-Vdcomp)./((abs(Vd)+abs(Vdcomp))./2+0.001)));
iterV=iterV+1;
end
toc
%---------------------------------------------------
% Replicating figures from Arellano (2008)
%---------------------------------------------------
% 1) Bond price schedule
% Values for Y
highY=mean(A)+std(A);
lowY=mean(A)-std(A);
iYhigh=16; %find(A>=highY,1)-1;
iYlow=10; %find(A>=lowY,1)-1;
% Values for B
highB=0;
lowB=-0.15;
% Indexes 
iBhigh=find(B>=highB,1);
iBlow=find(B>=lowB,1);
figure(1)
plot(B(1,[iBlow:iBhigh]),q([iBlow:iBhigh],iYlow),'Color',[48, 71, 94]/255,'LineWidth',2.5);
hold on;
plot(B(1,[iBlow:iBhigh]),q([iBlow:iBhigh],iYhigh),'Color',[234, 84, 85]/255,'LineWidth',2.5);
legend('y_{Low}','y_{High}','FontName','Serif','FontSize', 12,'Location','Best')
xlabel('B´','FontName','Serif','FontSize', 12)
title('Bond price schedule','FontName','Serif','FontSize', 14)
% 2) Savings function
% Values for B
highB=0;
lowB=-0.09;
% Indexes 
iBhigh=find(B>=highB,1);
iBlow=find(B>=lowB,1);
gB=[iBlow:iBhigh];
for i=1:length(gB)
    if def(gB(i),iYlow)==1
        BprimeL(1,i)=0;
    else
        BprimeL(1,i)=B(1,policy(gB(i),iYlow));
    end
end
for i=1:length(gB)
    if def(gB(i),iYhigh)==1
        BprimeH(1,i)=0;
    else
        BprimeH(1,i)=B(1,policy(gB(i),iYhigh));
    end
end
figure(2)
plot(B(1,[iBlow:iBhigh]),BprimeL,'Color',[48, 71, 94]/255,'LineWidth',2.5);
hold on;
plot(B(1,[iBlow:iBhigh]),BprimeH,'Color',[234, 84, 85]/255,'LineWidth',2.5);
legend('y_{Low}','y_{High}','FontName','Serif','FontSize', 12,'Location','Best')
xlabel('B','FontName','Serif','FontSize', 12)
title('Savings function','FontName','Serif','FontSize', 14)
% 3) Tax function
% Values for B
highB=0;
lowB=-0.042;
% Indexes 
iBhigh=find(B>=highB,1);
iBlow=find(B>=lowB,1);
gB=[iBlow:iBhigh];
for k=1:length(gB)
    if def(gB(k),iYlow)==1
        tL(1,k)=TmaxD(k);
    else
        tL(1,k)=TmaxR(iYlow,gB(k),policy(gB(k),iYlow));
    end
end
for k=1:length(gB)
    if def(gB(k),iYhigh)==1
        tH(1,k)=TmaxD(k);
    else
        tH(1,k)=TmaxR(iYhigh,gB(k),policy(gB(k),iYhigh));
    end
end
figure(3)
plot(B(1,[iBlow:iBhigh]),tL,'Color',[48, 71, 94]/255,'LineWidth',2.5);
hold on;
plot(B(1,[iBlow:iBhigh]),tH,'Color',[234, 84, 85]/255,'LineWidth',2.5);
legend('y_{Low}','y_{High}','FontName','Serif','FontSize', 12,'Location','Best')
xlabel('B','FontName','Serif','FontSize', 12)
title('Tax rate function','FontName','Serif','FontSize', 14)
