%-----------------------------------------------------------------------------------------------------
% Sovereign debt and fiscal policy in commodity-exporting economies
% Master thesis - UC3M Master in Economic Analysis
%-----------------------------------------------------------------------------------------------------
% (c) Carlos Rojas Quiroz
%-----------------------------------------------------------------------------------------------------
% This file requires the following scripts:
% discreteVAR.m (from Farmer and Toda, 2017)
% Utility.m
% This version: July 28th, 2021
%-----------------------------------------------------------------------------------------------------
tic
clear all; close all; clc;
% If you want to see the VAR estimation, de-comment the following line
% run VarEstimation.m;
%-----------------------------------------------------------------------------------------------------
% Calibration (parameters)
%-----------------------------------------------------------------------------------------------------
beta=0.964; % discount factor
sigma=2.0; % risk aversion
psi=0.4; % inverse of Frisch elasticity
theta=0.10; % re-entry probability
r=0.01; % free risk interest rate
phi=0.99; % default cost
pi=0.30; % gov. consumption on hh.'s preferences
tauP=0.0774; % average commodity income tax rate
%-----------------------------------------------------------------------------------------------------
% Discretizing B
%-----------------------------------------------------------------------------------------------------
Nb=150;
B=linspace(-0.4,0.0,Nb);
%-----------------------------------------------------------------------------------------------------
% Discretizing P and A
%-----------------------------------------------------------------------------------------------------
% I assume P and A have a VAR representation:
% P = C0(1) + B0(1,1) P(-1) + B0(1,2) A(-1) + e1
% A = C0(2) + B0(2,1) P(-1) + B0(2,2) A(-1) + e2
% Parameters obtained from run VarEstimation.m
C0=[0; 0];
B0=[0.9548 0; 0.0526 0.9650];
sigTOT=0.0433;
sigY=0.0116;
covTOTY=0.0100;
COVAR0=[sigTOT^2 0; 0 sigY^2];
Nstates=11;
[Pr,Xr] = discreteVAR(C0,B0,COVAR0,Nstates);
P=exp(Xr(1,:));
A=exp(Xr(2,:));
Ahat=min(phi*mean(A),A);
Phat=min(phi*mean(P),P);
%-----------------------------------------------------------------------------------------------------
% Calibration (commodity output)
%-----------------------------------------------------------------------------------------------------
Amean=mean(A);
tauC=0.15;
Lmean=(Amean/(1+tauC))^(1/psi);
ymean=Amean*Lmean; 
shyX=0.25;
yX=ymean*shyX/(1-shyX);
%-----------------------------------------------------------------------------------------------------
% Value functions
%-----------------------------------------------------------------------------------------------------
Ns=Nstates^2;
V=zeros(Nb,Ns);
Vc=zeros(Nb,Ns);
Vd=zeros(Ns,1);
q=ones(Nb,Ns).*(1/(1+r));
policy=zeros(Nb,Ns);
def=zeros(Nb,Ns);
%-----------------------------------------------------------------------------------------------------
% Parameters for the iteration
%-----------------------------------------------------------------------------------------------------
iterV=0;
distV=10;
distQ=10;
distVd=10;
tol=10e-06;
TmaxD=TmaxVal(Ahat,yX,Phat,0,0,q,sigma,pi,psi,tauP)';
TmaxR=TmaxVal(A,yX,P,B,B,q,sigma,pi,psi,tauP);
%-----------------------------------------------------------------------------------------------------
while distV>tol && distQ>tol && distVd>tol
    Vcomp=V;
    Qcomp=q;
    Vdcomp=Vd;
%-----------------------------------------------------------------------------------------------------
% Default Value
%-----------------------------------------------------------------------------------------------------
    for j=1:Ns
        ld(1,j)=(Ahat(1,j)/(1+TmaxD(1,j))).^(1/psi);
        cd(1,j)=(Ahat(1,j).*ld(1,j)+(1-tauP)*yX*Phat(1,j))./(1+TmaxD(1,j));
        gd(1,j)=max(TmaxD(1,j).*cd(1,j)+tauP*yX*Phat(1,j),0);
    end
zero_ind=find(B>=0,1);
Vd=Utility2(sigma,pi,psi,cd,ld,gd)'+beta*Pr*(theta*V(zero_ind,:)'+(1-theta)*Vd);
%-----------------------------------------------------------------------------------------------------
% Repayment and Total Values
%-----------------------------------------------------------------------------------------------------
    for j=1:Ns
        for i=1:Nb
            for h=1:Nb
            lr(j,i,h)=(A(1,j)/(1+TmaxR(j,i,h))).^(1/psi);
            cr(j,i,h)=(A(1,j).*lr(j,i,h)+(1-tauP)*yX*P(1,j))./(1+TmaxR(j,i,h));
            gr(j,i,h)=max(TmaxR(j,i,h).*cr(j,i,h)+tauP*yX*P(1,j)-B(1,h)*q(h,j)+B(1,i),0);
            m(j,i,h)=Utility2(sigma,pi,psi,cr(j,i,h),lr(j,i,h),gr(j,i,h))+beta*(Pr(j,:)*V(h,:)');
            end
        end
    end
    for j=1:Ns
        for i=1:Nb
        [Vc(i,j),policy(i,j)]=max(m(j,i,:));
        V(i,j)=(Vd(j,1)>Vc(i,j))*Vd(j,1)+(1-(Vd(j,1)>Vc(i,j)))*Vc(i,j);
        def(i,j)=(Vd(j,1)>Vc(i,j));    
        end
    end
%-----------------------------------------------------------------------------------------------------
% Prices
%-----------------------------------------------------------------------------------------------------
q=(1-def*Pr')./(1+r);
distV=max(max(abs(V-Vcomp)./((abs(V)+abs(Vcomp))./2+0.001)));
distQ=max(max(abs(q-Qcomp)./((abs(q)+abs(Qcomp))./2+0.001)));
distVd=max(max(abs(Vd-Vdcomp)./((abs(Vd)+abs(Vdcomp))./2+0.001)));
iterV=iterV+1;
end
toc
