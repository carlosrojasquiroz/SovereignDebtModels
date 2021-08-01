%---------------------------------------------------------------------------------------------
% Replication code for Arellano (2008)
%---------------------------------------------------------------------------------------------
% (c) Carlos Rojas Quiroz
%---------------------------------------------------------------------------------------------
% This file requires the following scripts:
% Tauchen.m (from Jan Hannes Lang' website)
% Utility.m
%---------------------------------------------------------------------------------------------
tic
clear all;
close all;
clc;
%---------------------------------------------------------------------------------------------
% Calibration
%---------------------------------------------------------------------------------------------
beta=0.953;
sigma=2.0;
r=0.017;
rho=0.945;
mue=0;
eta=0.025;
theta=0.282;
%---------------------------------------------------------------------------------------------
% Grids for B and Y
%---------------------------------------------------------------------------------------------
Nb=251;
Ny=21;
B=linspace(-0.4,0.0,Nb);
[Y,P]=Tauchen(mue,rho,eta,Ny);
Y=exp(Y');
Yhat=min(0.969*mean(Y),Y);
%---------------------------------------------------------------------------------------------
% Value functions
%---------------------------------------------------------------------------------------------
V=zeros(Nb,Ny);
Vc=zeros(Nb,Ny);
Vd=zeros(Ny,1);
q=ones(Nb,Ny).*(1/(1+r));
policy=zeros(Nb,Ny);
def=zeros(Nb,Ny);
%---------------------------------------------------------------------------------------------
% Parameters for the iteration
%---------------------------------------------------------------------------------------------
iterV=0;
distV=10;
tol=0.00001;

while distV>tol
    Vcomp=V;
%---------------------------------------------------------------------------------------------
% Default Value
%---------------------------------------------------------------------------------------------
zero_ind=find(B>=0,1);
Vd=Utility(sigma,Yhat)'+beta*P*(theta*V(zero_ind,:)'+(1-theta)*Vd);
%---------------------------------------------------------------------------------------------
% Repayment and Total Values
%---------------------------------------------------------------------------------------------
for j=1:Ny
    for i=1:Nb            
        c=max(Y(1,j)-B'.*q(:,j)+B(1,i),0);
        m=Utility(sigma,c)+beta*(P(j,:)*V')';
        [Vc(i,j),policy(i,j)]=max(m(:));
        V(i,j)=(Vd(j,1)>Vc(i,j))*Vd(j,1)+(1-(Vd(j,1)>Vc(i,j)))*Vc(i,j);
        def(i,j)=(Vd(j,1)>Vc(i,j));
    end
end
%---------------------------------------------------------------------------------------------
% Prices
%---------------------------------------------------------------------------------------------
q=(1-def*P')./(1+r);
distV=max(max(abs(V-Vcomp)./((abs(V)+abs(Vcomp))./2+0.001)));
iterV=iterV+1;
end
toc
%---------------------------------------------------------------------------------------------
%% Simulations 
%---------------------------------------------------------------------------------------------
% Number of individuals
Nid=300;
% Number of periods
Nt=10000;
% Initial indices
zero_index=find(B==0); 
Y_init_ind =randi(Ny,Nid,1); 
B_init_ind = randi(Nb,Nid,1);
% Markov Chain for Y
Y_sim_indices=MarkovChain(Nt+1,Nid,P);
%---------------------------------------------------------------------------------------------
% Matrices
%---------------------------------------------------------------------------------------------
Y_sim_val = zeros(Nid,Nt+1);
B_sim_val = zeros(Nid,Nt+1); 
q_sim_val = zeros(Nid,Nt+1);
B_sim_indices = zeros(Nid,Nt+1);
default_status = zeros(Nid,Nt+1);
C_sim_val = zeros(Nid,Nt+1); 
TB_sim_val = zeros(Nid,Nt+1); 
Spread_sim_val = zeros(Nid,Nt+1); 
B_sim_indices(:,1) = B_init_ind;
default_status(:,1) = 0;
Y_sim_val(:,1) = Y(Y_init_ind);
B_sim_va(:,1) = B(B_init_ind);
%---------------------------------------------------------------------------------------------
% Algorithm
%---------------------------------------------------------------------------------------------
for i=1:Nid
    for t=1:Nt
        % get today's indexes
        Yi = Y_sim_indices(i,t); 
        Bi = B_sim_indices(i,t);
        defstat = default_status(i,t);
        % if you are not in default
        if defstat==0
            default_today = Vc(Bi, Yi) < Vd(Yi);
            if default_today
                % default values
                default_status(i,t) = 1;
                default_status(i,t + 1) = 1;
                Y_sim_val(i,t) = Yhat(Y_sim_indices(i,t));
                B_sim_indices(i,t+1) = zero_index;
                B_sim_val(i,t+1) = 0;
                q_sim_val(i,t) = q(zero_index, Y_sim_indices(i,t));
                C_sim_val(i,t) = max(Y_sim_val(i,t)-B_sim_val(i,t+1)*q_sim_val(i,t)+B_sim_val(i,t),0);
                TB_sim_val(i,t) = Y_sim_val(i,t)-C_sim_val(i,t);
                Spread_sim_val(i,t)=((1/q_sim_val(i,t))^(4)-(1+r)^4)*100;
            else
                default_status(i,t) = 0;
                Y_sim_val(i,t) = Y(Y_sim_indices(i,t));
                B_sim_indices(i,t+1) = policy(Bi, Yi);
                B_sim_val(i,t+1) = B(B_sim_indices(i,t+1));
                q_sim_val(i,t) = q(B_sim_indices(i,t+1), Y_sim_indices(i,t));
                C_sim_val(i,t) = max(Y_sim_val(i,t)-B_sim_val(i,t+1)*q_sim_val(i,t)+B_sim_val(i,t),0);
                TB_sim_val(i,t) = Y_sim_val(i,t)-C_sim_val(i,t);
                Spread_sim_val(i,t)=((1/q_sim_val(i,t))^(4)-(1+r)^4)*100;
            end
        % if you are in default
        else
            B_sim_indices(i,t+1) = zero_index;
            B_sim_val(i,t+1) = 0;
            Y_sim_val(i,t) = Y(Y_sim_indices(i,t));
            q_sim_val(i,t) = q(zero_index, Y_sim_indices(i,t));
            C_sim_val(i,t) = max(Y_sim_val(i,t)-B_sim_val(i,t+1)*q_sim_val(i,t)+B_sim_val(i,t),0);
            TB_sim_val(i,t) = Y_sim_val(i,t)-C_sim_val(i,t);
            Spread_sim_val(i,t)=((1/q_sim_val(i,t))^(4)-(1+r)^4)*100;
            % with probability 'theta' exit default status
            default_status(i,t+1) = (rand() >= theta);
        end
    end
end
Y_sim_val=Y_sim_val(:,1:Nt); 
B_sim_val=B_sim_val(:,1:Nt); 
q_sim_val=q_sim_val(:,1:Nt);
C_sim_val=C_sim_val(:,1:Nt);
TB_sim_val=TB_sim_val(:,1:Nt);
default_status=default_status(:,1:Nt); 
%---------------------------------------------------------------------------------------------
% Computation of statistics
%---------------------------------------------------------------------------------------------
ii=1;
Bas=1000;
Tsample=74;
    for j=1:Nid
        Aux=default_status(j,Bas+1:end);
        idef=find(Aux==1,1,'first');
        if idef-Tsample>0
            Ysimul(ii,1:Tsample)=Y_sim_val(j,idef-Tsample:idef-1);
            Bsimul(ii,1:Tsample)=B_sim_val(j,idef-Tsample:idef-1);
            Qsimul(ii,1:Tsample)=q_sim_val(j,idef-Tsample:idef-1);
            Csimul(ii,1:Tsample)=C_sim_val(j,idef-Tsample:idef-1);
            TBsimul(ii,1:Tsample)=TB_sim_val(j,idef-Tsample:idef-1);
            Spreadsimul(ii,1:Tsample)=Spread_sim_val(j,idef-Tsample:idef-1);
            ii=ii+1;
        elseif isempty(idef)==1
            Ysimul(ii,1:Tsample)=Y_sim_val(j,1:Tsample);
            Bsimul(ii,1:Tsample)=B_sim_val(j,1:Tsample);
            Qsimul(ii,1:Tsample)=q_sim_val(j,1:Tsample);
            Csimul(ii,1:Tsample)=C_sim_val(j,1:Tsample);
            TBsimul(ii,1:Tsample)=TB_sim_val(j,1:Tsample);
            Spreadsimul(ii,1:Tsample)=Spread_sim_val(j,1:Tsample);
            ii=ii+1;   
        end
    end
B_Ysimul=Bsimul(:,2:end)./Ysimul(:,1:end-1)*100;
TB_Ysimul=TBsimul./Ysimul*100;
Ysimul=log(Ysimul)*100;
Csimul=log(Csimul)*100;
for jj=1:size(Ysimul,1)
    Ysimul(jj,:)=(detrend(Ysimul(jj,:)'))';
    Csimul(jj,:)=(detrend(Csimul(jj,:)'))';
end
Statistics(1,1)=mean(std(Spreadsimul,0,2));
Statistics(2,1)=mean(std(TB_Ysimul,0,2));
Statistics(3,1)=mean(std(Csimul,0,2));
Statistics(4,1)=mean(std(Ysimul,0,2));
[R,C]=size(Ysimul);
for j=1:R
Aux2(j,1)=corr(Spreadsimul(j,:)',Ysimul(j,:)');
end
Statistics(1,2)=mean(Aux2);
for j=1:R
Aux2(j,1)=corr(TB_Ysimul(j,:)',Ysimul(j,:)');
end
Statistics(2,2)=mean(Aux2);
for j=1:R
Aux2(j,1)=corr(Csimul(j,:)',Ysimul(j,:)');
end
Statistics(3,2)=mean(Aux2);
for j=1:R
Aux3(j,1)=corr(TB_Ysimul(j,:)',Spreadsimul(j,:)');
end
Statistics(2,3)=mean(Aux3);
for j=1:R
Aux3(j,1)=corr(Csimul(j,:)',Spreadsimul(j,:)');
end
Statistics(3,3)=mean(Aux3);
for j=1:R
Aux3(j,1)=corr(Ysimul(j,:)',Spreadsimul(j,:)');
end
Statistics(4,3)=mean(Aux3);
%---------------------------------------------------------------------------------------------
%% Graphs 
%---------------------------------------------------------------------------------------------
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
%---------------------------------------------------------------------------------------------
% 1) Bond price schedule
%---------------------------------------------------------------------------------------------
figure(1)
plot(B(1,[iBlow:iBhigh]),q([iBlow:iBhigh],iYlow),'Color',[48, 71, 94]/255,'LineWidth',2.5);
hold on;
plot(B(1,[iBlow:iBhigh]),q([iBlow:iBhigh],iYhigh),'Color',[234, 84, 85]/255,'LineWidth',2.5);
legend('y_{Low}','y_{High}','FontName','Serif','FontSize', 12,'Location','Best')
xlabel('B´','FontName','Serif','FontSize', 12)
title('Bond price schedule','FontName','Serif','FontSize', 14)
saveas(gcf,'BondPrice','epsc');
%---------------------------------------------------------------------------------------------
% Values for B
highB=0.15;
lowB=-0.30;
% Indexes
iBhigh=find(B>=highB,1);
iBlow=find(B>=lowB,1);
%---------------------------------------------------------------------------------------------
% 2) Value function
%---------------------------------------------------------------------------------------------
figure(2)
plot(B(1,[iBlow:iBhigh]),V([iBlow:iBhigh],iYlow),'Color',[48, 71, 94]/255,'LineWidth',2.5);
hold on;
plot(B(1,[iBlow:iBhigh]),V([iBlow:iBhigh],iYhigh),'Color',[234, 84, 85]/255,'LineWidth',2.5);
legend('y_{Low}','y_{High}','FontName','Serif','FontSize', 12,'Location','Best')
xlabel('B','FontName','Serif','FontSize', 12)
title('Value function','FontName','Serif','FontSize', 14)
saveas(gcf,'ValFn','epsc');
%---------------------------------------------------------------------------------------------
% 3) Savings function
%---------------------------------------------------------------------------------------------
figure(3)
plot(B(1,[iBlow:iBhigh]),B(1,policy([iBlow:iBhigh],iYlow)),'Color',[48, 71, 94]/255,'LineWidth',2.5);
hold on;
plot(B(1,[iBlow:iBhigh]),B(1,policy([iBlow:iBhigh],iYhigh)),'Color',[234, 84, 85]/255,'LineWidth',2.5);
legend('y_{Low}','y_{High}','FontName','Serif','FontSize', 12,'Location','Best')
xlabel('B','FontName','Serif','FontSize', 12)
title('Savings function','FontName','Serif','FontSize', 14)
saveas(gcf,'SavsFn','epsc');
%---------------------------------------------------------------------------------------------
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
%---------------------------------------------------------------------------------------------
% 4) Equilibrium interest rate
%---------------------------------------------------------------------------------------------
figure(4)
plot(B(1,[iBlow:iBhigh]),Rlow,'Color',[48, 71, 94]/255,'LineWidth',2.5);
hold on;
plot(B(1,[iBlow:iBhigh]),Rhigh,'Color',[234, 84, 85]/255,'LineWidth',2.5);
legend('y_{Low}','y_{High}','FontName','Serif','FontSize', 12,'Location','Best')
xlabel('B','FontName','Serif','FontSize', 12)
title('Equilibrium interest rate','FontName','Serif','FontSize', 14)
saveas(gcf,'EqInR','epsc');
