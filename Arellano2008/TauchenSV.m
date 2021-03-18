function [sY,sS, Pi]=TauchenSV(muY,muS,rhoY,rhoS,sigS,NY,NS,mS,mY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%       Discretize an AR(1) process with log AR(1) stochastic volatility
%       y_t = lambda*y_{t-1} + u_t
%       x_t = (1-rho)*mu + rho*x_{t-1} + epsilon_t
%       u_t ~ N(0,exp(x_t)); epsilon_t ~ N(0,sigma_e^2)
%
% Usage:
%       [P,yxGrids] = discreteSV(lambda,rho,sigmaU,sigmaE,Ny,Nx,method,nSigmaY)
%
% Inputs:
%       lambda    - persistence of y process
%       rho       - persistence of x process
%       sigmaU    - unconditional standard deviation of u_t
%       sigmaE    - standard deviation of epsilon_t
%       Ny        - number of grid points for y process
%       Nx        - number of grid points for x process
% This version: 18.03.2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First step: discretize the SV process
sS  = zeros(NS,1);
PiS = zeros(NS,NS);
sS(1)   = muS/(1-rhoS) - mS*sqrt(sigS^2/(1-rhoS^2));
sS(NS)   = muS/(1-rhoS) + mS*sqrt(sigS^2/(1-rhoS^2));
step    = (sS(NS)-sS(1))/(NS-1);
for i=2:(NS-1)
   sS(i) = sS(i-1) + step; 
end
% Second step: take the last value
SIG=exp(sS');
sigY=SIG(NS);
% Third step: discretize the Y process with sigY
sY  = zeros(NY,1);
PiS = zeros(NY,NY);
sY(1)   = muY/(1-rhoY) - mY*sqrt(sigY^2/(1-rhoY^2));
sY(NY)   = muY/(1-rhoY) + mY*sqrt(sigY^2/(1-rhoY^2));
step    = (sY(NY)-sY(1))/(NY-1);
for i=2:(NY-1)
   sY(i) = sY(i-1) + step; 
end
% Fourth step: obtaining the matrix P for each volatility level
for j = 1:NY
    for k = 1:NY
        for m=1:NS
            if k == 1
                Pi(j,k,m) = cdf_normal((sY(1) - muY - rhoY*sY(j) + step/2) / SIG(m));
            elseif k == NY
                Pi(j,k,m) = 1 - cdf_normal((sY(NY) - muY - rhoY*sY(j) - step/2) / SIG(m));
            else
                Pi(j,k,m) = cdf_normal((sY(k) - muY - rhoY*sY(j) + step/2) / SIG(m)) - ...
                      cdf_normal((sY(k) - muY - rhoY*sY(j) - step/2) / SIG(m));
            end
        end
    end
end

function c = cdf_normal(x)
    c = 0.5 * erfc(-x/sqrt(2));