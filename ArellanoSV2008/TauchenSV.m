function [sY,sS, Pi]=TauchenSV(muY,muS,rhoY,rhoS,sigS,NY,NS,mS,mY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function discretize a continuous log AR(1) process by Tauchen's
% method. The log AR(1) process is:
% y_t = muY + rhoY*y_{t-1} + exp(s_t)*e_t, e_t~N(0,1)
% The volatility process is:
% s_t = muS + rhoS*s_{t-1} + u_t, u_t~N(0,sigS^2)
% I used the function from Jan Hannes Lang's website
% This version: 18.03.2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First step: discretize the SV process
sS  = zeros(NS,1);
sS(1)   = muS/(1-rhoS) - mS*sqrt(sigS^2/(1-rhoS^2));
sS(NS)   = muS/(1-rhoS) + mS*sqrt(sigS^2/(1-rhoS^2));
stepS    = (sS(NS)-sS(1))/(NS-1);
for i=2:(NS-1)
   sS(i) = sS(i-1) + stepS; 
end
% Second step: take the last value
SIG=exp(sS');
sigY=SIG(NS);
% Third step: discretize the Y process with sigY
sY  = zeros(NY,1);
sY(1)   = muY/(1-rhoY) - mY*sqrt(sigY^2/(1-rhoY^2));
sY(NY)   = muY/(1-rhoY) + mY*sqrt(sigY^2/(1-rhoY^2));
stepY    = (sY(NY)-sY(1))/(NY-1);
for i=2:(NY-1)
   sY(i) = sY(i-1) + stepY; 
end
% Fourth step: obtaining the matrix P for each volatility level
for j = 1:NY
    for k = 1:NY
        for m=1:NS
            if k == 1
                Pi(j,k,m) = cdf_normal((sY(1) - muY - rhoY*sY(j) + stepY/2) / SIG(m));
            elseif k == NY
                Pi(j,k,m) = 1 - cdf_normal((sY(NY) - muY - rhoY*sY(j) - stepY/2) / SIG(m));
            else
                Pi(j,k,m) = cdf_normal((sY(k) - muY - rhoY*sY(j) + stepY/2) / SIG(m)) - ...
                      cdf_normal((sY(k) - muY - rhoY*sY(j) - stepY/2) / SIG(m));
            end
        end
    end
end

function c = cdf_normal(x)
    c = 0.5 * erfc(-x/sqrt(2));
