function Tmax=TmaxVal(A,B,Bprime,q,sigma,pi,psi)
Ny=length(A);
Nb=length(B);
funU=zeros(Ny,Nb,Nb);
Tmax=zeros(Ny,Nb,Nb);

for j=1:Ny % Productivity
        for i=1:Nb % Current debt
            for h=1:Nb % Future debt
                funU=@(T)-(pi*(max(T*(A(1,j)/(1+T)*((A(1,j)/(1+T))^(1/psi)))+B(1,i)-q(h,j)*Bprime(1,h),0)).^(1-sigma)./(1-sigma)+(1-pi)*(((A(1,j)/(1+T)*((A(1,j)/(1+T))^(1/psi)))-((A(1,j)/(1+T))^(1/psi)).^(1+psi)./(1+psi)).^(1-sigma)./(1-sigma)));
                T0=0.168;
                Tmax(j,i,h)=fminbnd(funU,0,1);
            end
        end
end
