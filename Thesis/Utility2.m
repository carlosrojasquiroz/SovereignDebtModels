function U=Utility2(sigma,pi,psi,c,l,g)
U=pi*(g.^(1-sigma)./(1-sigma))+(1-pi)*((c-l.^(1+psi)/(1+psi)).^(1-sigma)./(1-sigma));
end
