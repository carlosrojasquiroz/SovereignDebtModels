% Simulation
Nid=100;
Nt=250;
% Ys(id,t)
Ry=rand(Nid,Nt);
Cy=zeros(Ny,Ny);
% Financial exclusion
Re=rand(Nid,Nt);
Es=zeros(Nid,Nt);
% Simulation matrix for default decision and savings
defs=zeros(Nid,Nt);
Bprimes=zeros(Nid,Nt);
qs=zeros(Nid,Nt);
% Initial period
IndY=randi(Ny,Nid,1); % random assingments
IndB=randi(Nb,Nid,1); % random assingments
Ys=Y(IndY)';
Bs=B(IndB)';

for j=1:Nt
for i=1:Nid
    if j==1
        defs(i,j)=def(IndB(i,j),IndY(i,j));
        Bprimes(i,j)=B(policy(IndB(i,j),IndY(i,j)));
        qs(i,j)=((1-def(IndB(i,j),:))*P(IndY(i,j),:)')./(1+r);
    else 
        % Index of financial exclusion
        if defs(i,j-1)==1 || Es(i,j-1)==1
            if Re(i,j-1)>theta
                Es(i,j)=1;
            else
                Es(i,j)=0;
            end
        end
        % Output
        for k=1:Ny
            for l=1:Ny
                Cy(k,l)=sum(P(k,(1:l)));
            end
        end
        comp=10;
            for m=1:Ny
                Aux1(m,1)=Cy(IndY(i),m)-Ry(i,j);
                    if  Aux1(m,1)<0
                        Aux1(m,1)=100000;
                    end
                    if Aux1(m,1)<=comp
                        comp=Aux1(m,1);
                        jtilde=m;
                    end
            end
            Ys(i,j)=Y(jtilde);
    if Es(i,j)==1
        defs(i,j)=1;
        Bprimes(i,j)=B(zero_ind);
    elseif Es(i,j)==0
        IndB(i,j)=find(Bprimes(i,j-1)==B);
        IndY(i,j)=find(Ys(i,j)==Y);
        defs(i,j)=def(IndB(i,j),IndY(i,j));
        Bprimes(i,j)=B(policy(IndB(i,j),IndY(i,j)));
        qs(i,j)=((1-def(IndB(i,j),:))*P(IndY(i,j),:)')./(1+r);
        end
    end
end
end
