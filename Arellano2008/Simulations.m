Nid=300;
Nt=10000;
% get initial indices
zero_index=find(B==0); 
Y_init = median(Y);
B_init = median(B);
Y_init_ind =randi(Ny,Nid,1); % find(Y==Y_init); %
B_init_ind = randi(Nb,Nid,1); %find(B==B_init); %
% create a QE MarkovChain
Y_sim_indices=MarkovChain(Nt+1,Nid,P);
% allocate and fill output
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

%for j=1:Nid
%Y_sim_indices(j,:)=markovchain2(P,Nt+1)';
%end

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
