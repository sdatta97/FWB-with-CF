function SIM_SE = SIM_SE_cal(G,eta,p,mmse_combiner)
V = zeros(1,K);
Delta = zeros(1,K);
SINR = zeros(1,K);
for ue = 1:K
    for AP =1:M
        q(:,AP,ue)=1/sqrt(2)*eye(2)*rand(2,1);
    end
end
D=zeros(1,AP);
for AP=1:ap
    
    for ue=1:UE
        D(AP)= D(AP)+ p(AP,ue)*(q(:,ue,AP)'*(G(:,ue,AP)'*G(:,ue,AP))*q(:,ue,AP));
    end
    eta(AP) = rho_d*pinv(D(AP));
end

for ue = 1:K
    for AP =1:M
        for ue1=1:K
            if (ue1==ue)
                V(ue) = p(AP,ue)*eta(AP,ue)*(transpose(G(:,ue,AP))*mmse_combiner(:,ue,AP)*q(:,AP,ue))*(transpose(G(:,ue,AP))*mmse_combiner(:,ue,AP)*q(:,AP,ue))' +V(ue);
            else
                Delta(ue) =p(AP,ue)*eta(AP,ue)*(transpose(G(:,ue,AP))*mmse_combiner(:,ue1,AP)*q(:,AP,ue))*(transpose(G(:,ue,AP))*mmse_combiner(:,ue1,AP)*q(:,AP,ue))' +Delta(ue);
            end
        end
    end
    Delta(ue) = Delta(ue) + eye(2);
    SINR(ue) = log2(det(eye(2) + V(ue)*pinv(Delta(ue))));
end
SIM_SE = sum(SINR);

