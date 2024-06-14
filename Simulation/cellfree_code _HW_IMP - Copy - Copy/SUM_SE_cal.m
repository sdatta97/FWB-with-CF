function SIM_SE = SIM_SE_cal(G,eta,p,mmse_combiner)
V = zeros(1,K);
Delta = zeros(1,K);
SINR = zeros(1,K);
for ue = 1:K
    for AP =1:M
        q(:,AP,ue)=eye(2)*rand(2,1);
    end
end
for ue = 1:K
    for AP =1:M
        for ue1=1:K
            if (ue1==ue)
                V(ue) = p(m,k)*eta(m,k)*norm(transpose(G(:,AP,ue))*mmse_combiner(:,AP,ue)*q(:,AP,ue))^2 +V(ue);
            else
                Delta(ue) =p(m,k)*eta(m,k)*norm(transpose(G(:,AP,ue))*mmse_combiner(:,AP,ue1)*q(:,AP,ue))^2 +V(ue);
            end
    end
    end
    Delta(ue) = Delta(ue) + eye(2);
    SINR(ue) = log2(det(eye(2) + V(ue)*pinv(Delta(ue))));
end
SIM_SE = sum(SINR);

