function rate_dl = compute_link_rates_MIMO(params,channel_dl, channel_est_dl,ue_idx,sub6ConnectionState)
a = 1;
b = 1;
M = size(channel_dl,1);
Kmax = params.numUE_sub6_max;
Kd = size(channel_dl,2);
Kd_mmw = size(sub6ConnectionState,1);
K = params.num_sc_sub6; % no of subcarriers
scs = params.scs_sub6; %sub carrier spacing
BW = K*scs;
Ntx = size(channel_dl,3);
p_d = 1*Kd;
% perm_vec  = repmat(randperm(tau_p),1,2);
% phi_index = perm_vec(1:K);
% for k = 1:K
%     PHI(:,k) = PHI1(:,phi_index(k));
% end
% gamma_num = zeros(M,K);
% gamma_den = zeros(M,K);
% Gamma = zeros(M,K);
% for m = 1:M
%     for k = 1:K
%         gamma_num(m,k) = tau_p*p_p*(BETA(m,k)^2);
%         gamma_den_temp = zeros(1,K);
%         for j = 1:K
%             gamma_den_temp(j) = BETA(m,j)*(abs(PHI(:,j)'*PHI(:,k))^2);
%         end
% 
%         gamma_den(m,k) = tau_p*p_p*sum(gamma_den_temp)+1;
% 
%         Gamma(m,k) = gamma_num(m,k)/gamma_den(m,k);
%     end
% end
BETA = params.BETA;
[BETA_max, idx_max] = maxk(BETA,Kmax,2);
%% initialization of c
C_v = repmat(sqrt(1./(b*Ntx*sum(BETA(idx_max),2))),[1 Kd]);
idx_ue = zeros(M,Kd);
for i = 1:M
    for k = 1:Kd
        if any(idx_max(i,:) == k)
            idx_ue(i,k) = 1;
        end
    end
end
D = zeros(Kd,M);
for k = 1:Kd
    D(k,:) = reshape(sqrt(C_v(:,k).*sum(reshape(channel_dl(:,k,:),[M,Ntx]).*reshape(conj(channel_est_dl(:,k,:)),[M,Ntx]),2)),[1,M]).*(idx_ue(:,k)');
end
% E = zeros(M,Kd);
DS = zeros(Kd,1);
MUI = zeros(Kd,1);
% % UDI = zeros(Kd,1);
% % TQD = zeros(Kd,1);
N = abs(sqrt(0.5)*(randn(Kd,1) + 1j*randn(Kd,1))).^2;
snr_num = zeros(Kd,1);
snr_den = zeros(Kd,1);
rate_dl = zeros(1,Kd);
for k = 1:Kd_mmw
    if (sub6ConnectionState(k)==1 || k==ue_idx)
        DS(k) = a^2*p_d*abs(D(k,:)*D(k,:).')^2;
        MUI(k) = 0;
        for q = 1:Kd_mmw
             if (q~=k && sub6ConnectionState(q)==1)
               MUI(k) = MUI(k) + a^2*p_d*abs((C_v(:,q).*(idx_ue(:,q)==1))'*(reshape(sum(reshape(channel_dl(:,k,:),[M,Ntx]).*reshape(conj(channel_est_dl(:,q,:)),[M,Ntx]),2),[M,1])))^2;
             end
        end
        for q = 1+Kd_mmw:Kd
            % E(:,q) = sqrt(0.5*(b-a^2)*(C_v(:,q).^2))*(randn(M,1) + 1i*randn(M,1));
            % MUI(k) = MUI(k) + abs(C_v(:,q)'*reshape(sum(reshape(channel_dl(:,k,:),[M,Ntx]).*reshape(conj(channel_est_dl(:,q,:)),[M,Ntx]),2),[M,1]))^2;
            MUI(k) = MUI(k) + a^2*p_d*abs((C_v(:,q).*(idx_ue(:,q)==1))'*(reshape(sum(reshape(channel_dl(:,k,:),[M,Ntx]).*reshape(conj(channel_est_dl(:,q,:)),[M,Ntx]),2),[M,1]).*(idx_ue(:,q)==1)))^2;
        end
        % for l = 1:Ku
        %     UDI(k) = UDI(k) + abs(interUEchannel(k,l)*sqrt(p_u*Theta_v(l)))^2;
        % end
        snr_num(k) = DS(k);
        snr_den(k) = MUI(k) + N(k);
        rate_dl(k) = BW*TAU_FAC*log2(1+snr_num(k)/snr_den(k));
    end
end
for k = 1+Kd_mmw:Kd
    DS(k) = a^2*p_d*abs(D(k,:)*D(k,:).')^2;
    MUI(k) = 0;
    for q = 1:Kd_mmw
        if (sub6ConnectionState(q)==1)
           MUI(k) = MUI(k) + a^2*p_d*abs((C_v(:,q).*(idx_ue(:,q)==1))'*(reshape(sum(reshape(channel_dl(:,k,:),[M,Ntx]).*reshape(conj(channel_est_dl(:,q,:)),[M,Ntx]),2),[M,1])))^2;
        end
    end
    for q = 1+Kd_mmw:Kd
       if (q~=k)
           MUI(k) = MUI(k) + a^2*p_d*abs((C_v(:,q).*(idx_ue(:,q)==1))'*(reshape(sum(reshape(channel_dl(:,k,:),[M,Ntx]).*reshape(conj(channel_est_dl(:,q,:)),[M,Ntx]),2),[M,1])))^2;
       end
    end
    snr_num(k) = DS(k);
    snr_den(k) = MUI(k) + N(k);
    rate_dl(k) = BW*TAU_FAC*log2(1+snr_num(k)/snr_den(k));
end
end