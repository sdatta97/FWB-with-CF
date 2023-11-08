function rate_dl = compute_link_rates_MIMO(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,D,ue_idx,sub6ConnectionState)
M = size(channel_dl,1);
K = size(channel_dl,2);
K_mmW = size(sub6ConnectionState,1);
BW = params.Band;
Ntx = size(channel_dl,3);
N_mmW = size(channel_dl_mmW,4);
N_sub6 = size(channel_dl,4);
p_d = 1*K;
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
beta_uc = zeros(size(BETA));

%Prepare array to store the number of APs serving a specficic UE
La = zeros(K,1);
%Prepare cell to store the AP indices serving a specficic UE
Serv = cell(K,1);
%Prepare cell to store the AP indices not serving a specficic UE
NoServ = cell(K,1);
%Construc the above array and cells
for k = 1:K
    servingAPs = find(D(:,k)==1);
    NoservingAPs = find(D(:,k)==0);
    
    Serv{k} = servingAPs;
    NoServ{k} = NoservingAPs;
    
    La(k) = length(servingAPs);
    beta_uc(:,k) = BETA(:,k).*D(:,k);
end
%% initialization of c
eta_eq = zeros(L,K);
if (K_mmW == 0)
    for l = 1:L
        for k = 1:K
            if ismember(l,Serv{k})
                eta_eq(l,k) = 1./(N_AP*N_UE_sub6*sum(BETA(l,:)));
            end
        end
    end
else
    for l = 1:L
        for k = 1:K
            if ismember(l,Serv{k})
                if (k<=K_mmW)
                    eta_eq(l,k) = p_fac./(N_AP*(N_UE_mmW*p_fac*beta_uc(l,1:K_mmW)+N_UE_sub6*sum(beta_uc(l,2:K))));
                else
                    eta_eq(l,k) = 1./(N_AP*(N_UE_mmW*p_fac*beta_uc(l,1:K_mmW)+N_UE_sub6*sum(beta_uc(l,2:K))));
                end
            end
        end
    end
end
C_v = sqrt(eta_eq);
D_mmW = zeros(K_mmW,Ntx,Ntx);
D_sub6 = zeros(K-K_mmW,Ntx,Ntx);
for k = 1:K_mmW
    for m = 1:M
        D_mmW(k,:,:) = D_mmW(k,:,:) + sqrt(C_v(m,k))*reshape(channel_dl(m,k,:,:),[Ntx,Nrx])*reshape(conj(channel_est_dl(m,k,:,:)),[Ntx,Nrx])';
    end
end
for k = 1:K-K_mmW
    for m = 1:M
        D_sub6(k,:,:) = D_sub6(k,:,:) + sqrt(C_v(m,k+K_mmW))*reshape(channel_dl(m,k,:,:),[Ntx,Nrx])*reshape(conj(channel_est_dl(m,k,:,:)),[Ntx,Nrx])';
    end
end
DS_mmW = zeros(K_mmW,N_mmW);
MUI_mmW = zeros(K_mmW,N_mmW);
DS_sub6 = zeros(K-K_mmW,N_sub6);
MUI_sub6 = zeros(K-K_mmW,N_sub6);

noise_mmW = abs(sqrt(0.5)*(randn(K_mmW,N_mmW) + 1j*randn(K_mmW,N_mmW))).^2;
noise_sub6 = abs(sqrt(0.5)*(randn(K-K_mmW,N_sub6) + 1j*randn(K-K_mmW,N_sub6))).^2;
snr_num_mmW = zeros(K_mmW,N_mmW);
snr_den_mmW = zeros(K_mmW,N_mmW);
snr_num_sub6 = zeros(K-K_mmW,N_sub6);
snr_den_sub6 = zeros(K-K_mmW,N_sub6);
rate_dl = zeros(K,1);
for k = 1:K_mmW
    if (sub6ConnectionState(k)==1 || k==ue_idx)
        for n = 1:N_mmW
            DS_mmW(k,n) = p_d*abs(reshape(D_mmW(k,n,:),[1,N_mmW])*ue_sym)^2;
            for q = 1:K_mmW
                if (q~=k && sub6ConnectionState(q)==1)
                  MUI_mmW(k,n) = MUI_mmW(k,n) + p_d*abs(reshape(D_mmW(q,n,:),[1,N_mmW])*ue_sym)^2;
                end
            end
            for q = 1:K-K_mmW
               MUI_mmW(k,n) = MUI_mmW(k,n) + p_d*abs(reshape(D_sub6(q,n,:),[1,N_sub6])*ue_sym)^2;
            end
            snr_num_mmW(k,n) = DS_mmW(k,n);
            snr_den_mmW(k,n) = MUI_mmW(k,n) + noise_mmW(k,n);
            rate_dl(k) = rate_dl(k) + BW*TAU_FAC*log2(1+snr_num_mmW(k,n)/snr_den_mmW(k,n));
        end
    end
end
for k = 1:K-K_mmW
    for n = 1:N_mmW
        DS_sub6(k,n) = p_d*abs(reshape(D_sub6(k,n,:),[1,N_sub6])*ue_sym)^2;
        for q = 1:K_mmW
            if (q~=k && sub6ConnectionState(q)==1)
              MUI_sub6(k,n) = MUI_sub6(k,n) + p_d*abs(reshape(D_mmW(q,n,:),[1,N_mmW])*ue_sym)^2;
            end
        end
        for q = 1:K-K_mmW
           MUI_sub6(k,n) = MUI_sub6(k,n) + p_d*abs(reshape(D_sub6(q,n,:),[1,N_sub6])*ue_sym)^2;
        end
        snr_num_sub6(k,n) = DS_sub6(k,n);
        snr_den_sub6(k,n) = MUI_sub6(k,n) + noise_sub6(k,n);
        rate_dl(k+K_mmW) = rate_dl(k+K_mmW) + BW*TAU_FAC*log2(1+snr_num_sub6(k,n)/snr_den_sub6(k,n));
    end
end
end