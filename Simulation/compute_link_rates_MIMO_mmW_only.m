function rate_dl = compute_link_rates_MIMO_mmW_only(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,ue_idx,sub6ConnectionState)
M = size(channel_dl,1);
K_mmW = size(sub6ConnectionState,1);
K = K_mmW + size(channel_dl,2);
BW = params.Band;
TAU_FAC = params.preLogFactor;
Ntx = size(channel_dl,3);
N_UE_mmW = size(channel_dl_mmW,4);
N_UE_sub6 = size(channel_dl,4);
p_d = params.rho_tot; % 1*K;
D = params.D;
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
eta_eq = zeros(M,1);
N_AP = params.num_antennas_per_gNB;
if ((K_mmW == 0) || (sub6ConnectionState == zeros(K_mmW,1)))
    rate_dl = 0;
else
    for m = 1:M
       if ismember(m,Serv{ue_idx})
            if ((ue_idx<=K_mmW) && (sub6ConnectionState(ue_idx) == 1))
                eta_eq(m) = 1/(N_AP*N_UE_mmW*beta_uc(m,ue_idx));
            end
        end
    end
    D_mmW_mmW = zeros(N_UE_mmW,N_UE_mmW);
    for m = 1:M
        D_mmW_mmW = D_mmW_mmW + sqrt(eta_eq(m))*reshape(channel_dl_mmW(m,ue_idx,:,:),[Ntx,N_UE_mmW]).'*reshape(conj(channel_est_dl_mmW(m,ue_idx,:,:)),[Ntx,N_UE_mmW]);
    end        
    DS_mmW = zeros(N_UE_mmW,1);
    MSI_mmW = zeros(N_UE_mmW,1);
    
    noise_mmW = abs(sqrt(0.5)*(randn(N_UE_mmW,1) + 1j*randn(N_UE_mmW,1))).^2;
    snr_num_mmW = zeros(N_UE_mmW,1);
    snr_den_mmW = zeros(N_UE_mmW,1);
    rate_dl = 0;
    if (sub6ConnectionState(ue_idx)==1)
        for n = 1:N_UE_mmW
%             DS_mmW(k,n) = p_d*norm(reshape(D_mmW_mmW(k,k,n,:),[1,N_UE_mmW]))^2;
            DS_mmW(n) = p_d*(abs(D_mmW_mmW(n,n)))^2;
            for nn = 1:N_UE_mmW
%                 if (nn~=n)
                if (abs(D_mmW_mmW(nn,nn))<abs(D_mmW_mmW(n,n)))
%                     MSI_mmW(k,n) = MSI_mmW(k,n) + p_d*norm(reshape(D_mmW_mmW(k,k,nn,:),[1,N_UE_mmW]))^2;
                    MSI_mmW(n) = MSI_mmW(n) + p_d*(abs(D_mmW_mmW(n,nn)))^2;
                end
            end
            snr_num_mmW(n) = DS_mmW(n);
            snr_den_mmW(n) = MSI_mmW(n) + noise_mmW(n);
            rate_dl = rate_dl + BW*TAU_FAC*log2(1+snr_num_mmW(n)/snr_den_mmW(n));
        end
    end
end
end