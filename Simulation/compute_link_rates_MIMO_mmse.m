function rate_dl = compute_link_rates_MIMO_mmse(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,ue_idx,sub6ConnectionState)
M = size(channel_dl,1);
K_mmW = size(sub6ConnectionState,1);
K = K_mmW + size(channel_dl,2);
BW = params.Band;
TAU_FAC = params.preLogFactor;
Ntx = size(channel_dl,3);
N_UE_mmW = size(channel_dl_mmW,4);
N_UE_sub6 = size(channel_dl,4);
p_d = params.rho_tot; % 1*K;
p_fac = params.p_fac;
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
eta_eq = zeros(M,K);
N_AP = params.num_antennas_per_gNB;
if ((K_mmW == 0) || (sub6ConnectionState == zeros(K_mmW,1)))
    for m = 1:M
        for k = 1+K_mmW:K
            if ismember(m,Serv{k})
                eta_eq(m,k) = 1/(N_AP*N_UE_sub6*sum(beta_uc(m,:)));
            end
        end
    end
else
    for m = 1:M
        for k = 1:K
            if ismember(m,Serv{k})
                if ((k<=K_mmW) && (sub6ConnectionState(k) == 1))
                    eta_eq(m,k) = 1/(N_AP*(N_UE_mmW*(beta_uc(m,1:K_mmW).*sub6ConnectionState)+N_UE_sub6*sum(beta_uc(m,:))));
                elseif (k>K_mmW)
                    eta_eq(m,k) = 1/(N_AP*(N_UE_mmW*(beta_uc(m,1:K_mmW).*sub6ConnectionState)+N_UE_sub6*sum(beta_uc(m,:))));
                end
            end
        end
    end
end
D_mmW_mmW = zeros(K_mmW,K_mmW,N_UE_mmW,N_UE_mmW);
D_mmW_sub6 = zeros(K_mmW,K-K_mmW,N_UE_mmW,N_UE_sub6);
D_sub6_mmW = zeros(K-K_mmW,K_mmW,N_UE_sub6,N_UE_mmW);
D_sub6_sub6 = zeros(K-K_mmW,K-K_mmW,N_UE_sub6,N_UE_sub6);
dl_mmse_precoder_mmW = zeros(size(channel_est_dl_mmW));
dl_mmse_precoder = zeros(size(channel_est_dl));
for m = 1:M
    for k = 1:K_mmW
        inv_matrix = zeros(Ntx,N_UE_mmW);
        for q = 1:K_mmW
            if ismember(m,Serv{q})
                inv_matrix = inv_matrix + p_d*reshape(channel_dl_mmW(m,q,:,:),[Ntx,N_UE_mmW])*reshape(channel_dl_mmW(m,q,:,:),[Ntx,N_UE_mmW])';
            end
        end
        for q = 1:K-K_mmW
            if ismember(m,Serv{q+K_mmW})
                inv_matrix = inv_matrix + p_d*reshape(channel_dl(m,q,:,:),[Ntx,N_UE_sub6])*reshape(channel_dl(m,q,:,:),[Ntx,N_UE_sub6])';
            end
        end
        inv_matrix = inv_matrix + eye(Ntx);
        dl_mmse_precoder_mmW(m,k,:,:) = reshape(dl_mmse_precoder_mmW(m,k,:,:),[Ntx,N_UE_mmW]) + p_d*inv_matrix\(D(m,k)*reshape(channel_dl_mmW(m,k,:,:),[Ntx,N_UE_mmW]));
    end
    for k = 1:K-K_mmW
        inv_matrix = zeros(Ntx,N_UE_sub6);
        for q = 1:K_mmW
            if ismember(m,Serv{q})
                inv_matrix = inv_matrix + p_d*reshape(channel_dl_mmW(m,q,:,:),[Ntx,N_UE_mmW])*reshape(channel_dl_mmW(m,q,:,:),[Ntx,N_UE_mmW])';
            end
        end
        for q = 1:K-K_mmW
            if ismember(m,Serv{q+K_mmW})
                inv_matrix = inv_matrix +  p_d*reshape(channel_dl(m,q,:,:),[Ntx,N_UE_sub6])*reshape(channel_dl(m,q,:,:),[Ntx,N_UE_sub6])';
            end
        end
        inv_matrix = inv_matrix + eye(Ntx);
        dl_mmse_precoder(m,k,:,:) = reshape(dl_mmse_precoder(m,k,:,:),[Ntx,N_UE_sub6]) + p_d*inv_matrix\(D(m,k+K_mmW)*reshape(channel_dl(m,k,:,:),[Ntx,N_UE_sub6]));
    end
end
% dl_mmse_precoder_mmW = conj(channel_est_dl_mmW);
% dl_mmse_precoder = conj(channel_est_dl);
for k = 1:K_mmW
    for q = 1:K_mmW
        for m = 1:M
            D_mmW_mmW(k,q,:,:) = reshape(D_mmW_mmW(k,q,:,:),[N_UE_mmW,N_UE_mmW]) + sqrt(eta_eq(m,q))*reshape(channel_dl_mmW(m,k,:,:),[Ntx,N_UE_mmW]).'*reshape(dl_mmse_precoder_mmW(m,q,:,:),[Ntx,N_UE_mmW]);
        end
    end
    for q = 1:K-K_mmW
        for m = 1:M
            D_mmW_sub6(k,q,:,:) = reshape(D_mmW_sub6(k,q,:,:),[N_UE_mmW,N_UE_sub6]) + sqrt(eta_eq(m,q+K_mmW))*reshape(channel_dl_mmW(m,k,:,:),[Ntx,N_UE_mmW]).'*reshape(dl_mmse_precoder(m,q,:,:),[Ntx,N_UE_sub6]);
        end
    end
end
for k = 1:K-K_mmW
    for q = 1:K_mmW
        for m = 1:M
            D_sub6_mmW(k,q,:,:) = reshape(D_sub6_mmW(k,q,:,:),[N_UE_sub6,N_UE_mmW]) + sqrt(eta_eq(m,q))*reshape(channel_dl(m,k,:,:),[Ntx,N_UE_sub6]).'*reshape(dl_mmse_precoder_mmW(m,q,:,:),[Ntx,N_UE_mmW]);
        end
    end
    for q = 1:K-K_mmW
        for m = 1:M
            D_sub6_sub6(k,q,:,:) = reshape(D_sub6_sub6(k,q,:,:),[N_UE_sub6,N_UE_sub6]) + sqrt(eta_eq(m,q+K_mmW))*reshape(channel_dl(m,k,:,:),[Ntx,N_UE_sub6]).'*reshape(dl_mmse_precoder(m,q,:,:),[Ntx,N_UE_sub6]);
        end
    end
end
DS_mmW = zeros(K_mmW,N_UE_mmW);
MSI_mmW = zeros(K_mmW,N_UE_mmW);
MUI_mmW = zeros(K_mmW,N_UE_mmW);
DS_sub6 = zeros(K-K_mmW,N_UE_sub6);
MSI_sub6 = zeros(K-K_mmW,N_UE_sub6);
MUI_sub6 = zeros(K-K_mmW,N_UE_sub6);

noise_mmW = abs(sqrt(0.5)*(randn(K_mmW,N_UE_mmW) + 1j*randn(K_mmW,N_UE_mmW))).^2;
noise_sub6 = abs(sqrt(0.5)*(randn(K-K_mmW,N_UE_sub6) + 1j*randn(K-K_mmW,N_UE_sub6))).^2;
snr_num_mmW = zeros(K_mmW,N_UE_mmW);
snr_den_mmW = zeros(K_mmW,N_UE_mmW);
snr_num_sub6 = zeros(K-K_mmW,N_UE_sub6);
snr_den_sub6 = zeros(K-K_mmW,N_UE_sub6);
rate_dl = zeros(K,1);
for k = 1:K_mmW
%     if (sub6ConnectionState(k)==1 || k==ue_idx)
    if (sub6ConnectionState(k)==1)
        for n = 1:N_UE_mmW
%             DS_mmW(k,n) = p_d*norm(reshape(D_mmW_mmW(k,k,n,:),[1,N_UE_mmW]))^2;
            DS_mmW(k,n) = p_d*(abs(D_mmW_mmW(k,k,n,n)))^2;
            for nn = 1:N_UE_mmW
%                 if (nn~=n)
                if (abs(D_mmW_mmW(k,k,nn,nn))<abs(D_mmW_mmW(k,k,n,n)))
%                     MSI_mmW(k,n) = MSI_mmW(k,n) + p_d*norm(reshape(D_mmW_mmW(k,k,nn,:),[1,N_UE_mmW]))^2;
                    MSI_mmW(k,n) = MSI_mmW(k,n) + p_d*(abs(D_mmW_mmW(k,k,n,nn)))^2;
                end
            end
            for q = 1:K_mmW
                if (q~=k && sub6ConnectionState(q)==1)
                  MUI_mmW(k,n) = MUI_mmW(k,n) + p_d*norm(reshape(D_mmW_mmW(k,q,n,:),[1,N_UE_mmW]))^2;
                end
            end
            for q = 1:K-K_mmW
               MUI_mmW(k,n) = MUI_mmW(k,n) + p_d*norm(reshape(D_mmW_sub6(k,q,n,:),[1,N_UE_sub6]))^2;
            end
            snr_num_mmW(k,n) = DS_mmW(k,n);
            snr_den_mmW(k,n) = MUI_mmW(k,n) + noise_mmW(k,n);
            rate_dl(k) = rate_dl(k) + BW*TAU_FAC*log2(1+snr_num_mmW(k,n)/snr_den_mmW(k,n));
        end
    end
end
for k = 1:K-K_mmW
    for n = 1:N_UE_sub6
%         DS_sub6(k,n) = p_d*norm(reshape(D_sub6_sub6(k,k,n,:),[1,N_UE_sub6]))^2;
        DS_sub6(k,n) = p_d*(abs(D_sub6_sub6(k,k,n,n)))^2;
        for nn = 1:N_UE_sub6
%             if (nn~=n)
%             if (nn>n)
            if (abs(D_sub6_sub6(k,k,nn,nn))<abs(D_sub6_sub6(k,k,n,n)))
%                 MSI_sub6(k,n) = MSI_sub6(k,n) + p_d*norm(reshape(D_sub6_sub6(k,k,nn,:),[1,N_UE_sub6]))^2;
                MSI_sub6(k,n) = MSI_sub6(k,n) + p_d*(abs(D_sub6_sub6(k,k,n,nn)))^2;
            end
        end
        for q = 1:K_mmW
            if (q~=k && sub6ConnectionState(q)==1)
              MUI_sub6(k,n) = MUI_sub6(k,n) + p_d*norm(reshape(D_sub6_mmW(k,q,n,:),[1,N_UE_mmW]))^2;
            end
        end
        for q = 1:K-K_mmW
           MUI_sub6(k,n) = MUI_sub6(k,n) + p_d*norm(reshape(D_sub6_sub6(k,q,n,:),[1,N_UE_sub6]))^2;
        end
        snr_num_sub6(k,n) = DS_sub6(k,n);
        snr_den_sub6(k,n) = MUI_sub6(k,n) + noise_sub6(k,n);
        rate_dl(k+K_mmW) = rate_dl(k+K_mmW) + BW*TAU_FAC*log2(1+snr_num_sub6(k,n)/snr_den_sub6(k,n));
    end
end
end