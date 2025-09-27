function [phy_channel_sub6, phy_channel_sub6_est, phy_channel_mmW, phy_channel_mmW_est] = computePhysicalChannels_sub6_MIMO(params)
K = params.numUE+params.numUE_sub6;
K_mmW = params.numUE;
Ntx = params.num_antennas_per_gNB;
N_mmW = params.N_UE_mmW;
N_sub6 = params.N_UE_sub6;
M = params.numGNB;
BETA = params.BETA;
D = params.D;
R_gNB = params.R_gNB;
R_ue_mmW = params.R_ue_mmW;
R_ue_sub6 = params.R_ue_sub6;
phy_channel_mmW = zeros(M,K_mmW,Ntx,N_mmW);
phy_channel_sub6 = zeros(M,K-K_mmW,Ntx,N_sub6);
for m = 1:M
    for k = 1:K_mmW
%         phy_channel_mmW (m,k,:,:) = sqrt(0.5*BETA(m,k))*(randn(Ntx,N_mmW) + 1i*randn(Ntx,N_mmW));        
        phy_channel_mmW (m,k,:,:) = sqrt(0.5)*sqrtm(R_gNB(:,:,m,k,1))*(randn(Ntx,N_mmW) + 1i*randn(Ntx,N_mmW))*sqrtm(R_ue_mmW(:,:,m,k,1));        
    end
    for k = 1:K-K_mmW
%         phy_channel_sub6 (m,k,:,:) = sqrt(0.5*BETA(m,k+K_mmW))*(randn(Ntx,N_sub6) + 1i*randn(Ntx,N_sub6));        
        phy_channel_sub6 (m,k,:,:) = sqrt(0.5)*sqrtm(R_gNB(:,:,m,k+K_mmW,1))*(randn(Ntx,N_sub6) + 1i*randn(Ntx,N_sub6))*sqrtm(R_ue_sub6(:,:,m,k,1));        
    end 
end
if params.Perf_CSI
    phy_channel_mmW_est = phy_channel_mmW;
    phy_channel_sub6_est = phy_channel_sub6;
else
    phy_channel_mmW_est = zeros(M,K_mmW,Ntx,N_mmW);
    phy_channel_sub6_est = zeros(M,K-K_mmW,Ntx,N_sub6);
    tau = params.tau_p;
    rho = params.pilot_pow;
    PHI1    = orth(rand(tau));   % generate an orthonormal matrix of dimension tau_p
    PHI     = zeros(size(PHI1));
    % perm_vec  = repmat(randperm(tau),1,2);
    % phi_index = perm_vec(1:K);
    for k = 1:K_mmW
        PHI(:,k) = PHI1(:,k);
    end
    for k = (1+K_mmW):K
        % PHI(:,k) = PHI1(:,randi(K_mmW));
        [~,idx] = min(sum(BETA(find(D(:,k)),1:K_mmW),1));
        PHI(:,k) = PHI1(:,idx);
    end
    for m = 1:M
        W_tx_mmW = sqrt(0.5)*(randn(Ntx, tau, N_mmW)+1i*randn(Ntx, tau, N_mmW));
        W_tx_sub6 = sqrt(0.5)*(randn(Ntx, tau, N_sub6)+1i*randn(Ntx, tau, N_sub6));
        for k = 1:K
            c_dl = sqrt(rho*tau)*BETA(m,k)/(1+rho*tau*BETA(m,k));
            for n = 1:Ntx
                if (k<=K_mmW)
                    phy_channel_mmW_est (m,k,n,:) = c_dl*(sqrt(rho*tau)*phy_channel_mmW(m,k,n,:)+ reshape(permute(W_tx_mmW(n,:,:),[1,3,2]),[N_mmW,tau])*PHI(:,k));
                else
                    phy_channel_sub6_est (m,k-K_mmW,n,:) = c_dl*(sqrt(rho*tau)*phy_channel_sub6(m,k-K_mmW,n,:)+ reshape(permute(W_tx_sub6(n,:,:),[1 3 2]),[N_sub6,tau])*PHI(:,k));
                end
            end        
        end
    end
end
end