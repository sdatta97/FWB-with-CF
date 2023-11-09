function [phy_channel_sub6, phy_channel_sub6_est, phy_channel_mmW, phy_channel_mmW_est] = computePhysicalChannels_sub6_MIMO(params)
K = params.numUE+params.numUE_sub6;
K_mmW = params.numUE;
Ntx = params.num_antennas_per_gNB;
N_mmW = params.N_UE_mmW;
N_sub6 = params.N_UE_sub6;
% M = params.numGNB;
M = params.numGNB_sub6;
BETA = params.BETA;
phy_channel_mmW = zeros(M,K_mmW,Ntx,N_mmW);
phy_channel_sub6 = zeros(M,K-K_mmW,Ntx,N_sub6);
for m = 1:M
    for k = 1:K_mmW
        phy_channel_mmW (m,k,:,:) = sqrt(0.5*BETA(m,k))*(randn(Ntx,N_mmW) + 1i*randn(Ntx,N_mmW));        
    end
    for k = 1:K-K_mmW
        phy_channel_sub6 (m,k,:,:) = sqrt(0.5*BETA(m,k+K_mmW))*(randn(Ntx,N_sub6) + 1i*randn(Ntx,N_sub6));        
    end 
end
phy_channel_mmW_est = phy_channel_mmW;
phy_channel_sub6_est = phy_channel_sub6;
% phy_channel_sub6_est = zeros(M,Kd,Ntx);
% PHI1    = orth(rand(tau));   % generate an orthonormal matrix of dimension tau_p
% PHI     = zeros(size(PHI1));
% perm_vec  = repmat(randperm(tau),1,2);
% phi_index = perm_vec(1:K);
% for k = 1:K
%     PHI(:,k) = PHI1(:,phi_index(k));
% end
% for m = 1:M
%     W_tx = sqrt(0.5)*(randn(Ntx, tau)+1i*randn(Ntx,tau));
%     % W_rx = sqrt(0.5)*(randn(Nrx, tau)+1i*randn(Nrx,tau));
%     for k = 1:Kd
%         phy_channel_sub6 (m,k,:) = sqrt(0.5*BETA(m,k))*(randn(Ntx,1) + 1i*randn(Ntx,1));
%         % channel_bar_dl (m,k,:) = sqrt(0.5*BETA(m,k))*(randn(Nrx,1) + 1i*randn(Nrx,1));
%         c_dl (m,k) = sqrt(rho*tau)*BETA(m,k)/(1+rho*tau*BETA(m,k));
%         for n = 1:Ntx
%             phy_channel_sub6_est (m,k,n) = c_dl(m,k)*(sqrt(rho*tau)*phy_channel_sub6 (m,k,n)+ W_tx(n,:)*PHI(:,k));
%         end        
%     end
% end
end