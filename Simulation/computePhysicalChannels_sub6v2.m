% function [BETA, phy_channel_sub6, phy_channel_sub6_est] = computePhysicalChannels_sub6v2(params)
% function [BETA, phy_channel_sub6, phy_channel_sub6_est] = computePhysicalChannels_sub6v2(params,link)
% function [phy_channel_sub6, phy_channel_sub6_est] = computePhysicalChannels_sub6v2(params,link,BETA,ricianFactor)
function [phy_channel_sub6, phy_channel_sub6_est] = computePhysicalChannels_sub6v2(params,link)
% Kd = params.numUE;
Kd = params.numUE+params.numUE_sub6;
% tau = Kd;
% rho = 1;
Ntx = params.num_antennas_per_gNB;
% M = params.numGNB;
M = params.numGNB_sub6;
K = params.num_sc_sub6; % no of subcarriers
% scs = 15e3; %sub carrier spacing
scs = params.scs_sub6;
BW = K*scs;
BETA = params.BETA;
ricianFactor = params.ricianFactor;
% [BETA, ricianFactor] = channel_cellfree_2v2(params.numUE,params.numGNB,BW, params.locationsBS_sub6, params.UE_locations);
% [BETA, ricianFactor] = channel_cellfree_2v2(Kd,M,BW, params.locationsBS_sub6, [params.UE_locations;params.UE_locations_sub6]);
% [BETA, ricianFactor] = channel_cellfree_2v2(Kd,M,BW, params.locationsBS, [params.UE_locations;params.UE_locations_sub6]);
phy_channel_sub6 = zeros(M,Kd,Ntx);
% channel_bar_dl = zeros(M,Kd,Nrx);
% channel_ul = zeros(M,Ku,Nrx);
% channel_bar_ul = zeros(M,Ku,Ntx);
% c_dl = zeros(M,Kd);
% c_ul = zeros(M,Ku); 

for m = 1:M
    for k = 1:Kd
%         phy_channel_sub6 (m,k,:) = sqrt(0.5*BETA(m,k))*(randn(Ntx,1) + 1i*randn(Ntx,1));
        % channel_bar_dl (m,k,:) = sqrt(0.5*BETA(m,k))*(randn(Nrx,1) + 1i*randn(Nrx,1));
%         c_dl (m,k) = sqrt(rho*tau)*BETA(m,k)/(1+rho*tau*BETA(m,k));
%         for n = 1:Ntx
%             phy_channel_sub6_est (m,k,n) = c_dl(m,k)*(sqrt(rho*tau)*phy_channel_sub6 (m,k,n)+ W_tx(n,:)*PHI(:,k));
%         end      
        % phy_channel_sub6 (m,k,:) = sqrt(ricianFactor(m,k)/(1 + ricianFactor(m,k)))*sqrt(BETA(m,k)) + (sqrt(0.5*BETA(m,k))/sqrt(1 + ricianFactor(m,k)))*(randn(Ntx,1) + 1i*randn(Ntx,1));
        % if k<=params.numUE
        % isBSRecovered = link{k,m}.blockageStatus;
        % isDiscovered = link{k,m}.discovery_state;
    %         if isDiscovered
    %             if isBSRecovered
    %                 phy_channel_sub6 (m,k,:) = sqrt(ricianFactor(m,k)/(1 + ricianFactor(m,k)))*sqrt(0.5*BETA(m,k)) + (sqrt(0.5*BETA(m,k))/sqrt(1 + ricianFactor(m,k)))*(randn(Ntx,1) + 1i*randn(Ntx,1));
    %             else
    %                 phy_channel_sub6 (m,k,:) = (sqrt(0.5*BETA(m,k))/sqrt(1 + ricianFactor(m,k)))*(randn(Ntx,1) + 1i*randn(Ntx,1));
    %                 BETA(m,k) = BETA(m,k)/(1 + ricianFactor(m,k));
    %             end
    %         else
    %             BETA(m,k) = 0;
    %         end
        % if isBSRecovered && isDiscovered
        % if (k==2) %sub-6 link
            % phy_channel_sub6 (m,k,:) = sqrt(ricianFactor(m,k)/(1 + ricianFactor(m,k)))*sqrt(0.5*BETA(m,k)) + (sqrt(0.5*BETA(m,k))/sqrt(1 + ricianFactor(m,k)))*(randn(Ntx,1) + 1i*randn(Ntx,1));
            phy_channel_sub6 (m,k,:) = sqrt(0.5*BETA(m,k))*(randn(Ntx,1) + 1i*randn(Ntx,1));
        % elseif (k==1) %orig. mmW link got blocked now
            % phy_channel_sub6 (m,k,:) = (sqrt(0.5*BETA(m,k))/sqrt(1 + ricianFactor(m,k)))*(randn(Ntx,1) + 1i*randn(Ntx,1));
            % BETA(m,k) = BETA(m,k)/(1 + ricianFactor(m,k));
        % end
        % else
            % phy_channel_sub6 (m,k,:) = sqrt(ricianFactor(m,k)/(1 + ricianFactor(m,k)))*sqrt(BETA(m,k)) + (sqrt(0.5*BETA(m,k))/sqrt(1 + ricianFactor(m,k)))*(randn(Ntx,1) + 1i*randn(Ntx,1));
        % end
    end
    % for l = 1:Ku
    %     channel_ul (m,l,:) = sqrt(0.5*BETA(m,l+Kd))*(randn(Nrx,1)+1i*randn(Nrx,1));
    %     channel_bar_ul (m,l,:) = sqrt(0.5*BETA(m,l+Kd))*(randn(Ntx,1)+1i*randn(Ntx,1));
    %     c_ul (m,l) = sqrt(rho*tau)*BETA(m,l+Kd)/(1+rho*tau*BETA(m,l+Kd));
    %     for n = 1:Nrx
    %         channel_est_ul (m,l,n) = c_ul(m,l)*(sqrt(rho*tau)*channel_ul(m,l,n)+ W_rx(n,:)*PHI(:,l+Kd));
    %     end 
    % end
end
phy_channel_sub6_est = phy_channel_sub6;
% phy_channel_sub6_est = zeros(M,Kd,Ntx);
% channel_est_ul = zeros(M,Ku,Nrx);
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
%     % for l = 1:Ku
%     %     channel_ul (m,l,:) = sqrt(0.5*BETA(m,l+Kd))*(randn(Nrx,1)+1i*randn(Nrx,1));
%     %     channel_bar_ul (m,l,:) = sqrt(0.5*BETA(m,l+Kd))*(randn(Ntx,1)+1i*randn(Ntx,1));
%     %     c_ul (m,l) = sqrt(rho*tau)*BETA(m,l+Kd)/(1+rho*tau*BETA(m,l+Kd));
%     %     for n = 1:Nrx
%     %         channel_est_ul (m,l,n) = c_ul(m,l)*(sqrt(rho*tau)*channel_ul(m,l,n)+ W_rx(n,:)*PHI(:,l+Kd));
%     %     end 
%     % end
% end
end