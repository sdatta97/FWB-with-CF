% function [SE_k, SE_num_m, SE_den_m, HI_UE_rx_m,HI_AP_tr_m] = function_monte_carlo(L,K,K_mmW,N,eta,h,h_hat,K_AP_TR,K_UE_RX,no_of_rea,plos)
% function [SE_k, SE_num_m, SE_den_m, HI_UE_rx_m,HI_AP_tr_m] = function_monte_carlo(L,K,K_mmW,N,N_mmW,eta,h_mmW,h_hat_mmW,h,h_hat,K_AP_TR,K_UE_RX,no_of_rea,plos)
function [SE_k, SE_num_m, SE_den_m, HI_UE_rx_m,HI_AP_tr_m] = function_monte_carlo(L,K,K_mmW,N,N_mmW,eta,h_mmW,h_hat_mmW,h,h_hat,K_AP_TR,K_UE_RX,no_of_rea,plos)

% SE_k = zeros(K,1);
% SE_num_m = zeros(K,1);
% SE_den_m = zeros(K,1);
% HI_UE_rx_m = zeros(K,1);
% HI_AP_tr_m = zeros(K,1);
SE_k = zeros(K-K_mmW,1);
SE_num_m = zeros(K-K_mmW,1);
SE_den_m = zeros(K-K_mmW,1);
HI_UE_rx_m = zeros(K-K_mmW,1);
HI_AP_tr_m = zeros(K-K_mmW,1);

% HI_UE_rx_monte_mmW =zeros(K_mmW,no_of_rea,L,K_mmW);
% HI_AP_tx_monte_mmW =zeros(K_mmW,no_of_rea);
HI_UE_rx_monte_sub6 =zeros(K-K_mmW,no_of_rea,L,K-K_mmW);
HI_AP_tx_monte_sub6 =zeros(K-K_mmW,no_of_rea);
%% ---
% for k = 1:K
for k = 1+K_mmW:K
    % if ((k > K_mmW) || (k<=K_mmW && plos(k) > 0))
    SE_k_monte = zeros(no_of_rea,1);
    SE_num_m1 = zeros(no_of_rea,1);
    SE_den_m1 = zeros(no_of_rea,1);
    HI_UE_rx_m1 = zeros(no_of_rea,1);
    HI_AP_tx_m1 = zeros(no_of_rea,1);
    
    for ch = 1:no_of_rea    
        %numerator -- monte   h_hat(:,chreali,ap,k)
        num_monte = 0;
        for ap =1:L %a-th AP
            % if k<=K_mmW
            %     num_monte = num_monte + sqrt(eta(ap,k))*h_mmW(:,ch,ap,k)'*h_hat_mmW(:,ch,ap,k);
            % else
            num_monte = num_monte + sqrt(eta(ap,k))*h(:,ch,ap,k-K_mmW)'*h_hat(:,ch,ap,k-K_mmW);
            % end
            % HI_AP_variance_mmW = zeros(N_mmW,N_mmW);
            HI_AP_variance_sub6 = zeros(N,N);
            warning('off','all');
%             if (k<=K_mmW)
%                 for kk = 1:K_mmW
%                 %----------MONTE-CARLO
%                     HI_UE_rx_monte_mmW(k,ch,ap,kk) = K_UE_RX*(eta(ap,kk)*(h_mmW(:,ch,ap,k)' *(h_hat_mmW(:,ch,ap,kk)*h_hat_mmW(:,ch,ap,kk)'*h_mmW(:,ch,ap,k))) ...
%                         + K_AP_TR*eta(ap,kk)*(h_mmW(:,ch,ap,k)'*diag(diag(h_hat_mmW(:,ch,ap,kk)*h_hat_mmW(:,ch,ap,kk)'))*h_mmW(:,ch,ap,k)));   
%                 end
%                 HI_AP_variance_mmW = HI_AP_variance_mmW + K_AP_TR*eta(ap,kk)*diag(diag(h_hat_mmW(:,ch,ap,kk)*h_hat_mmW(:,ch,ap,kk)')); % sum over UEs -- a-th AP HW 
% %                 HI_AP_tx_monte(k,ch,ap,kk) =K_AP_TR*eta(ap,kk)*(h(:,ch,ap,k)'*diag(diag(h_hat(:,ch,ap,kk)*h_hat(:,ch,ap,kk)'))*h(:,ch,ap,k)); %closed-frominsome
%                 HI_AP_tx_monte_mmW(k,ch) = HI_AP_tx_monte_mmW(k,ch) + (h_mmW(:,ch,ap,k)'* sqrtm(HI_AP_variance_mmW)*sqrt(0.5)*(randn(N_mmW,1)+1i*randn(N_mmW,1)) ); %sum over APs
%             else
            for kk = 1:K-K_mmW
            %----------MONTE-CARLO
                HI_UE_rx_monte_sub6(k-K_mmW,ch,ap,kk) = K_UE_RX*(eta(ap,kk)*(h(:,ch,ap,k-K_mmW)' *(h_hat(:,ch,ap,kk)*h_hat(:,ch,ap,kk)'*h(:,ch,ap,k-K_mmW))) ...
            + K_AP_TR*eta(ap,kk)*(h(:,ch,ap,k-K_mmW)'*diag(diag(h_hat(:,ch,ap,kk)*h_hat(:,ch,ap,kk)'))*h(:,ch,ap,k-K_mmW)));                    
%                 HI_AP_tx_monte(k,ch,ap,kk) =K_AP_TR*eta(ap,kk)*(h(:,ch,ap,k)'*diag(diag(h_hat(:,ch,ap,kk)*h_hat(:,ch,ap,kk)'))*h(:,ch,ap,k)); %closed-frominsome
            end
            HI_AP_variance_sub6 = HI_AP_variance_sub6 + K_AP_TR*eta(ap,kk)*diag(diag(h_hat(:,ch,ap,kk)*h_hat(:,ch,ap,kk)')); % sum over UEs -- a-th AP HW 
            HI_AP_tx_monte_sub6(k-K_mmW,ch) = HI_AP_tx_monte_sub6(k-K_mmW,ch) + (h(:,ch,ap,k-K_mmW)'* sqrtm(HI_AP_variance_sub6)*sqrt(0.5)*(randn(N,1)+1i*randn(N,1)) ); %sum over APs
            % end
        end
        % if (k<=K_mmW)
        %     HI_UE_RX_MC1 = sum(sum(HI_UE_rx_monte_mmW,4),3); 
        %     HI_UE_rx_m1(ch) = abs( sqrt(0.5*HI_UE_RX_MC1(k,ch))*(randn(1,1)+1i*randn(1,1)))^2;
        %     HI_AP_tx_m1(ch) = abs( HI_AP_tx_monte_mmW(k,ch))^2;
        % else
        HI_UE_RX_MC1 = sum(sum(HI_UE_rx_monte_sub6,4),3);
        HI_UE_rx_m1(ch) = abs( sqrt(0.5*HI_UE_RX_MC1(k-K_mmW,ch))*(randn(1,1)+1i*randn(1,1)))^2;
        HI_AP_tx_m1(ch) = abs( HI_AP_tx_monte_sub6(k-K_mmW,ch))^2;
        % end
    
        snr_num_monte_c = abs(num_monte)^2;
        % num_montec1=num_montec1+num_monte;
        
        interference_kj_monte_sum=0;
        % if (k<=K_mmW)
        %     for kd = 1:K_mmW
        %         if kd~=k
        %             interference_kj_monte =0;
        %             for ap2 = 1:L
        %                 % interference_kj_monte = interference_kj_monte + sqrt(eta(ap2,kd))*h(:,ch,ap2,k)'*h_hat(:,ch,ap2,kd);
        %                 interference_kj_monte = interference_kj_monte + sqrt(eta(ap2,kd))*h_mmW(:,ch,ap2,k)'*h_hat_mmW(:,ch,ap2,kd);
        %             end
        %             interference_kj_monte = abs(interference_kj_monte)^2;
        %             interference_kj_monte_sum = interference_kj_monte_sum + plos(kd)*interference_kj_monte;
        %         end
        %     end
        % else
        for kd = 1:K_mmW
            interference_kj_monte =0;
            for ap2 = 1:L
                N_UE = size(h_hat_mmW,2);
                interference_kj_monte = interference_kj_monte + sqrt(eta(ap2,kd))*h(:,ch,ap2,k-K_mmW)'*(h_hat_mmW(:,:,ch,ap2,kd)*ones(N_UE,1));
            end
            interference_kj_monte = abs(interference_kj_monte)^2;
            interference_kj_monte_sum = interference_kj_monte_sum + interference_kj_monte;
        end
        for kd = 1+K_mmW:K
            if kd~=k
                interference_kj_monte =0;
                for ap2 = 1:L
                    % interference_kj_monte = interference_kj_monte + sqrt(eta(ap2,kd))*h(:,ch,ap2,k)'*h_hat(:,ch,ap2,kd);
                    interference_kj_monte = interference_kj_monte + sqrt(eta(ap2,kd))*h(:,ch,ap2,k-K_mmW)'*h_hat(:,ch,ap2,kd-K_mmW);
                end
                interference_kj_monte = abs(interference_kj_monte)^2;
                interference_kj_monte_sum = interference_kj_monte_sum + interference_kj_monte;
            end
        end
        % end
        %SE for each realization
        SE_k_monte(ch) = log2(1+ snr_num_monte_c/(1+HI_AP_tx_m1(ch)+ HI_UE_rx_m1(ch) + interference_kj_monte_sum ));        
        %to check values
        SE_num_m1(ch) = snr_num_monte_c;
        SE_den_m1(ch) = HI_AP_tx_m1(ch)+ HI_UE_rx_m1(ch) + interference_kj_monte_sum;
        
    end
    
    % SE_num_m(k)= mean(SE_num_m1);
    % SE_den_m(k) = mean(SE_den_m1);
    % HI_UE_rx_m(k) = mean(HI_UE_rx_m1);
    % HI_AP_tr_m(k) = mean(HI_AP_tx_m1);    
    % SE_k(k)= mean(SE_k_monte);   % SE
    SE_num_m(k-K_mmW)= mean(SE_num_m1);
    SE_den_m(k-K_mmW) = mean(SE_den_m1);
    HI_UE_rx_m(k-K_mmW) = mean(HI_UE_rx_m1);
    HI_AP_tr_m(k-K_mmW) = mean(HI_AP_tx_m1);    
    SE_k(k-K_mmW)= mean(SE_k_monte);   % SE
end
end