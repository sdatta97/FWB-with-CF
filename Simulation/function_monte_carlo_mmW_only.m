% function [SE_k, SE_num_m, SE_den_m, HI_UE_rx_m,HI_AP_tr_m] = function_monte_carlo(L,K,K_mmW,N,eta,h,h_hat,K_AP_TR,K_UE_RX,no_of_rea,plos)
function [SE_k, SE_num_m, SE_den_m, HI_UE_rx_m,HI_AP_tr_m] = function_monte_carlo_mmW_only(L,K,K_mmW,N,N_mmW,N_UE,eta,h_mmW,h_hat_mmW,h,h_hat,K_AP_TR,K_UE_RX,no_of_rea,plos)

SE_k = zeros(K_mmW,1);
SE_num_m = zeros(K_mmW,1);
SE_den_m = zeros(K_mmW,1);
HI_UE_rx_m = zeros(K_mmW,1);
HI_AP_tr_m = zeros(K_mmW,1);

HI_UE_rx_monte_mmW =zeros(K_mmW,no_of_rea,L,K_mmW);
HI_AP_tx_monte_mmW =zeros(K_mmW,no_of_rea);
%% ---
for k = 1:K_mmW
    if (plos(k) > 0)
        SE_k_monte = zeros(no_of_rea,1);
        SE_num_m1 = zeros(no_of_rea,1);
        SE_den_m1 = zeros(no_of_rea,1);
        HI_UE_rx_m1 = zeros(no_of_rea,1);
        HI_AP_tx_m1 = zeros(no_of_rea,1);
        
        for ch = 1:no_of_rea    
            %numerator -- monte   h_hat(:,chreali,ap,k)
            num_monte = 0;
            for ap =1:L %a-th AP
                for n1 = 1:N_UE
                    for n2 = 1:N_UE
                        num_monte = num_monte + sqrt(eta(ap,k))*h_mmW(:,n1,ch,ap,k)'*h_hat_mmW(:,n2,ch,ap,k);
                    end
                end
                HI_AP_variance_mmW = zeros(N_mmW,N_mmW);
                warning('off','all');
                for kk = 1:K_mmW
                %----------MONTE-CARLO
                    HI_UE_rx_monte_mmW(k,ch,ap,kk) = K_UE_RX*eta(ap,kk)*sum(sum((h_mmW(:,:,ch,ap,k)' *(h_hat_mmW(:,:,ch,ap,kk)*h_hat_mmW(:,:,ch,ap,kk)'*h_mmW(:,:,ch,ap,k))),2),1) ...
                        + K_AP_TR*eta(ap,kk)*sum(sum((h_mmW(:,:,ch,ap,k)'*diag(diag(h_hat_mmW(:,:,ch,ap,kk)*h_hat_mmW(:,:,ch,ap,kk)'))*h_mmW(:,:,ch,ap,k)),2),1);   
                end
                HI_AP_variance_mmW = HI_AP_variance_mmW + K_AP_TR*eta(ap,kk)*diag(diag(h_hat_mmW(:,:,ch,ap,kk)*h_hat_mmW(:,:,ch,ap,kk)')); % sum over UEs -- a-th AP HW 
%                 HI_AP_tx_monte(k,ch,ap,kk) =K_AP_TR*eta(ap,kk)*(h(:,ch,ap,k)'*diag(diag(h_hat(:,ch,ap,kk)*h_hat(:,ch,ap,kk)'))*h(:,ch,ap,k)); %closed-frominsome
                HI_AP_tx_monte_mmW(k,ch) = HI_AP_tx_monte_mmW(k,ch) + sum(h_mmW(:,:,ch,ap,k)'* sqrtm(HI_AP_variance_mmW)*sqrt(0.5)*(randn(N_mmW,1)+1i*randn(N_mmW,1))); %sum over APs
            end
            HI_UE_RX_MC1 = sum(sum(HI_UE_rx_monte_mmW,4),3); 
            HI_UE_rx_m1(ch) = abs( sqrt(0.5*HI_UE_RX_MC1(k,ch))*(randn(1,1)+1i*randn(1,1)))^2;
            HI_AP_tx_m1(ch) = abs( HI_AP_tx_monte_mmW(k,ch))^2;
                    
            snr_num_monte_c = abs(num_monte)^2;
            % num_montec1=num_montec1+num_monte;
            
            interference_kj_monte_sum=0;
            for kd = 1:K_mmW
                if kd~=k
                    interference_kj_monte =0;
                    for ap2 = 1:L
                        % interference_kj_monte = interference_kj_monte + sqrt(eta(ap2,kd))*h(:,ch,ap2,k)'*h_hat(:,ch,ap2,kd);
                        for n1 = 1:N_UE
                            for n2 = 1:N_UE
                                interference_kj_monte = interference_kj_monte + sqrt(eta(ap2,kd))*h_mmW(:,n1,ch,ap2,k)'*h_hat_mmW(:,n2,ch,ap2,kd);
                            end
                        end
                    end
                    interference_kj_monte = abs(interference_kj_monte)^2;
                    interference_kj_monte_sum = interference_kj_monte_sum + plos(kd)*interference_kj_monte;
                end
            end
            for kd = 1+K_mmW:K
                interference_kj_monte =0;
                for ap2 = 1:L
                    % interference_kj_monte = interference_kj_monte + sqrt(eta(ap2,kd))*h(:,ch,ap2,k)'*h_hat(:,ch,ap2,kd);
                    for n = 1:N_UE
                        interference_kj_monte = interference_kj_monte + sqrt(eta(ap2,kd))*h_mmW(:,n,ch,ap2,k)'*h_hat(:,ch,ap2,kd-K_mmW);
                    end
                end
                interference_kj_monte = abs(interference_kj_monte)^2;
                interference_kj_monte_sum = interference_kj_monte_sum + interference_kj_monte;
            end
            %SE for each realization
            SE_k_monte(ch) = log2(1+ snr_num_monte_c/(1+HI_AP_tx_m1(ch)+ HI_UE_rx_m1(ch) + interference_kj_monte_sum ));        
            %to check values
            SE_num_m1(ch) = snr_num_monte_c;
            SE_den_m1(ch) = HI_AP_tx_m1(ch)+ HI_UE_rx_m1(ch) + interference_kj_monte_sum;
            
        end
        
        SE_num_m(k)= mean(SE_num_m1);
        SE_den_m(k) = mean(SE_den_m1);
        HI_UE_rx_m(k) = mean(HI_UE_rx_m1);
        HI_AP_tr_m(k) = mean(HI_AP_tx_m1);    
        SE_k(k)= mean(SE_k_monte);   % SE
    end
end
end