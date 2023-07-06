function [SE_k7, snr_numerator1, snr_denominator,HI_term_UE_rx1, HI_term_AP_tx1,beamforming,INTERFERENCE_UAV_GUE ] = function_LB_impCSI(K_mmW,K,L,N,eta,h_LOS,R,psi,eta_p,PHI,K_AP_TX,K_UE_RX,gamma, gamma_MAT, beta_actual, beta_actual_MAT,C_ERR,GAMMA_NLOS, plos, plos2)

%function [SE_k,snr_numerator1,snr_denominator,Interfsum,beamforming,HI_term_AP_tx1,HI_term_UE_rx1] = function_LB_impCSI(K,L,N,eta,h_LOS,R,psi,eta_p,PHI,k_t2,k_r2_UE_rx,gamma )
%NUMERATOR

%% CHECK UE RX imapriment %% Not multiplied by HW factors, and TRANSMIT POWERS
INTERFERENCE_UAV_GUE = zeros(K,K);
HI_UE_rx_CL =zeros(K,K,L);
HI_AP_tx_CL =zeros(K,K,L);
% for k = 1:K  %k-th USER
%     for ap =1:L %a-th AP
%         for kk = 1:K
%             %----------COSED-FORM  K_AP_TX,K_UE_RX
%             if kk~= k
%                 HI_UE_rx_CL(k,kk,ap) = K_UE_RX*eta(ap,kk)*(trace(beta_actual_MAT(:,:,ap,k)*gamma_MAT(:,:,ap,kk))+ K_AP_TX*trace(beta_actual_MAT(:,:,ap,k)*diag(diag(gamma_MAT(:,:,ap,kk))))); %correct %WITHOUT TR.ap hw
%                 HI_AP_tx_CL(k,kk,ap) = K_AP_TX*eta(ap,kk)*trace(beta_actual_MAT(:,:,ap,k)*diag(diag(gamma_MAT(:,:,ap,kk))));
%             else
%                 HI_UE_rx_CL(k,kk,ap) = K_UE_RX*(eta(ap,kk)*(abs(trace(GAMMA_NLOS(:,:,ap,kk)))^2 + trace(GAMMA_NLOS(:,:,ap,kk)^2) + abs(h_LOS(:,ap,kk)'*h_LOS(:,ap,kk))^2 ...
%                     + 2*h_LOS(:,ap,kk)'*GAMMA_NLOS(:,:,ap,kk)*h_LOS(:,ap,kk) + 2*real(trace(GAMMA_NLOS(:,:,ap,kk))*h_LOS(:,ap,kk)'*h_LOS(:,ap,kk)) +trace(C_ERR(:,:,ap,kk)*gamma_MAT(:,:,ap,kk)))...
%                     + K_AP_TX*eta(ap,kk)*((2/N)*trace(GAMMA_NLOS(:,:,ap,kk))^2 + (4/N)*trace(GAMMA_NLOS(:,:,ap,kk))*h_LOS(:,ap,kk)'*h_LOS(:,ap,kk)...
%                     + sum(abs(abs(h_LOS(:,ap,kk)).^2).^2) +trace(C_ERR(:,:,ap,kk)*diag(diag(gamma_MAT(:,:,ap,kk))))));
%                 HI_AP_tx_CL(k,kk,ap) = K_AP_TX*eta(ap,kk)*((2/N)*trace(GAMMA_NLOS(:,:,ap,kk))^2 + (4/N)*trace(GAMMA_NLOS(:,:,ap,kk))*h_LOS(:,ap,kk)'*h_LOS(:,ap,kk)...
%                     + sum(abs(abs(h_LOS(:,ap,kk)).^2).^2) +trace(C_ERR(:,:,ap,kk)*diag(diag(gamma_MAT(:,:,ap,kk)))));
%             end
%         end
%     end
% end


SE_k7 = zeros(K,1);
snr_numerator1 = zeros(K,1);
snr_denominator= zeros(K,1);
Interfsum = zeros(K,1);
beamforming = zeros(K,1);
HI_term_UE_rx1 = zeros(K,1);
HI_term_AP_tx1 = zeros(K,1);
for k=1:K
    snr_numerator = 0;
    HI_term_AP_tr =0;
    HI_UE_rx_temp = 0;
    for ap =1:L
        snr_numerator = snr_numerator + sqrt(eta(ap,k))*gamma(ap,k); %norm(h_LOS(:,ap,k))^2 + eta_p*trace(R(:,:,ap,k)*psi(:,:,ap,k)*R(:,:,ap,k)) );
        for kk = 1:K
            %----------COSED-FORM  K_AP_TX,K_UE_RX
            if  round( PHI(:,kk)'* PHI(:,k)) ==1   % same pilots j in P_k \k %kk~= k
                HI_UE_rx_CL(k,kk,ap) = K_UE_RX*(eta(ap,kk)*(abs(trace(sqrtm(GAMMA_NLOS(:,:,ap,kk))*sqrtm(GAMMA_NLOS(:,:,ap,k)) ))^2 + trace(GAMMA_NLOS(:,:,ap,kk)*GAMMA_NLOS(:,:,ap,k)) + abs(h_LOS(:,ap,kk)'*h_LOS(:,ap,k))^2 ...
                    + h_LOS(:,ap,k)'*GAMMA_NLOS(:,:,ap,kk)*h_LOS(:,ap,k) + h_LOS(:,ap,kk)'*GAMMA_NLOS(:,:,ap,k)*h_LOS(:,ap,kk) + 2*real(trace(sqrtm(GAMMA_NLOS(:,:,ap,kk))*sqrtm(GAMMA_NLOS(:,:,ap,k)))*h_LOS(:,ap,kk)'*h_LOS(:,ap,k)) +trace(C_ERR(:,:,ap,k)*gamma_MAT(:,:,ap,kk)))...
                    + K_AP_TX*eta(ap,kk)*((2/N)*trace(GAMMA_NLOS(:,:,ap,kk))*trace(GAMMA_NLOS(:,:,ap,k)) + (2/N)*sqrt(trace(GAMMA_NLOS(:,:,ap,kk)))*sqrt(trace(GAMMA_NLOS(:,:,ap,k)))*h_LOS(:,ap,kk)'*h_LOS(:,ap,k)...
                    + trace(GAMMA_NLOS(:,:,ap,k))*h_LOS(:,ap,kk)'*h_LOS(:,ap,kk) + trace(GAMMA_NLOS(:,:,ap,kk))*h_LOS(:,ap,k)'*h_LOS(:,ap,k) + sum(abs((h_LOS(:,ap,kk).*h_LOS(:,ap,k)).^2)) +trace(C_ERR(:,:,ap,k)*diag(diag(gamma_MAT(:,:,ap,kk)))) ));               
                HI_AP_tx_CL(k,kk,ap) = K_AP_TX*eta(ap,kk)*((2/N)*trace(GAMMA_NLOS(:,:,ap,kk))*trace(GAMMA_NLOS(:,:,ap,k)) + (2/N)*sqrt(trace(GAMMA_NLOS(:,:,ap,kk)))*sqrt(trace(GAMMA_NLOS(:,:,ap,k)))*h_LOS(:,ap,kk)'*h_LOS(:,ap,k)...
                     + trace(GAMMA_NLOS(:,:,ap,k))*h_LOS(:,ap,kk)'*h_LOS(:,ap,kk) + trace(GAMMA_NLOS(:,:,ap,kk))*h_LOS(:,ap,k)'*h_LOS(:,ap,k) + sum(abs((h_LOS(:,ap,kk).*h_LOS(:,ap,k)).^2)) +trace(C_ERR(:,:,ap,k)*diag(diag(gamma_MAT(:,:,ap,kk)))));
               
            else
                HI_UE_rx_CL(k,kk,ap) = K_UE_RX*eta(ap,kk)*(trace(beta_actual_MAT(:,:,ap,k)*gamma_MAT(:,:,ap,kk))+ K_AP_TX*trace(beta_actual_MAT(:,:,ap,k)*diag(diag(gamma_MAT(:,:,ap,kk))))); %correct %WITHOUT TR.ap hw
                HI_AP_tx_CL(k,kk,ap) = K_AP_TX*eta(ap,kk)*trace(beta_actual_MAT(:,:,ap,k)*diag(diag(gamma_MAT(:,:,ap,kk))));
            end
        end
    end
    %%-------------------------  
    HI_UE_rx_CL3 = sum(HI_UE_rx_CL(k,:,:),3);
    HI_UE_rx_CL2 = sum(HI_UE_rx_CL3);

    HI_AP_tx_CL3 = sum(HI_AP_tx_CL(k,:,:),3);
    HI_AP_tx_CL2 = sum(HI_AP_tx_CL3);
    
    snr_numerator = abs(snr_numerator)^2;
    HI_term_AP_tx1(k) = HI_AP_tx_CL2;
    HI_term_UE_rx1(k) = HI_UE_rx_CL2;
    
    %DENOMINATOR
    % psi(N,N,L,K)   (N,N,L,K)   R(N,N,L,K)  gamma(L,K)  ,h_LOS(N,L,K)
    interference_sum =0;
    if (k <= K_mmW)
        for kd = 1:K_mmW
            if (kd~=k)
                interference_temp1 =0;
                interference_temp2 =0;
                for ap =1:L
                    if round( PHI(:,kd)'* PHI(:,k)) ==1   % same pilots j in P_k \k
                        interference_temp1 = interference_temp1 +  eta(ap,kd)*( h_LOS(:,ap,k)'*h_LOS(:,ap,kd)*h_LOS(:,ap,kd)'*h_LOS(:,ap,k)...
                            + eta_p*h_LOS(:,ap,k)'*h_LOS(:,ap,kd)*trace(psi(:,:,ap,k)*R(:,:,ap,kd)*R(:,:,ap,k)) +eta_p*h_LOS(:,ap,k)'*R(:,:,ap,kd)*psi(:,:,ap,k)*R(:,:,ap,kd)*h_LOS(:,ap,k)...
                            + eta_p*trace(psi(:,:,ap,k)*R(:,:,ap,k)*h_LOS(:,ap,kd)*h_LOS(:,ap,kd)'*R(:,:,ap,k)) +eta_p*trace( psi(:,:,ap,k)*R(:,:,ap,k)*R(:,:,ap,kd))*h_LOS(:,ap,kd)'*h_LOS(:,ap,k)...
                            + (eta_p^2)*trace( psi(:,:,ap,k)*R(:,:,ap,k)*R(:,:,ap,kd))^2 + (eta_p^2)*trace( psi(:,:,ap,k)*R(:,:,ap,k)*R(:,:,ap,kd)*psi(:,:,ap,k)*R(:,:,ap,kd)*R(:,:,ap,k))...
                            + trace( (R(:,:,ap,k)-eta_p*R(:,:,ap,k)*psi(:,:,ap,k)*R(:,:,ap,k))*(h_LOS(:,ap,kd)* h_LOS(:,ap,kd)' + eta_p*R(:,:,ap,kd)*psi(:,:,ap,k)*R(:,:,ap,kd)) ) ) ; % p + q terms est, error
                        
                    else
                        interference_temp1 = interference_temp1 + plos(kd)*eta(ap,kd)*trace((h_LOS(:,ap,k)*h_LOS(:,ap,k)'+R(:,:,ap,k))*(h_LOS(:,ap,kd)*h_LOS(:,ap,kd)' + eta_p*R(:,:,ap,kd)*psi(:,:,ap,kd)*R(:,:,ap,kd)));
                    end
                    %interference_temp1 = interference_temp1 +  eta(ap,kd)*( (norm(h_LOS(:,ap,k))^2 *norm(h_LOS(:,ap,kd))^2) + (norm(h_LOS(:,ap,k))^2 * trace(R(:,:,ap,kd)))...
                    %   +(norm(h_LOS(:,ap,kd))^2 * trace(R(:,:,ap,k))) + trace(R(:,:,ap,k))*trace(R(:,:,ap,kd))  );
                    for ap2 = 1:L
                        if ap2 ~= ap
                            if round( PHI(:,kd)'* PHI(:,k)) ==1   % same pilots j in P_k \k
                                interference_temp2 = interference_temp2 + sqrt(eta(ap,kd)*eta(ap2,kd))* (trace(h_LOS(:,ap,kd)*h_LOS(:,ap,k)' + eta_p*R(:,:,ap,kd)*psi(:,:,ap,kd)*R(:,:,ap,k))...
                                    *trace( h_LOS(:,ap2,k)*h_LOS(:,ap2,kd)' + eta_p*R(:,:,ap2,k)*psi(:,:,ap2,kd)*R(:,:,ap2,kd)) );
                            else
                                interference_temp2 = interference_temp2 + plos(kd)*sqrt(eta(ap,kd)*eta(ap2,kd))* ( h_LOS(:,ap,k)' *h_LOS(:,ap,kd)*h_LOS(:,ap2,kd)'*h_LOS(:,ap2,k) );  % diff pilots case
                            end
                        end
                    end
                end
                temp = interference_temp1 + interference_temp2;
                INTERFERENCE_UAV_GUE(k,kd) = abs(temp);
                interference_sum = interference_sum + abs(temp);   % total interference
            end
        end
        for kd = 1+K_mmW:K
            interference_temp1 =0;
            interference_temp2 =0;
            for ap =1:L
                if round( PHI(:,kd)'* PHI(:,k)) ==1   % same pilots j in P_k \k
                    interference_temp1 = interference_temp1 +  eta(ap,kd)*( h_LOS(:,ap,k)'*h_LOS(:,ap,kd)*h_LOS(:,ap,kd)'*h_LOS(:,ap,k)...
                        + eta_p*h_LOS(:,ap,k)'*h_LOS(:,ap,kd)*trace(psi(:,:,ap,k)*R(:,:,ap,kd)*R(:,:,ap,k)) +eta_p*h_LOS(:,ap,k)'*R(:,:,ap,kd)*psi(:,:,ap,k)*R(:,:,ap,kd)*h_LOS(:,ap,k)...
                        + eta_p*trace(psi(:,:,ap,k)*R(:,:,ap,k)*h_LOS(:,ap,kd)*h_LOS(:,ap,kd)'*R(:,:,ap,k)) +eta_p*trace( psi(:,:,ap,k)*R(:,:,ap,k)*R(:,:,ap,kd))*h_LOS(:,ap,kd)'*h_LOS(:,ap,k)...
                        + (eta_p^2)*trace( psi(:,:,ap,k)*R(:,:,ap,k)*R(:,:,ap,kd))^2 + (eta_p^2)*trace( psi(:,:,ap,k)*R(:,:,ap,k)*R(:,:,ap,kd)*psi(:,:,ap,k)*R(:,:,ap,kd)*R(:,:,ap,k))...
                        + trace( (R(:,:,ap,k)-eta_p*R(:,:,ap,k)*psi(:,:,ap,k)*R(:,:,ap,k))*(h_LOS(:,ap,kd)* h_LOS(:,ap,kd)' + eta_p*R(:,:,ap,kd)*psi(:,:,ap,k)*R(:,:,ap,kd)) ) ) ; % p + q terms est, error
                    
                else
                    interference_temp1 = interference_temp1 + eta(ap,kd)*trace((h_LOS(:,ap,k)*h_LOS(:,ap,k)'+R(:,:,ap,k))*(h_LOS(:,ap,kd)*h_LOS(:,ap,kd)' + eta_p*R(:,:,ap,kd)*psi(:,:,ap,kd)*R(:,:,ap,kd)));
                end
                %interference_temp1 = interference_temp1 +  eta(ap,kd)*( (norm(h_LOS(:,ap,k))^2 *norm(h_LOS(:,ap,kd))^2) + (norm(h_LOS(:,ap,k))^2 * trace(R(:,:,ap,kd)))...
                %   +(norm(h_LOS(:,ap,kd))^2 * trace(R(:,:,ap,k))) + trace(R(:,:,ap,k))*trace(R(:,:,ap,kd))  );
                for ap2 = 1:L
                    if ap2 ~= ap
                        if round( PHI(:,kd)'* PHI(:,k)) ==1   % same pilots j in P_k \k
                            interference_temp2 = interference_temp2 + sqrt(eta(ap,kd)*eta(ap2,kd))* (trace(h_LOS(:,ap,kd)*h_LOS(:,ap,k)' + eta_p*R(:,:,ap,kd)*psi(:,:,ap,kd)*R(:,:,ap,k))...
                                *trace( h_LOS(:,ap2,k)*h_LOS(:,ap2,kd)' + eta_p*R(:,:,ap2,k)*psi(:,:,ap2,kd)*R(:,:,ap2,kd)) );
                        else
                            interference_temp2 = interference_temp2 + sqrt(eta(ap,kd)*eta(ap2,kd))* ( h_LOS(:,ap,k)' *h_LOS(:,ap,kd)*h_LOS(:,ap2,kd)'*h_LOS(:,ap2,k) );  % diff pilots case
                        end
                    end
                end
            end
            temp = interference_temp1 + interference_temp2;
            INTERFERENCE_UAV_GUE(k,kd) = abs(temp);
            interference_sum = interference_sum + abs(temp);   % total interference
        end
    else 
        for kd = 1:K_mmW
            interference_temp1 =0;
            interference_temp2 =0;
            for ap =1:L
                if round( PHI(:,kd)'* PHI(:,k)) ==1   % same pilots j in P_k \k
                    interference_temp1 = interference_temp1 +  eta(ap,kd)*( h_LOS(:,ap,k)'*h_LOS(:,ap,kd)*h_LOS(:,ap,kd)'*h_LOS(:,ap,k)...
                        + eta_p*h_LOS(:,ap,k)'*h_LOS(:,ap,kd)*trace(psi(:,:,ap,k)*R(:,:,ap,kd)*R(:,:,ap,k)) +eta_p*h_LOS(:,ap,k)'*R(:,:,ap,kd)*psi(:,:,ap,k)*R(:,:,ap,kd)*h_LOS(:,ap,k)...
                        + eta_p*trace(psi(:,:,ap,k)*R(:,:,ap,k)*h_LOS(:,ap,kd)*h_LOS(:,ap,kd)'*R(:,:,ap,k)) +eta_p*trace( psi(:,:,ap,k)*R(:,:,ap,k)*R(:,:,ap,kd))*h_LOS(:,ap,kd)'*h_LOS(:,ap,k)...
                        + (eta_p^2)*trace( psi(:,:,ap,k)*R(:,:,ap,k)*R(:,:,ap,kd))^2 + (eta_p^2)*trace( psi(:,:,ap,k)*R(:,:,ap,k)*R(:,:,ap,kd)*psi(:,:,ap,k)*R(:,:,ap,kd)*R(:,:,ap,k))...
                        + trace( (R(:,:,ap,k)-eta_p*R(:,:,ap,k)*psi(:,:,ap,k)*R(:,:,ap,k))*(h_LOS(:,ap,kd)* h_LOS(:,ap,kd)' + eta_p*R(:,:,ap,kd)*psi(:,:,ap,k)*R(:,:,ap,kd)) ) ) ; % p + q terms est, error
                    
                else
                    interference_temp1 = interference_temp1 + plos(kd)*eta(ap,kd)*trace((h_LOS(:,ap,k)*h_LOS(:,ap,k)'+R(:,:,ap,k))*(h_LOS(:,ap,kd)*h_LOS(:,ap,kd)' + eta_p*R(:,:,ap,kd)*psi(:,:,ap,kd)*R(:,:,ap,kd)));
                end
                %interference_temp1 = interference_temp1 +  eta(ap,kd)*( (norm(h_LOS(:,ap,k))^2 *norm(h_LOS(:,ap,kd))^2) + (norm(h_LOS(:,ap,k))^2 * trace(R(:,:,ap,kd)))...
                %   +(norm(h_LOS(:,ap,kd))^2 * trace(R(:,:,ap,k))) + trace(R(:,:,ap,k))*trace(R(:,:,ap,kd))  );
                for ap2 = 1:L
                    if ap2 ~= ap
                        if round( PHI(:,kd)'* PHI(:,k)) ==1   % same pilots j in P_k \k
                            interference_temp2 = interference_temp2 + sqrt(eta(ap,kd)*eta(ap2,kd))* (trace(h_LOS(:,ap,kd)*h_LOS(:,ap,k)' + eta_p*R(:,:,ap,kd)*psi(:,:,ap,kd)*R(:,:,ap,k))...
                                *trace( h_LOS(:,ap2,k)*h_LOS(:,ap2,kd)' + eta_p*R(:,:,ap2,k)*psi(:,:,ap2,kd)*R(:,:,ap2,kd)) );
                        else
                            interference_temp2 = interference_temp2 + plos(kd)*sqrt(eta(ap,kd)*eta(ap2,kd))* ( h_LOS(:,ap,k)' *h_LOS(:,ap,kd)*h_LOS(:,ap2,kd)'*h_LOS(:,ap2,k) );  % diff pilots case
                        end
                    end
                end
            end
            temp = interference_temp1 + interference_temp2;
            INTERFERENCE_UAV_GUE(k,kd) = abs(temp);
            interference_sum = interference_sum + abs(temp);   % total interference
        end
        for kd = 1+K_mmW:K
            if (kd~=k)
                interference_temp1 =0;
                interference_temp2 =0;
                for ap =1:L
                    if round( PHI(:,kd)'* PHI(:,k)) ==1   % same pilots j in P_k \k
                        interference_temp1 = interference_temp1 +  eta(ap,kd)*( h_LOS(:,ap,k)'*h_LOS(:,ap,kd)*h_LOS(:,ap,kd)'*h_LOS(:,ap,k)...
                            + eta_p*h_LOS(:,ap,k)'*h_LOS(:,ap,kd)*trace(psi(:,:,ap,k)*R(:,:,ap,kd)*R(:,:,ap,k)) +eta_p*h_LOS(:,ap,k)'*R(:,:,ap,kd)*psi(:,:,ap,k)*R(:,:,ap,kd)*h_LOS(:,ap,k)...
                            + eta_p*trace(psi(:,:,ap,k)*R(:,:,ap,k)*h_LOS(:,ap,kd)*h_LOS(:,ap,kd)'*R(:,:,ap,k)) +eta_p*trace( psi(:,:,ap,k)*R(:,:,ap,k)*R(:,:,ap,kd))*h_LOS(:,ap,kd)'*h_LOS(:,ap,k)...
                            + (eta_p^2)*trace( psi(:,:,ap,k)*R(:,:,ap,k)*R(:,:,ap,kd))^2 + (eta_p^2)*trace( psi(:,:,ap,k)*R(:,:,ap,k)*R(:,:,ap,kd)*psi(:,:,ap,k)*R(:,:,ap,kd)*R(:,:,ap,k))...
                            + trace( (R(:,:,ap,k)-eta_p*R(:,:,ap,k)*psi(:,:,ap,k)*R(:,:,ap,k))*(h_LOS(:,ap,kd)* h_LOS(:,ap,kd)' + eta_p*R(:,:,ap,kd)*psi(:,:,ap,k)*R(:,:,ap,kd)) ) ) ; % p + q terms est, error
                        
                    else
                        interference_temp1 = interference_temp1 + eta(ap,kd)*trace((h_LOS(:,ap,k)*h_LOS(:,ap,k)'+R(:,:,ap,k))*(h_LOS(:,ap,kd)*h_LOS(:,ap,kd)' + eta_p*R(:,:,ap,kd)*psi(:,:,ap,kd)*R(:,:,ap,kd)));
                    end
                    %interference_temp1 = interference_temp1 +  eta(ap,kd)*( (norm(h_LOS(:,ap,k))^2 *norm(h_LOS(:,ap,kd))^2) + (norm(h_LOS(:,ap,k))^2 * trace(R(:,:,ap,kd)))...
                    %   +(norm(h_LOS(:,ap,kd))^2 * trace(R(:,:,ap,k))) + trace(R(:,:,ap,k))*trace(R(:,:,ap,kd))  );
                    for ap2 = 1:L
                        if ap2 ~= ap
                            if round( PHI(:,kd)'* PHI(:,k)) ==1   % same pilots j in P_k \k
                                interference_temp2 = interference_temp2 + sqrt(eta(ap,kd)*eta(ap2,kd))* (trace(h_LOS(:,ap,kd)*h_LOS(:,ap,k)' + eta_p*R(:,:,ap,kd)*psi(:,:,ap,kd)*R(:,:,ap,k))...
                                    *trace( h_LOS(:,ap2,k)*h_LOS(:,ap2,kd)' + eta_p*R(:,:,ap2,k)*psi(:,:,ap2,kd)*R(:,:,ap2,kd)) );
                            else
                                interference_temp2 = interference_temp2 + sqrt(eta(ap,kd)*eta(ap2,kd))* ( h_LOS(:,ap,k)' *h_LOS(:,ap,kd)*h_LOS(:,ap2,kd)'*h_LOS(:,ap2,k) );  % diff pilots case
                            end
                        end
                    end
                end
                temp = interference_temp1 + interference_temp2;
                INTERFERENCE_UAV_GUE(k,kd) = abs(temp);
                interference_sum = interference_sum + abs(temp);   % total interference
            end
        end
    end
    beamforming_uncer_cvx = 0;  beamforming_uncer_cvx2=0;
    for ap = 1:L
        beamforming_uncer_cvx = beamforming_uncer_cvx + eta(ap,k)*(trace(GAMMA_NLOS(:,:,ap,k)^2) ...
            + 2*h_LOS(:,ap,k)'*GAMMA_NLOS(:,:,ap,k)*h_LOS(:,ap,k) + trace(C_ERR(:,:,ap,k)*gamma_MAT(:,:,ap,k)) );
        beamforming_uncer_cvx2 = beamforming_uncer_cvx2 + eta(ap,k)*( abs(trace(GAMMA_NLOS(:,:,ap,k)))^2 + trace(GAMMA_NLOS(:,:,ap,k)^2) + abs(h_LOS(:,ap,k)'*h_LOS(:,ap,k))^2 ...
            + 2*h_LOS(:,ap,k)'*GAMMA_NLOS(:,:,ap,k)*h_LOS(:,ap,k) + 2*real(trace(GAMMA_NLOS(:,:,ap,k))*h_LOS(:,ap,k)'*h_LOS(:,ap,k)) + trace(C_ERR(:,:,ap,k)*gamma_MAT(:,:,ap,k)) - gamma(ap,k)^2 );

    end
    
    Interfsum(k) = interference_sum ;
    beamforming(k) = abs(beamforming_uncer_cvx);
    
    snr_numerator1(k) = snr_numerator;
    snr_denominator(k) = abs(interference_sum +  HI_term_AP_tx1(k) + HI_term_UE_rx1(k) + beamforming_uncer_cvx  +1);
    SE_k7(k) = log2(1+ snr_numerator/snr_denominator(k));
end
end