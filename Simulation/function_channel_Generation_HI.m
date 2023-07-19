function [h_mmW,h_hat_HI_mmW,psi_HI_mmW,h_sub6,h_hat_HI_sub6,psi_HI_sub6]= function_channel_Generation_HI(N,N_mmW,L,K,K_mmW,R_mmW,h_LOS_mmW,R,h_LOS,PHI,tau_p,p,k_r2,k_t2_UE,no_of_rea)

    %% channel generation
warning('off','all');
    %-------------------ACTUAL CHANNEL
    h_mmW = zeros(N_mmW,no_of_rea,L,K);
    % h = zeros(N,no_of_rea,L,K);
    h_sub6 = zeros(N,no_of_rea,L,K);
    for ch_re = 1: no_of_rea
        for l = 1:L
            for k1 = 1:K
                h_mmW(:,ch_re,l,k1) = h_LOS_mmW(:,l,k1)+ sqrtm(R_mmW(:,:,l,k1))*sqrt(0.5)*(randn(N_mmW,1)+1i*randn(N_mmW,1));    %if consider spatial correlation,....  for all ch_realizations all at once
                h_sub6(:,ch_re,l,k1) = h_LOS(:,l,k1)+ sqrtm(R(:,:,l,k1))*sqrt(0.5)*(randn(N,1)+1i*randn(N,1));    %if consider spatial correlation,....  for all ch_realizations all at once
            end
        end
    end
    
    % CORRELATION of actual CHANNEL  ---- E( g_ka g_ka^H )
    CORR_actual_channel_mmW = zeros(N_mmW,N_mmW,L,K);
    CORR_actual_channel_sub6 = zeros(N,N,L,K);
    for ap =1:L
        for ue =1:K
            CORR_actual_channel_mmW(:,:,ap,ue) = h_LOS_mmW(:,ap,ue)*h_LOS_mmW(:,ap,ue)' + R_mmW(:,:,ap,ue);
            CORR_actual_channel_sub6(:,:,ap,ue) = h_LOS(:,ap,ue)*h_LOS(:,ap,ue)' + R(:,:,ap,ue);
        end
    end
    % COVARIENCE OF EST. CHANNEL (Y_a * PHI(:,k) )

    psi_HI_mmW = zeros(N_mmW,N_mmW,L,K);   %cov matrix of pilot signal with HARDWARE IMPARIMENTS (only RX- side distortion)
    C_HI_mmW = zeros(N_mmW,N_mmW,L,K);     % COV OF CHANNEL ERROR with H_impari
    psi_HI_sub6 = zeros(N,N,L,K);   %cov matrix of pilot signal with HARDWARE IMPARIMENTS (only RX- side distortion)
    C_HI_sub6 = zeros(N,N,L,K);     % COV OF CHANNEL ERROR with H_impari   
    for ap = 1:L
        for ue = 1:K
            psi_ka_mmW = zeros(N_mmW,N_mmW);
            RX_Distortion_cov_mmW = zeros(N_mmW,N_mmW);  % H.Impairment
            user_Distortion_cov_mmW = zeros(N_mmW,N_mmW); 
            psi_ka_sub6 = zeros(N,N);
            RX_Distortion_cov_sub6 = zeros(N,N);  % H.Impairment
            user_Distortion_cov_sub6 = zeros(N,N);             
            for ue_d = 1:K
                if round(PHI(:,ue)'*PHI(:,ue_d)) ==1
                    % psi_ka = psi_ka +  (tau_p*p)*R(:,:,ap,ue_d);  %SAME AS SDEEP
                    psi_ka_mmW = psi_ka_mmW +  (tau_p*p)*R_mmW(:,:,ap,ue_d);  %SAME AS SDEEP
                    psi_ka_sub6 = psi_ka_sub6 +  (tau_p*p)*R(:,:,ap,ue_d);  %SAME AS SDEEP
                end
                RX_Distortion_cov_mmW = RX_Distortion_cov_mmW + k_r2*(p+k_t2_UE)*diag( diag(CORR_actual_channel_mmW(:,:,ap,ue_d))); 
                user_Distortion_cov_mmW = user_Distortion_cov_mmW + k_t2_UE*(p)*(h_LOS_mmW(:,ap,ue_d)*h_LOS_mmW(:,ap,ue_d)' + R_mmW(:,:,ap,ue_d));
                RX_Distortion_cov_sub6 = RX_Distortion_cov_sub6 + k_r2*(p+k_t2_UE)*diag( diag(CORR_actual_channel_sub6(:,:,ap,ue_d))); 
                user_Distortion_cov_sub6 = user_Distortion_cov_sub6 + k_t2_UE*(p)*(h_LOS(:,ap,ue_d)*h_LOS(:,ap,ue_d)' + R(:,:,ap,ue_d));
            end
            psi_ka_mmW = psi_ka_mmW + eye(N_mmW);  % tau_p*eye(N); --- is NOISE %SAME AS SDEEP
            psi_ka_sub6 = psi_ka_sub6 + eye(N);  % tau_p*eye(N); --- is NOISE %SAME AS SDEEP
            psi_HI_mmW(:,:,ap,ue) =  inv(psi_ka_mmW + RX_Distortion_cov_mmW + user_Distortion_cov_mmW);  %SAME AS SDEEP 
            C_HI_mmW(:,:,ap,ue) = R_mmW(:,:,ap,ue)- (tau_p*p)*R_mmW(:,:,ap,ue)*psi_HI_mmW(:,:,ap,ue)*R_mmW(:,:,ap,ue);
            psi_HI_sub6(:,:,ap,ue) =  inv(psi_ka_sub6 + RX_Distortion_cov_sub6 + user_Distortion_cov_sub6);  %SAME AS SDEEP 
            C_HI_sub6(:,:,ap,ue) = R(:,:,ap,ue)- (tau_p*p)*R(:,:,ap,ue)*psi_HI_sub6(:,:,ap,ue)*R(:,:,ap,ue);
         end
    end
    
    %% ESTIMATED CHANNEL  //   % h (N,no_of_rea,AP,UE)
    
    y_p_mmW = zeros(N_mmW,no_of_rea,L,K);
    y_p_mean_mmW = zeros(N_mmW,no_of_rea,L,K);
    %----------
    y_p_HI_mmW = zeros(N_mmW,no_of_rea,L,K);
    h_hat_HI_mmW = zeros(N_mmW,no_of_rea,L,K);
    HI_UE_tr_mult_mmW = zeros(N_mmW,K,K);
    y_p_sub6 = zeros(N,no_of_rea,L,K);
    y_p_mean_sub6 = zeros(N,no_of_rea,L,K);
    %----------
    y_p_HI_sub6 = zeros(N,no_of_rea,L,K);
    h_hat_HI_sub6 = zeros(N,no_of_rea,L,K);
    HI_UE_tr_mult_sub6 = zeros(N,K,K);
    %--------
    
    for chreali = 1:no_of_rea
        for ap = 1:L
            %--------------------with impairment
            if k_r2 ~=0
                COV_rx_distortion_mmW = zeros(N_mmW);
                COV_rx_distortion_sub6 = zeros(N);
                for ue1 =1:K
                    COV_rx_distortion_mmW = COV_rx_distortion_mmW + k_r2*(tau_p*p + k_t2_UE)* diag( diag(h_mmW(:,chreali,ap,ue1)*h_mmW(:,chreali,ap,ue1)' ));  % to cal Upper Bound
                    COV_rx_distortion_sub6 = COV_rx_distortion_sub6 + k_r2*(tau_p*p + k_t2_UE)* diag( diag( h_sub6(:,chreali,ap,ue1)*h_sub6(:,chreali,ap,ue1)' ));  % to cal Upper Bound
                end
                RX_distortion_mmW =  sqrt(0.5)*sqrtm(COV_rx_distortion_mmW)*( randn(N_mmW,1) +1i*randn(N_mmW,1));   %N x Tau_p
                RX_distortion_sub6 =  sqrt(0.5)*sqrtm(COV_rx_distortion_sub6)*( randn(N,1) +1i*randn(N,1));   %N x Tau_p
            else
                RX_distortion_mmW = zeros(N_mmW,1);
                RX_distortion_sub6 = zeros(N,1);
            end
            %------------------
            for k = 1:K
                for i =1:K
                    y_p_mmW(:,chreali,ap,k) = y_p_mmW(:,chreali,ap,k) + tau_p*sqrt(p)*h_mmW(:,chreali,ap,i)*(PHI(:,i)'*PHI(:,k))  ;            % pilot signal with PHI(:,k)
                    y_p_mean_mmW(:,chreali,ap,k) = y_p_mean_mmW(:,chreali,ap,k) + tau_p*sqrt(p)*h_LOS_mmW(:,ap,i)*(PHI(:,i)'*PHI(:,k));
                    HI_UE_tr_mult_mmW(:,k,i) = h_mmW(:,chreali,ap,i)*sqrt(k_t2_UE*tau_p*p)*(randn(1,1)+1j*randn(1,1))*sqrt(0.5);    %                  
                    y_p_sub6(:,chreali,ap,k) = y_p_sub6(:,chreali,ap,k) + tau_p*sqrt(p)*h_sub6(:,chreali,ap,i)*(PHI(:,i)'*PHI(:,k))  ;            % pilot signal with PHI(:,k)
                    y_p_mean_sub6(:,chreali,ap,k) = y_p_mean_sub6(:,chreali,ap,k) + tau_p*sqrt(p)*h_LOS(:,ap,i)*(PHI(:,i)'*PHI(:,k));
                    HI_UE_tr_mult_sub6(:,k,i) = h_sub6(:,chreali,ap,i)*sqrt(k_t2_UE*tau_p*p)*(randn(1,1)+1j*randn(1,1))*sqrt(0.5);    %
                end
                Noise_phi_mmW = sqrtm(tau_p*eye(N_mmW))*sqrt(0.5)*(randn(N_mmW,1)+1j*randn(N_mmW,1));
                Noise_phi_sub6 = sqrtm(tau_p*eye(N))*sqrt(0.5)*(randn(N,1)+1j*randn(N,1));
                y_p_HI_mmW(:,chreali,ap,k)= y_p_mmW(:,chreali,ap,k)+ RX_distortion_mmW+ sum(HI_UE_tr_mult_mmW(:,k,:),3) + Noise_phi_mmW;   %HImparim
                h_hat_HI_mmW(:,chreali,ap,k) = h_LOS_mmW(:,ap,k)+ sqrt(p)*R_mmW(:,:,ap,k)*psi_HI_mmW(:,:,ap,k)*( y_p_HI_mmW(:,chreali,ap,k)-y_p_mean_mmW(:,chreali,ap,k));
                y_p_HI_sub6(:,chreali,ap,k)= y_p_sub6(:,chreali,ap,k)+ RX_distortion_sub6+ sum(HI_UE_tr_mult_sub6(:,k,:),3) + Noise_phi_sub6;   %HImparim
                h_hat_HI_sub6(:,chreali,ap,k) = h_LOS(:,ap,k)+ sqrt(p)*R(:,:,ap,k)*psi_HI_sub6(:,:,ap,k)*( y_p_HI_sub6(:,chreali,ap,k)-y_p_mean_sub6(:,chreali,ap,k));
            end
        end
    end
     
 