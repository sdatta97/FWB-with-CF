function [h_mmW,h_hat_HI_mmW,psi_HI_mmW]= function_channel_Generation_HI_mmW_only(N_mmW,N_UE,L,K_mmW,R_mmW,h_LOS_mmW,PHI,tau_p,p,k_r2,k_t2_UE,no_of_rea)

    %% channel generation
warning('off','all');
    %-------------------ACTUAL CHANNEL
    h_mmW = zeros(N_mmW,N_UE,no_of_rea,L,K_mmW);
    % h = zeros(N,no_of_rea,L,K);
    for ch_re = 1: no_of_rea
        for l = 1:L
            for k1 = 1:K_mmW
                h_mmW(:,:,ch_re,l,k1) = h_LOS_mmW(:,:,l,k1)+ sqrtm(R_mmW(:,:,l,k1))*sqrt(0.5)*(randn(N_mmW,N_UE)+1i*randn(N_mmW,N_UE));    %if consider spatial correlation,....  for all ch_realizations all at once
            end
        end
    end
    
    % CORRELATION of actual CHANNEL  ---- E( g_ka g_ka^H )
    CORR_actual_channel_mmW = zeros(N_mmW,N_mmW,L,K_mmW);
    for ap =1:L
        for ue = 1:K_mmW
            CORR_actual_channel_mmW(:,:,ap,ue) = h_LOS_mmW(:,:,ap,ue)*h_LOS_mmW(:,:,ap,ue)' + R_mmW(:,:,ap,ue);
        end
    end
    % COVARIENCE OF EST. CHANNEL (Y_a * PHI(:,k) )

    psi_HI_mmW = zeros(N_mmW,N_mmW,L,K_mmW);   %cov matrix of pilot signal with HARDWARE IMPARIMENTS (only RX- side distortion)
    C_HI_mmW = zeros(N_mmW,N_mmW,L,K_mmW);     % COV OF CHANNEL ERROR with H_impari
    for ap = 1:L
        for ue = 1:K_mmW
            psi_ka_mmW = zeros(N_mmW,N_mmW);
            RX_Distortion_cov_mmW = zeros(N_mmW,N_mmW);  % H.Impairment
            user_Distortion_cov_mmW = zeros(N_mmW,N_mmW); 
            for ue_d = 1:K_mmW
                if round(PHI(:,ue)'*PHI(:,ue_d)) ==1
                    % psi_ka = psi_ka +  (tau_p*p)*R(:,:,ap,ue_d);  %SAME AS SDEEP
                    psi_ka_mmW = psi_ka_mmW +  (tau_p*p)*R_mmW(:,:,ap,ue_d);  %SAME AS SDEEP
                end
                RX_Distortion_cov_mmW = RX_Distortion_cov_mmW + k_r2*(p+k_t2_UE)*diag( diag(CORR_actual_channel_mmW(:,:,ap,ue_d))); 
                user_Distortion_cov_mmW = user_Distortion_cov_mmW + k_t2_UE*(p)*(h_LOS_mmW(:,:,ap,ue_d)*h_LOS_mmW(:,:,ap,ue_d)' + R_mmW(:,:,ap,ue_d));
            end
            psi_ka_mmW = psi_ka_mmW + eye(N_mmW);  % tau_p*eye(N); --- is NOISE %SAME AS SDEEP
            psi_HI_mmW(:,:,ap,ue) =  inv(psi_ka_mmW + RX_Distortion_cov_mmW + user_Distortion_cov_mmW);  %SAME AS SDEEP 
            C_HI_mmW(:,:,ap,ue) = R_mmW(:,:,ap,ue)- (tau_p*p)*R_mmW(:,:,ap,ue)*psi_HI_mmW(:,:,ap,ue)*R_mmW(:,:,ap,ue);
        end
    end
    
    %% ESTIMATED CHANNEL  //   % h (N,no_of_rea,AP,UE)
    
    y_p_mmW = zeros(N_mmW,N_UE,no_of_rea,L,K_mmW);
    y_p_mean_mmW = zeros(N_mmW,N_UE,no_of_rea,L,K_mmW);
    %----------
    y_p_HI_mmW = zeros(N_mmW,N_UE,no_of_rea,L,K_mmW);
    h_hat_HI_mmW = zeros(N_mmW,N_UE,no_of_rea,L,K_mmW);
    HI_UE_tr_mult_mmW = zeros(N_mmW,N_UE,K_mmW,K_mmW);
    %--------
    
    for chreali = 1:no_of_rea
        for ap = 1:L
            %--------------------with impairment
            if k_r2 ~=0
                COV_rx_distortion_mmW = zeros(N_mmW);
                for ue1 =1:K_mmW
                    COV_rx_distortion_mmW = COV_rx_distortion_mmW + k_r2*(tau_p*p + k_t2_UE)* diag( diag(h_mmW(:,:,chreali,ap,ue1)*h_mmW(:,:,chreali,ap,ue1)' ));  % to cal Upper Bound
                end
                RX_distortion_mmW =  sqrt(0.5)*sqrtm(COV_rx_distortion_mmW)*( randn(N_mmW,1) +1i*randn(N_mmW,1));   %N x Tau_p
            else
                RX_distortion_mmW = zeros(N_mmW,1);
            end
            %------------------
            for k = 1:K_mmW
                for i =1:K_mmW
                    y_p_mmW(:,:,chreali,ap,k) = y_p_mmW(:,:,chreali,ap,k) + tau_p*sqrt(p)*h_mmW(:,:,chreali,ap,i)*(PHI(:,i)'*PHI(:,k))  ;            % pilot signal with PHI(:,k)
                    y_p_mean_mmW(:,:,chreali,ap,k) = y_p_mean_mmW(:,:,chreali,ap,k) + tau_p*sqrt(p)*h_LOS_mmW(:,:,ap,i)*(PHI(:,i)'*PHI(:,k));
                    HI_UE_tr_mult_mmW(:,:,k,i) = h_mmW(:,:,chreali,ap,i)*sqrt(k_t2_UE*tau_p*p)*(randn(1,1)+1j*randn(1,1))*sqrt(0.5);    %                  
                end
                Noise_phi_mmW = sqrtm(tau_p*eye(N_mmW))*sqrt(0.5)*(randn(N_mmW,N_UE)+1j*randn(N_mmW,N_UE));
                y_p_HI_mmW(:,:,chreali,ap,k)= y_p_mmW(:,:,chreali,ap,k)+ RX_distortion_mmW+ sum(HI_UE_tr_mult_mmW(:,:,k,:),4) + Noise_phi_mmW;   %HImparim
                h_hat_HI_mmW(:,:,chreali,ap,k) = h_LOS_mmW(:,:,ap,k)+ sqrt(p)*R_mmW(:,:,ap,k)*psi_HI_mmW(:,:,ap,k)*( y_p_HI_mmW(:,:,chreali,ap,k)-y_p_mean_mmW(:,:,chreali,ap,k));
            end
        end
    end
end
     
 