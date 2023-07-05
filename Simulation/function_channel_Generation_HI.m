function [h,h_hat_HI,psi_HI]= function_channel_Generation_HI(N,L,K,R,h_LOS,PHI,tau_p,p,k_r2,k_t2_UE,no_of_rea)

    %% channel generation
warning('off','all');
    %-------------------ACTUAL CHANNEL
    h = zeros(N,no_of_rea,L,K);
    for ch_re = 1: no_of_rea
        for l = 1:L
            for k1 = 1:K
                h(:,ch_re,l,k1) = h_LOS(:,l,k1)+ sqrtm(R(:,:,l,k1))*sqrt(0.5)*(randn(N,1)+1i*randn(N,1));    %if consider spatial correlation,....  for all ch_realizations all at once
            end
        end
    end
    
    % CORRELATION of actual CHANNEL  ---- E( g_ka g_ka^H )
    CORR_actual_channel = zeros(N,N,L,K);
    for ap =1:L
        for ue =1:K
            CORR_actual_channel(:,:,ap,ue) = h_LOS(:,ap,ue)*h_LOS(:,ap,ue)' + R(:,:,ap,ue);
        end
    end
    
    % COVARIENCE OF EST. CHANNEL (Y_a * PHI(:,k) )

    psi_HI = zeros(N,N,L,K);   %cov matrix of pilot signal with HARDWARE IMPARIMENTS (only RX- side distortion)
    C_HI = zeros(N,N,L,K);     % COV OF CHANNEL ERROR with H_impari
    
    for ap = 1:L
        for ue = 1:K
            psi_ka = zeros(N,N);
            RX_Distortion_cov = zeros(N,N);  % H.Impairment
            user_Distortion_cov = zeros(N,N); 
            
            for ue_d = 1:K
                if round(PHI(:,ue)'*PHI(:,ue_d)) ==1
                    psi_ka = psi_ka +  (tau_p*p)*R(:,:,ap,ue_d);  %SAME AS SDEEP
                end
                RX_Distortion_cov = RX_Distortion_cov + k_r2*(p+k_t2_UE)*diag( diag(CORR_actual_channel(:,:,ap,ue_d))); 
                  user_Distortion_cov = user_Distortion_cov + k_t2_UE*(p)*(h_LOS(:,ap,ue_d)*h_LOS(:,ap,ue_d)' + R(:,:,ap,ue_d));
%                 user_Distortion_cov = user_Distortion_cov + k_t2_UE*(tau_p*p)* CORR_actual_channel(:,:,ap,ue_d);
            end
            
            psi_ka = psi_ka + eye(N);  % tau_p*eye(N); --- is NOISE %SAME AS SDEEP
%             psi_ka = psi_ka + eye(N);  %eye(N); --- is NOISE
            
            psi_HI(:,:,ap,ue) =  inv(psi_ka + RX_Distortion_cov + user_Distortion_cov);  %SAME AS SDEEP 
%             psi_HI(:,:,ap,ue) = pinv(psi_ka + RX_Distortion_cov + user_Distortion_cov);
           
            C_HI(:,:,ap,ue) = R(:,:,ap,ue)- (tau_p*p)*R(:,:,ap,ue)*psi_HI(:,:,ap,ue)*R(:,:,ap,ue);
        end
    end
    
    %% ESTIMATED CHANNEL  //   % h (N,no_of_rea,AP,UE)
    
    y_p = zeros(N,no_of_rea,L,K);
    y_p_mean = zeros(N,no_of_rea,L,K);
    %----------
    y_p_HI = zeros(N,no_of_rea,L,K);
    h_hat_HI = zeros(N,no_of_rea,L,K);
    HI_UE_tr_mult = zeros(N,K,K);

    %--------
    
    for chreali = 1:no_of_rea
        for ap = 1:L
            %--------------------with impairment
            if k_r2 ~=0
                COV_rx_distortion = zeros(N);
                for ue1 =1:K
                    COV_rx_distortion = COV_rx_distortion + k_r2*(tau_p*p + k_t2_UE)* diag( diag( h(:,chreali,ap,ue1)*h(:,chreali,ap,ue1)' ));  % to cal Upper Bound
                end
                RX_distortion =  sqrt(0.5)*sqrtm(COV_rx_distortion)*( randn(N,1) +1i*randn(N,1));   %N x Tau_p
            else
                RX_distortion = zeros(N,1);
            end
            %------------------
            for k = 1:K
                Noise_phi = sqrtm(tau_p*eye(N))*sqrt(0.5)*(randn(N,1)+1j*randn(N,1));
                for i =1:K
                    y_p(:,chreali,ap,k) = y_p(:,chreali,ap,k) + tau_p*sqrt(p)*h(:,chreali,ap,i)*(PHI(:,i)'*PHI(:,k))  ;            % pilot signal with PHI(:,k)
                    y_p_mean(:,chreali,ap,k) = y_p_mean(:,chreali,ap,k) + tau_p*sqrt(p)*h_LOS(:,ap,i)*(PHI(:,i)'*PHI(:,k));
                    %y_p_mean(:,chreali,ap,k) = y_p_mean(:,chreali,ap,k) + sqrt(eta_p)*h_LOS(:,ap,i)*(PHI(:,i)'*PHI(:,k));
                    % HI_UE_tr = HI_UE_tr + sqrt(k_t2_UE*p)*h(:,chreali,ap,i)*randn(1,tau_p);
                    HI_UE_tr_mult(:,k,i) = h(:,chreali,ap,i)*sqrt(k_t2_UE*tau_p*p)*(randn(1,1)+1j*randn(1,1))*sqrt(0.5);    %
                end
                
                
%                 HI_UE_tr_mult(:,:,ap) = sum(h(:,chreali,ap,:),4)*sqrt(k_t2_UE*p)*randn(1,tau_p);    %wrong  -- randn(1,tau_p) === randn(1,tau_p) +1j*randn(1,tau_p)  
                % K-times adding, wrong, because already included K-users
                
                y_p_HI(:,chreali,ap,k)= y_p(:,chreali,ap,k)+ RX_distortion+ sum(HI_UE_tr_mult(:,k,:),3) + Noise_phi;   %HImparim

                
                h_hat_HI(:,chreali,ap,k) = h_LOS(:,ap,k)+ sqrt(p)*R(:,:,ap,k)*psi_HI(:,:,ap,k)*( y_p_HI(:,chreali,ap,k)-y_p_mean(:,chreali,ap,k));
%                 h_hat_HI(:,chreali,ap,k) = h_LOS(:,ap,k)+ sqrt(tau_p*p)*R(:,:,ap,k)*psi_HI(:,:,ap,k)*( y_p_HI(:,chreali,ap,k)-y_p_mean(:,chreali,ap,k));
            end
        end
    end
     
 