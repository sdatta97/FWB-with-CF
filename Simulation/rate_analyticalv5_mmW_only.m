% function rate_dl = rate_analytical(params, plos2, plos)
% function rate_dl = rate_analyticalv4(params, plos2, plos, R_GUE,h_LOS_GUE, PLOS_GUE)
function rate_dl = rate_analyticalv5_mmW_only(params, ue_idx)
 %%  define
N = params.num_antennas_per_gNB;  % antennas per AP
N_mmW = params.num_antennas_per_gNB_mmW;  % antennas per AP
N_UE = params.num_antennas_per_UE_mmW;
L = params.numGNB_sub6;
K = params.numUE + params.numUE_sub6;  % --Ground UEs
K_mmW = params.numUE;
snr_db = params.snr_db;
snr_db_mmw = params.snr_db_mmw;
LOOP = length(params.snr_db);
asd_length = length(params.ASD_VALUE);
hi_length = length(params.Kt_Kr_vsUE);
ASD_VALUE = params.ASD_VALUE;
ASD_CORR = params.ASD_CORR;
Kt_Kr_vsUE = params.Kt_Kr_vsUE;
K_Factor = params.K_Factor;
RAYLEIGH=params.RAYLEIGH;   %1= rayleigh, % 0=rician
Perf_CSI = params.Perf_CSI;
cov_area = params.cov_area;
%%
TAU_P_K_by_two = params.TAU_P_K_by_two;  
CH_estimation = params.CH_estimation;  
%%
LB = params.LB;  %Lower bound
UB = params.UB;  %Upper bound
no_of_rea = params.no_of_rea;     % no.of channel realizations
%%
pilot_pow = params.pilot_pow; 
noiseFigure = params.noiseFigure;
sigma_sf = params.sigma_sf;
Band = params.Band; %Communication bandwidth
tau_c = params.tau_c;      % coherence block length
NMSE = zeros(LOOP,1);
SE_UB_each = zeros(K,LOOP,asd_length,hi_length);  %UPPER-bound
SE_LB_each = zeros(K,LOOP,asd_length,hi_length);  %lower-bound
SE_LB = zeros(LOOP,asd_length,hi_length);  %lower-bound
SE_monte_impCSI = zeros(LOOP,asd_length,hi_length);  %upper bound

SNR_NUM_LB7 = zeros(K,LOOP,asd_length,hi_length);  %lower-bou
SNR_DEN_LB7= zeros(K,LOOP,asd_length,hi_length);
HI_UE_rx7= zeros(K,LOOP,asd_length,hi_length);
HI_AP_tx7= zeros(K,LOOP,asd_length,hi_length);
BU7= zeros(K,LOOP,asd_length,hi_length);
INTERFERENCE_UAV_GUE_EACH7= zeros(K,K,LOOP,asd_length,hi_length);     
% [channelGain_GUE,R_GUE,h_LOS_GUE,K_Rician,PLOS_GUE] = channel_cellfree_GUE3(K,L,N,ASD_VALUE,ASD_CORR,RAYLEIGH,0,K_Factor,cov_area,Band);
channelGain_GUE = params.BETA;
K_Rician = params.ricianFactor;
R10_mmW = params.R_GUE_mmW;
h_LOSall_mmW = params.h_LOS_GUE_mmW;
R10 = params.R_GUE;
h_LOSall = params.h_LOS_GUE;
%% PLOS cal
%         % GUE
%         CAL_PLOS_GUE = (PLOS_GUE==1);
%         PLOS_GUE_ones = sum(CAL_PLOS_GUE(:))/(L*K)*100;
%
%         PLOS_GUE_g0p9 = (PLOS_GUE>0.9);
%         PLOS_GUE_g0p9 = sum(PLOS_GUE_g0p9(:))/(L*K)*100;
%
%         PLOS_GUE_g0p8 = (PLOS_GUE>0.8);
%         PLOS_GUE_g0p8 = sum(PLOS_GUE_g0p8(:))/(L*K)*100;
%
%         % GUE rician modelling
%         CAL_Krici_GUE = (K_Rician>5); %plos=0.833
%         Krici_GUE_g_5 = sum(CAL_Krici_GUE(:))/(L*K)*100;

%%
for iter = 1:LOOP    
    for iHI = 1:hi_length          
        Kt_Kr_val = Kt_Kr_vsUE (iHI);  % K_t,K-r
        k_t2 = Kt_Kr_val;       %AP tr. H.impairment factor
        k_r2 = Kt_Kr_val;       %AP rx. H.impairment factor
        k_t2_UE =  Kt_Kr_val;    %UE tr HW
        k_r2_UE =  Kt_Kr_val;    %UE rx HW
        
        for iASD = 1: asd_length
            %%
            R_mmW= R10_mmW(:,:,:,:,iASD);
            disp(['ASD:'   num2str(iASD)]);
            h_LOS_mmW = h_LOSall_mmW;
            R= R10(:,:,:,:,iASD);
            disp(['ASD:'   num2str(iASD)]);
            h_LOS = h_LOSall(:,:,:);            % orthogonal pilot seq
            if TAU_P_K_by_two == 1
                tau_p = K/2;  %no. of pilots per coher block
            else
                tau_p = K;
            end
            eta_p = tau_p*pilot_pow;
            tau_factor = 1-(tau_p/tau_c);
            
            PHI = zeros(tau_p,K);
            PHI1   =  orth(rand(tau_p));
            %% GUE-GUE pilot contamination
            if TAU_P_K_by_two == 1
                phi_indexUG = repmat(randperm(tau_p),1,2); % take reaming half pilots for GUEs
            else
                phi_indexUG = randperm(tau_p);
            end
            for k=1:K
                PHI(:,k)=PHI1(:,phi_indexUG(k));   %PHI -- tau_p X (K-columns)               % orthogonal pilot seq
            end
            %%
            % CHANNEL GENERATION, ESTIMATION
            [h_mmW,h_hat_HI_mmW,psi_HI_mmW,h_sub6,h_hat_HI_sub6,psi_HI_sub6]= function_channel_Generation_HI(N,N_mmW,L,K,K_mmW,R_mmW,h_LOS_mmW,R,h_LOS,PHI,tau_p,pilot_pow,k_r2,k_t2_UE,no_of_rea);%(N,N_mmW,L,K,K_mmW,R,h_LOS,PHI,tau_p,pilot_pow,k_r2,k_t2_UE,no_of_rea);
            %% EST_CHANNEL GAIN
            gamma_mmW = zeros(L,K_mmW);
            GAMMA_NLOS_mmW = zeros(N_mmW,N_mmW,L,K_mmW);
            gamma_MAT_mmW = zeros(N_mmW,N_mmW,L,K_mmW);
            beta_actual_MAT_mmW = zeros(N_mmW,N_mmW,L,K_mmW);
            beta_actual_mmW  = zeros(L,K_mmW);
            C_ERR_mmW = zeros(N_mmW,N_mmW,L,K_mmW);
            h_NORMsq_mmW = zeros(N_mmW,L,K_mmW);
            gamma_sub6 = zeros(L,K-K_mmW);
            GAMMA_NLOS_sub6 = zeros(N,N,L,K-K_mmW);
            gamma_MAT_sub6 = zeros(N,N,L,K-K_mmW);
            beta_actual_MAT_sub6 = zeros(N,N,L,K-K_mmW);
            beta_actual_sub6 = zeros(L,K-K_mmW);
            C_ERR_sub6 = zeros(N,N,L,K-K_mmW);
            h_NORMsq_sub6 = zeros(N,L,K-K_mmW);
            if CH_estimation == 1
                for ap=1:L
                    for ue=1:K_mmW
                        GAMMA_NLOS_mmW(:,:,ap,ue) = eta_p*R_mmW(:,:,ap,ue)*psi_HI_mmW(:,:,ap,ue)*R_mmW(:,:,ap,ue);
                        gamma_MAT_mmW(:,:,ap,ue) = h_LOS_mmW(:,ap,ue)*h_LOS_mmW(:,ap,ue)' + eta_p*R_mmW(:,:,ap,ue)*psi_HI_mmW(:,:,ap,ue)*R_mmW(:,:,ap,ue);
                        gamma_mmW(ap,ue) = abs(trace(gamma_MAT_mmW(:,:,ap,ue)));
                        beta_actual_MAT_mmW(:,:,ap,ue) = h_LOS_mmW(:,ap,ue)*h_LOS_mmW(:,ap,ue)' + R_mmW(:,:,ap,ue);
                        beta_actual_mmW(ap,ue) = abs(trace(beta_actual_MAT_mmW(:,:,ap,ue)));
                        C_ERR_mmW(:,:,ap,ue) = beta_actual_MAT_mmW(:,:,ap,ue)-gamma_MAT_mmW(:,:,ap,ue);
                        for n2=1:N_mmW
                            h_NORMsq_mmW(n2,ap,ue) = norm(h_LOS_mmW(n2,ap,ue))^2;
                        end
                    end
                    for ue=1:K-K_mmW
                        GAMMA_NLOS_sub6(:,:,ap,ue) = eta_p*R(:,:,ap,ue)*psi_HI_sub6(:,:,ap,ue)*R(:,:,ap,ue);
                        gamma_MAT_sub6(:,:,ap,ue) = h_LOS(:,ap,ue)*h_LOS(:,ap,ue)' + eta_p*R(:,:,ap,ue)*psi_HI_sub6(:,:,ap,ue)*R(:,:,ap,ue);
                        gamma_sub6(ap,ue) = abs(trace(gamma_MAT_sub6(:,:,ap,ue)));
                        beta_actual_MAT_sub6(:,:,ap,ue) = h_LOS(:,ap,ue)*h_LOS(:,ap,ue)' + R(:,:,ap,ue);
                        beta_actual_sub6(ap,ue) = abs(trace(beta_actual_MAT_sub6(:,:,ap,ue)));
                        C_ERR_sub6(:,:,ap,ue) = beta_actual_MAT_sub6(:,:,ap,ue)-gamma_MAT_sub6(:,:,ap,ue);
                        for n2=1:N
                            h_NORMsq_sub6(n2,ap,ue) = norm(h_LOS(n2,ap,ue))^2;
                        end
                    end
                end
            else
                % NO CHANNEL -ESTIMATOIN CASE
                psi_HI_mmW = zeros(N_mmW,N_mmW,L,K_mmW);
                R_mmW = zeros(N_mmW,N_mmW,L,K_mmW);
                psi_HI_sub6 = zeros(N,N,L,K-K_mmW);
                R = zeros(N,N,L,K-K_mmW);
                for ap=1:L
                    for ue=1:K_mmW
                        GAMMA_NLOS_mmW(:,:,ap,ue) = eta_p*R(:,:,ap,ue)*psi_HI_mmW(:,:,ap,ue)*R(:,:,ap,ue);
                        gamma_MAT_mmW(:,:,ap,ue) = h_LOS_mmW(:,ap,ue)*h_LOS_mmW(:,ap,ue)' + eta_p*R_mmW(:,:,ap,ue)*psi_HI_mmW(:,:,ap,ue)*R_mmW(:,:,ap,ue);
                        gamma_mmW(ap,ue) = abs(trace(gamma_MAT_mmW(:,:,ap,ue)));
                        
                        beta_actual_MAT_mmW(:,:,ap,ue) = h_LOS_mmW(:,ap,ue)*h_LOS_mmW(:,ap,ue)' + R_mmW(:,:,ap,ue);  % true channel-- both LoS+NLoS
                        beta_actual_mmW(ap,ue) = abs(trace(beta_actual_MAT_mmW(:,:,ap,ue)));
                        C_ERR_mmW(:,:,ap,ue) = R_mmW(:,:,ap,ue); % beta_actual_MAT(:,:,ap,ue)-gamma_MAT(:,:,ap,ue);
                        for n2=1:N_mmW
                            h_NORMsq_mmW(n2,ap,ue) = norm(h_LOS_mmW(n2,ap,ue))^2;
                        end
                    end
                    for ue=1:K-K_mmW
                        GAMMA_NLOS_sub6(:,:,ap,ue) = eta_p*R(:,:,ap,ue)*psi_HI_sub6(:,:,ap,ue)*R(:,:,ap,ue);
                        gamma_MAT_sub6(:,:,ap,ue) = h_LOS(:,ap,ue)*h_LOS(:,ap,ue)' + eta_p*R(:,:,ap,ue)*psi_HI_sub6(:,:,ap,ue)*R(:,:,ap,ue);
                        gamma_sub6(ap,ue) = abs(trace(gamma_MAT_sub6(:,:,ap,ue)));
                        
                        beta_actual_MAT_sub6(:,:,ap,ue) = h_LOS(:,ap,ue)*h_LOS(:,ap,ue)' + R(:,:,ap,ue);  % true channel-- both LoS+NLoS
                        beta_actual_sub6(ap,ue) = abs(trace(beta_actual_MAT_sub6(:,:,ap,ue)));
                        C_ERR_sub6(:,:,ap,ue) = R(:,:,ap,ue); % beta_actual_MAT(:,:,ap,ue)-gamma_MAT(:,:,ap,ue);
                        for n2=1:N
                            h_NORMsq_sub6(n2,ap,ue) = norm(h_LOS(n2,ap,ue))^2;
                        end
                    end
                end
            end
            %%
            if Perf_CSI == 1
                GAMMA_NLOS_mmW = R_mmW;
                gamma_MAT_mmW = beta_actual_MAT_mmW;
                gamma_mmW = beta_actual_mmW;
                h_hat_HI_mmW=h_mmW;
                GAMMA_NLOS_sub6 = R;
                gamma_MAT_sub6 = beta_actual_MAT_sub6;
                gamma_sub6 = beta_actual_sub6;
                h_hat_HI_sub6=h_sub6;
            end
            
            %%
            snr_linear = db2pow(snr_db(iter));
            disp([' SNR: '   num2str(snr_db(iter))]);
            snr = snr_linear;
            snr_linear_mmw = db2pow(snr_db_mmw(iter));
            disp([' SNR mmW: '   num2str(snr_db_mmw(iter))]);
            snr_mmw = snr_linear_mmw;            % power allocation
            eta = zeros(L,K);
            %% allocate equal power to all
            for ap = 1:L
                summ = sum(gamma_mmW(ap,:));
                for k=1:K_mmW
                    if gamma_mmW(ap,k)== 0
                        eta(ap,k) =  0;
                    else
                        eta(ap,k) =  snr_mmw/summ;
                        %                             eta(ap,k) =  snr; %snr/K;
                    end
                end
                summ = sum(gamma_sub6(ap,:));
                for k=1:K-K_mmW
                    if gamma_sub6(ap,k)== 0
                        eta(ap,k) =  0;
                    else
                        eta(ap,k) =  snr/summ;
                        %                             eta(ap,k) =  snr; %snr/K;
                    end
                end
            end
            %% monte--carlo
            % perfect-CSI monte-carlo (Upper bound)
            % Imp-CSI monte-carlo (Upper bound)
            % [SE_UB, SE_num_UB, SE_den_UB, HI_UE_rx_UB, HI_AP_tr_UB ] = function_monte_carlo(L,K,K_mmW, N,eta,h,h_hat_HI,k_t2,k_r2_UE,no_of_rea,plos);
            [SE_UB, SE_num_UB, SE_den_UB, HI_UE_rx_UB, HI_AP_tr_UB ] = function_monte_carlo(L,K,K_mmW,N,N_mmW,eta,h_mmW,h_hat_HI_mmW,h_sub6,h_hat_HI_sub6,k_t2,k_r2_UE,no_of_rea,plos);
            SE_monte_impCSI(iter,iASD,iHI) = tau_factor*sum(SE_UB);
            SE_UB_each(1:K,iter,iASD,iHI) = tau_factor*SE_UB;
            % ImpCSI -- LB, closed form iCSI
            % [SE_LB_ALL, SNR_NUM_LB7(1:K,iter,iASD,iHI), SNR_DEN_LB7(1:K,iter,iASD,iHI), HI_UE_rx7(1:K,iter,iASD,iHI), HI_AP_tx7(1:K,iter,iASD,iHI), BU7(1:K,iter,iASD,iHI),INTERFERENCE_UAV_GUE_EACH7(1:K,1:K,iter,iASD,iHI)] = function_LB_impCSI(K_mmW,K,L,N,eta,h_LOS,R,psi_HI,eta_p,PHI,k_t2,k_r2_UE,gamma, gamma_MAT, beta_actual, beta_actual_MAT, C_ERR, GAMMA_NLOS, plos, plos2);            
            % SE_LB_each(1:K,iter,iASD,iHI) = tau_factor*SE_LB_ALL;
            % SE_LB(iter,iASD,iHI) = tau_factor*sum(SE_LB_ALL);                                       
        end  
    end
end
SNR_NUM_LB9 = mean(SNR_NUM_LB7,5);     NUMM = squeeze(SNR_NUM_LB9);
SNR_DEN_LB9 = mean(SNR_DEN_LB7,5);     denn = squeeze(SNR_DEN_LB9);
HI_UE_rx9 = mean(HI_UE_rx7,5);         HI_UEE = squeeze(HI_UE_rx9); 
HI_AP_tx9 = mean(HI_AP_tx7,5);        HI_APP = squeeze(HI_AP_tx9);
BU9  = mean(BU7,5);                     BUU = squeeze(BU9);
%% GUE
NUMM_GUE_avg9 = NUMM(1:K,:);        NUMM_GUE_avg = mean(NUMM_GUE_avg9,1);
DENN_GUE_avg9 = denn(1:K,:);        DENN_GUE_avg = mean(DENN_GUE_avg9,1);
HI_UEE_GUE_avg9 = HI_UEE(1:K,:);    HI_UEE_GUE_avg = mean(HI_UEE_GUE_avg9,1);
HI_APP_GUE_avg9 = HI_APP(1:K,:);    HI_APP_GUE_avg = mean(HI_APP_GUE_avg9,1);
BU_GUE_avg9 = BUU(1:K,:);           BU_GUE_avg = mean(BU_GUE_avg9,1);

INTERFERENCE_UAV_GUE_EACH9 = mean(INTERFERENCE_UAV_GUE_EACH7,6);
INTERFERENCE_GUE_EACH9 = INTERFERENCE_UAV_GUE_EACH9(1:K,:,:,:,:);
INTERFERENCE_GUE_from_GUEs9 = squeeze(INTERFERENCE_GUE_EACH9(:,1:K,:,:,:));
%%  IMPORTANT -- for each iteration-- HW
INTERFERENCE_GUE_from_GUEs99 = sum(INTERFERENCE_GUE_from_GUEs9,2);
INTERFERENCE_GUE_from_GUEs = mean(squeeze(INTERFERENCE_GUE_from_GUEs99),1);
Total_interference_GUEs_add = INTERFERENCE_GUE_from_GUEs;  %just interference-- NO-HW-NO BU
%% 
SUM_SE_setup = squeeze(SE_LB);
SE_LB_avg = mean(SE_LB,4);
sum_SE_LB=squeeze(SE_LB_avg); %sum over K, % imp CSI LB

SE_UB_avg = mean(SE_monte_impCSI,4);
sum_SE_UB=squeeze(SE_UB_avg);
% rate_dl = Band*sum_SE_LB/K;
% rate_dl = Band*mean(SE_LB_each,2);
rate_dl = Band*mean(SE_UB_each,2);
end