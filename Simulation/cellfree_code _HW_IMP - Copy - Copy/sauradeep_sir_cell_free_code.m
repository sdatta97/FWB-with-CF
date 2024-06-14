%% Copyright @2018 Ekant and Dheeraj
% UPLINK NMSE FOR CELL FREE WITH RF IMPAIRMENTS AND DYNAMIC ADC/DAC IMPAIRMENTS
clc;
%clear all;
%close all;

%% CONSTANTS

UE_SETUP      = 1;     % No of random user location setups
CHAN_ITER     = 1;     % No of monte-carlo realizations for each user setup.
B             = 20e6;   % Bandwidth
T_C           = 200;    % Coherence interval


% Correlation Paramenters
ASDdeg    = 10;
ZETA      = 0.0;

% Select the value of large scale fading coefficients
Sigma_option = 1; % if, Beta_option = 0 then, beta is all ones vector, if Beta_option = 1 then beta is generated from practical system values

% For perfect/imperfect CSI
pft_csi=0; % pft_csi=1 for perfect CSI case and 0 for imperfect CSI case

% Select the correlation matrix model
CORR_option = 1; % 0 implies uncorrelated; 1 implies Local scattering model;

% Select whether all users are Rician or Rayleigh
Rician  = 1;

%Select uplink NMSE or downlink NMSE
UPLINK_NMSE   =0;
DOWNLINK_NMSE =0;
RATE          =1;
GEE_EPA       =0;
GEE_OPA       =0;

PARTIAL_RATE  =1;
TOTAL_RATE    =0;

UP_RATE       =0;
DOWN_RATE_WCH =1;
DOWN_RATE_CH  =0;
UP_EE         =0;
DOWN_EE_WCH   =0;
DOWN_EE_CH    =0;

% Select what you want to plot i.e., NMSE Vs Pilot Power or Vs Pilot length
PILOT_POWER  =0;
PILOT_LENGTH =0;
AP_ANT       =0;
TX_SNR       =1;
M_AP         =0;
K_USER       =0;


% Select FD or HD
FD = 0; % Select FD=1 if Full duplex and FD=0 for Half Duplex

% Pre-log-factor
if FD==1
    pre_log_fact = 1;
elseif FD==0
    pre_log_fact = 0.5;
end

% Loop_interference suppression factor
Gamma_LI = 10^(-20/10);

%% Variances
noiseFigure      = 9; %Noise figure at the AP (in dB)
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%% Hardware Impairment parameters
ALPHA_DU    = 1-0.0;    %User DAC Impairments
ALPHA_AU    = 1-0.0;    %User ADC Impairments
ALPHA_AAP   = 1-0.03534;    %AP ADC Impairment
ALPHA_DAP   = 1-0.03534;    %AP DAC Impairment
KAPPA_TU    = 0.1; % Transmit residual HW Impairments at user
KAPPA_RU    = 0.1; % Receive residual HW Impairments at user
KAPPA_RAP   = 0.1; % Receive residual HW Impairments at AP
KAPPA_TAP   = 0.1; % Transmit residual HW Impairments at AP
SIGMA_NAP   = 1;  % Noise Variance at the APs
SIGMA_NU    = 1;  %Noise variance at the users
SIGMA_NC    = 1;  %Noise variance at the CPU
%% Simulation parameters initialization
if UPLINK_NMSE==1
    M_PLOT        = 64;  M=M_PLOT;   % Number of APs
    N_PLOT        = 10;  N=N_PLOT;   % Number of antennas at AP
    K_PLOT        = 5;   K=K_PLOT;   % Number of user
    if (PILOT_POWER == 1)
        P_P_UP_dB_PLOT  =[-20:10:60];         % Uplink Pilot power in dB
        TAU_UP_PLOT   =K;
        loop          =length(P_P_UP_dB_PLOT);
    elseif (PILOT_LENGTH == 1)
        P_P_UP_dB_PLOT  = 10;
        TAU_UP_PLOT     = [1:1:10].*K;        % Uplink Pilot length
        loop            = length(TAU_UP_PLOT);
    end
elseif DOWNLINK_NMSE==1
    M_PLOT          = 64;  M=M_PLOT;   % Number of APs
    N_PLOT          = 10;  N=N_PLOT;   % Number of antennas at AP
    K_PLOT          = 5;   K=K_PLOT;   % Number of user
    P_P_UP_dB_PLOT  = 10;
    TAU_UP_PLOT     = K;
    if (PILOT_POWER == 1)
        P_P_DOWN_dB_PLOT  = [-40:20:40];         % Downlink Pilot power in dB
        TAU_DOWN_PLOT     = K;
        loop              = length(P_P_DOWN_dB_PLOT);
    elseif (PILOT_LENGTH == 1)
        P_P_DOWN_dB_PLOT   = 10;
        TAU_DOWN_PLOT      = [1:1:10].*K;        % Downlink Pilot length
        loop               = length(TAU_DOWN_PLOT);
    end
elseif RATE==1||GEE_EPA==1||GEE_OPA==1
    P_P_UP_dB_PLOT   = 40;
    P_P_DOWN_dB_PLOT = 40;
    if AP_ANT==1
        N_PLOT          = [10:20:90]; %No. of AP antennas
        M_PLOT          = 64;
        K_PLOT          = 5;
        P_DATA_dB_PLOT  = 10;
        loop            = length(N_PLOT);
    elseif M_AP==1
        N_PLOT         = 4;
        M_PLOT         = [16,32,64,128]; %No. of APs
        K_PLOT         = 5;
        P_DATA_dB_PLOT = 10;
        loop           = length(M_PLOT);
    elseif K_USER==1
        N_PLOT         = 4;
        M_PLOT         = 64;
        K_PLOT         = [5:5:25]; %No. of Users
        P_DATA_dB_PLOT = 10;
        loop           = length(K_PLOT);
    elseif TX_SNR==1
        N_PLOT         = 4;
        M_PLOT         = 64;
        K_PLOT         = 5;
        P_DATA_dB_PLOT = [-60:10:20]; %Transmit power
        loop           = length(P_DATA_dB_PLOT);
    elseif PILOT_POWER==1
        N_PLOT           = 4;
        M_PLOT           = 64;
        K_PLOT           = 5;
        P_DATA_dB_PLOT   = 20; 
        P_P_UP_dB_PLOT   = [-40:10:20]; %Pilot power
        P_P_DOWN_dB_PLOT = P_P_UP_dB_PLOT;
        loop             = length(P_P_UP_dB_PLOT);        
    end
TAU_UP_PLOT    = round(max(K_PLOT))-1;
TAU_DOWN_PLOT  = TAU_UP_PLOT-1;
end

%%
%initialization of matrices to store the NMSE values for different user setups
Nmax=max(N_PLOT);
Mmax=max(M_PLOT);
Kmax=max(K_PLOT);

UP_NMSE_TH_SETUP     = zeros(UE_SETUP, loop);   % Theoritical
UP_NMSE_MC_SETUP     = zeros(UE_SETUP, loop);   % Monte-Carlo
DOWN_NMSE_TH_SETUP   = zeros(UE_SETUP, loop);   % Theoritical
DOWN_NMSE_MC_SETUP   = zeros(UE_SETUP, loop);   % Monte-Carlo
SUM_SE_TH_SETUP      = zeros(UE_SETUP, loop);   % Theoritical
SUM_SE_MC_SETUP      = zeros(UE_SETUP, loop);   % Monte-Carlo
GEE_TH_SETUP         = zeros(UE_SETUP, loop);   % Theoritical
GEE_OPT_TH_SETUP    = zeros(UE_SETUP, loop);   % Monte-Carlo


R_G_SETUP       = zeros(Nmax,Nmax,Kmax,Mmax,UE_SETUP);   %Stores channel covariance matrices for max (N,K,M) and all user setups
Gmean_SETUP     = zeros(Nmax,Kmax,Mmax,UE_SETUP);        %Stores channel means for max (N,K,M) and all user setups
R_H_SETUP       = zeros(Nmax,Nmax,Kmax,Mmax,UE_SETUP);   %Stores channel covariance matrices for max (N,K,M) and all user setups
Hmean_SETUP     = zeros(Nmax,Kmax,Mmax,UE_SETUP);        %Stores channel means for max (N,K,M) and all user setups
R_LR_SETUP      = zeros(Nmax,Nmax,Mmax,UE_SETUP);        %Stores Rx inter-AP covariance matrices
R_LT_SETUP      = zeros(Nmax,Nmax,Mmax,UE_SETUP);        %Stores Tx inter-AP covariance matrices
SIGMA_LAP_SETUP = zeros(Mmax,Mmax,UE_SETUP);             %Stores inter-AP channel gains
SIGMA_LU_SETUP  = zeros(Kmax,Kmax,UE_SETUP);             %Stors inter-UE channel gains
for ue_iters = 1:UE_SETUP
    disp(['UE_SETUP = ' num2str(ue_iters)]);
    
    % Initializating vectors to store values per each setup
    UP_NMSE_MC   = zeros(1,loop);
    UP_NMSE_TH   = zeros(1,loop);
    DOWN_NMSE_MC = zeros(1,loop);
    DOWN_NMSE_TH = zeros(1,loop);
    SUM_SE_TH    = zeros(1,loop);
    SUM_SE_MC    = zeros(1,loop);
    SE_TH        = zeros(loop,Kmax);
    UP_SE_TH     = zeros(loop,Kmax);
    DOWN_SE_TH   = zeros(loop,Kmax);
    GEE_TH       = zeros(loop,1);
    GEE_OPT_TH  = zeros(loop,1);
    
    tic;
    
    
    
    % Generation of Practical channel
    [R_G_tmp, Gmean_tmp, R_H_tmp, Hmean_tmp, K_R, SIGMA_G, R_LT_tmp, R_LR_tmp, SIGMA_LAP_tmp, SIGMA_LU_tmp] = FD_channel_cellfree_rician(Mmax,Nmax,Kmax,ASDdeg,Rician,CORR_option,Sigma_option, noiseVariancedBm, FD);
    %This function returns all the channel correlation matrices (NxNxMxK) and channel means (NxMxK) (with respect to each user, each AP)
    %The FD loop interfernce parameters are: R_LT_tmp, R_LR_tmp (NxNxMxM)(Tx and Rx AP covariance matrices), SIGMA_LAP (MxM), SIGMA_LU (KxK) loop interference large scale fading coefficients
    SIGMA_H=SIGMA_G;
    % The large scale fading coefficients remain same for the plink and the downlink channels
    
    R_G_SETUP(:,:,:,:,ue_iters)   = R_G_tmp;       %Storing channel covariance matrices for each UE setup
    Gmean_SETUP(:,:,:,ue_iters)   = Gmean_tmp;     %Storing channel mean for each UE setup
    R_H_SETUP(:,:,:,:,ue_iters)   = R_H_tmp;       %Storing channel covariance matrices for each UE setup
    Hmean_SETUP(:,:,:,ue_iters)   = Hmean_tmp;     %Storing channel mean for each UE setup
    R_LR_SETUP(:,:,:,ue_iters)    = R_LT_tmp;      %Storing Rx inter-AP covariance matrices for each setup
    R_LT_SETUP(:,:,:,ue_iters)    = R_LR_tmp;      %Storing Tx inter-AP covariance matrices for each setup
    SIGMA_LAP_SETUP(:,:,ue_iters) = Gamma_LI*SIGMA_LAP_tmp*FD; %Storing inter-AP channel gains
    SIGMA_LU_SETUP(:,:,ue_iters)  = SIGMA_LU_tmp*FD;  %Storing inter-UE channel gains
    
    
    for n = 1:loop % going through all values of SNR/N/K in each setup
        %Setting the plot parameter
        if UPLINK_NMSE==1
            if (PILOT_POWER == 1)
                P_P_UP     = 10.^(P_P_UP_dB_PLOT(n)./10);  % Pilot power linear
                TAU_UP   = TAU_UP_PLOT;                % Pilot length (To estimate one of the channels)
            elseif (PILOT_LENGTH == 1)
                P_P_UP     = 10.^(P_P_UP_dB_PLOT./10);      % Pilot power linear
                TAU_UP   = TAU_UP_PLOT(n);              % Pilot length
            end
        elseif DOWNLINK_NMSE==1
            P_P_UP     = 10.^(P_P_UP_dB_PLOT./10);  % Uplink Pilot power linear
            TAU_UP   = TAU_UP_PLOT;             % Uplink Pilot length
            if (PILOT_POWER == 1)
                P_P_DOWN     = 10.^(P_P_DOWN_dB_PLOT(n)./10);  % Downlink Pilot power linear
                TAU_DOWN   = TAU_DOWN_PLOT;                % Downlink Pilot length
            elseif (PILOT_LENGTH == 1)
                P_P_DOWN     = 10.^(P_P_DOWN_dB_PLOT./10);      % Downlink Pilot power linear
                TAU_DOWN   = TAU_DOWN_PLOT(n);              % Downlink Pilot length
            end
        elseif RATE==1||GEE_EPA==1||GEE_OPA==1
            P_P_UP      = 10.^(P_P_UP_dB_PLOT./10);    % Uplink Pilot power linear
            TAU_UP      = TAU_UP_PLOT;                 % Uplink Pilot length
            P_P_DOWN   = 10.^(P_P_DOWN_dB_PLOT./10);  % Downlink Pilot power linear
            TAU_DOWN    = TAU_DOWN_PLOT;               % Downlink Pilot length
            if AP_ANT==1
                N      = N_PLOT(n);            % No. of AP antennas
                M      = M_PLOT;               % No. of APs
                K      = K_PLOT;               % No. of users
                P_DATA = 10^(P_DATA_dB_PLOT./10); % Transmit SNR
            elseif M_AP==1
                N      = N_PLOT;                % No. of AP antennas
                M      = M_PLOT(n);             % No. of APs
                K      = K_PLOT;                % No. of users
                P_DATA = 10^(P_DATA_dB_PLOT./10); % Transmit SNR
            elseif K_USER==1
                N      = N_PLOT;                % No. of AP antennas
                M      = M_PLOT;                % No. of APs
                K      = K_PLOT(n);             % No. of users
                P_DATA = 10^(P_DATA_dB_PLOT./10); % Transmit SNR
            elseif TX_SNR==1
                N      = N_PLOT;                   % No. of AP antennas
                M      = M_PLOT;                   % No. of APs
                K      = K_PLOT;                   % No. of users
                P_DATA = 10^(P_DATA_dB_PLOT(n)./10); % Transmit SNR
            elseif PILOT_POWER==1
                N      = N_PLOT;                   % No. of AP antennas
                M      = M_PLOT;                   % No. of APs
                K      = K_PLOT;                   % No. of users
                P_DATA = 10^(P_DATA_dB_PLOT./10);  % Transmit SNR
                P_P_UP     = 10.^(P_P_UP_dB_PLOT(n)./10);    % Uplink Pilot power linear
                P_P_DOWN   = 10.^(P_P_DOWN_dB_PLOT(n)./10);  % Downlink Pilot power linear
            end
        end
        
        TAU_TUP = 2*TAU_UP;               % Total pilot length for uplink channel estimation (for both G and H)
        if DOWN_RATE_WCH==1||DOWN_EE_WCH==1
            TAU_DOWN=0;
        end
        TAU_D   = T_C-(TAU_TUP+TAU_DOWN); % Data length
        
        
        %Select channel mean and channel covariance matrices
        Gmean     = Gmean_SETUP(1:N,1:K,1:M,ue_iters);
        R_G       = R_G_SETUP(1:N,1:N,1:K,1:M,ue_iters);
        Hmean     = Hmean_SETUP(1:N,1:K,1:M,ue_iters);
        R_H       = R_H_SETUP(1:N,1:N,1:K,1:M,ue_iters);
        R_LR      = R_LR_SETUP(1:N,1:N,1:M,ue_iters);
        R_LT      = R_LT_SETUP(1:N,1:N,1:M,ue_iters);
        SIGMA_LAP = SIGMA_LAP_SETUP(1:M,1:M,ue_iters);
        SIGMA_LUE = SIGMA_LU_SETUP(1:K,1:K,ue_iters);
        
        
        
        %% ADC/DAC impairments
        b_AP   = 2*ones(N,M); %ADC/DAC resolution at APs
        b_U    = 2*ones(K,1); %ADC/DAC resolution at Users
        A_A_AP = zeros(N,N,M);
        A_D_AP = zeros(N,N,M);
        for l = 1:M
            A_A_AP(:,:,l) = ALPHA_AAP*eye(N);  %ADC impairments for all the APs
            A_D_AP(:,:,l) = ALPHA_DAP*eye(N);  %DAC impairments for all the APs
        end
        A_A_U  = ALPHA_AU*ones(K,1); %ADC impairments for all the users
        A_D_U  = ALPHA_DU*ones(K,1); %DAC impairments for all the users
        
        %% Pilot power coefficient
        BETA_UP_P    = ones(K,1); %Equal uplink pilot power coefficient
        BETA_DOWN_P  = ones(K,M)./(M*K); %Equal downlink pilot power coefficient (This is same as Downlink data power coefficient)
        BETA_UP_DATA = ones(K,1)./(K); %Equal uplink data power coefficient
        
        
        %% Theoritical UPLINK/DOWNLINK NMSE Calculation
        [R_G_HAT, C_G, CHI_G, LAMBDA_G, R_G_HAT_S, R_H_HAT, C_H, CHI_H, LAMBDA_H, R_H_HAT_S] = FD_UP_CSI_statistics_cellfree( M, K, N, TAU_UP, P_P_UP, BETA_UP_P, R_G, Gmean, R_H, Hmean, A_D_U, A_A_AP, KAPPA_TU, KAPPA_RAP, SIGMA_NAP, pft_csi);
        %This function calculates the covariance matrix of estimated channel and error covariance matrices
        %This function estimates both the uplink and downlink channels (G and H) usnig the uplink pilots
        
        % Downlink power coefficients
        P_constraint=zeros(M,1);
        for ap=1:M
            for ue_1=1:K
                for ue_2=1:K
                    if rem((ue_1-ue_2),TAU_UP)==0
                        P_constraint(ap)= P_constraint(ap) + (Hmean(:,ue_1,ap)'*Hmean(:,ue_2,ap)+trace(R_H_HAT(:,:,ue_1,ap)*LAMBDA_H(:,:,ue_1,ue_2,ap)));
                    else
                        P_constraint(ap)= P_constraint(ap) + (Hmean(:,ue_1,ap)'*Hmean(:,ue_2,ap));
                    end
                end
            end
        end
        for ap=1:M
            for ue=1:K
                BETA_DOWN_P(ue,ap)=1/(P_constraint(ap));
            end
        end
        %---------------------------------------------------
        
        
      
                    [a_bar, R_a, C_ya, C_yy, R_a_HAT, C_a, C_D_AP, R_a_11_sum, R_a_12_sum, R_a_13_sum, R_a_14_sum] = FD_DOWN_CSI_statistics_cellfree(TAU_UP,K,M,N, A_D_AP,A_A_U,KAPPA_TAP,KAPPA_RU,SIGMA_NU, P_P_DOWN,BETA_DOWN_P, Hmean,R_H_HAT,R_H_HAT_S,C_H,LAMBDA_H, pft_csi);
                    %This fynction calculates the variance of the estimated channel and the error variance
                    [SINR_TH] = FD_Down_ch_Theoritical_RATE(K, P_DATA, A_A_U, A_D_U, KAPPA_RU, KAPPA_TU, SIGMA_NU, SIGMA_LUE, TAU_DOWN, a_bar, R_a_HAT, C_a, R_a_11_sum, R_a_12_sum, R_a_13_sum, R_a_14_sum, BETA_UP_DATA,FD);
                    %This function returns the Theoritical SINR at all the K users
                               
          
        
        
        if UPLINK_NMSE==1
            for k=1:K
                for m=1:M
                    UP_NMSE_TH(n) = UP_NMSE_TH(n) + trace(C_G(:,:,k,m))/(trace(Gmean(:,k,m)*Gmean(:,k,m)'+R_G(:,:,k,m))*M*K);%Evaluate Theoritical NMSE using the formula for each P_p/tau_p
                end
            end
        elseif DOWNLINK_NMSE==1
            for ue_1=1:K
                for ue_2=1:K
                    DOWN_NMSE_TH(n) = DOWN_NMSE_TH(n) + C_a(ue_1,ue_2)/(a_bar(ue_1,ue_2)*conj(a_bar(ue_1,ue_2))+R_a(ue_1,ue_2,ue_2));%Evaluate Theoritical NMSE using the formula for each P_p/tau_p
                end
            end
        elseif RATE==1
            if PARTIAL_RATE==1
                for ue=1:K
                    SE_TH(n,ue) = pre_log_fact*(TAU_D/(T_C))*log2(1+SINR_TH(ue));
                end
            elseif TOTAL_RATE==1
                for ue=1:K
                    SE_TH(n,ue) = pre_log_fact*(TAU_D/(T_C))*(log2(1+UP_SINR_TH(ue))+log2(1+DOWN_SINR_TH(ue)));
                end
            end
            SUM_SE_TH(n) = sum(SE_TH(n,:));
        end
        
        %-------------------------------------------------GEE------------------------------------------------
        
        
        
        %% Monte-carlo UPLINK NMSE Calculation
        UP_NMSE_SIM   = zeros(CHAN_ITER,1);
        DOWN_NMSE_SIM = zeros(CHAN_ITER,1);
        SINR_MC       = zeros(CHAN_ITER,K);
        SE_MC         = zeros(CHAN_ITER,K);
        for ch = 1:CHAN_ITER
            [G,H, G_LAP]    = FD_random_channel_generator(M,K,N, R_G,Gmean, R_H, Hmean,R_LR,R_LT,SIGMA_LAP); %Actual channel generated
            %This function generates original channels using the R and Gmean
            if pft_csi==0
                [G_HAT,H_HAT] = FD_UP_MMSE_channel_estimator(TAU_UP, M,K,N, P_P_UP,BETA_UP_P, G,H,Gmean,Hmean,R_G,R_H, KAPPA_TU,KAPPA_RAP,A_D_U,A_A_AP,SIGMA_NAP, CHI_G,CHI_H);
                %This function generates the uplink channel estimates of G and H
            elseif pft_csi==1
                G_HAT = G;
                H_HAT = H;
            end
            
            
            
            
               
                    [a_hat, a] = FD_DOWN_MMSE_channel_estimator(M,K,N, BETA_DOWN_P,P_P_DOWN, H,H_HAT, KAPPA_TAP,KAPPA_RU,A_D_AP,A_A_U,C_D_AP,TAU_UP,SIGMA_NU, a_bar,C_ya,C_yy);
                    %This function generates the downlink channel estimates
                    if pft_csi==1
                        a_hat=a;
                    end
                    [SINR_MC(ch,:)] = FD_Down_wch_Montecarlo_Rate(N, K, M, TAU_UP, P_DATA, BETA_DOWN_P, A_A_U, KAPPA_RU, C_D_AP, SIGMA_NU, H,H_HAT, a, SIGMA_LUE, A_D_U, KAPPA_TU,BETA_UP_DATA, FD);
                    % This function returns the monte-carlo SINR of all the users
               
                
                
            
            %Calculation of Uplink NMSE Monte-carlo
            if UPLINK_NMSE==1
                for k=1:K
                    for m=1:M
                        UP_NMSE_SIM(ch) = UP_NMSE_SIM(ch) + norm(G(:,k,m)-G_HAT(:,k,m))^2/(norm(G(:,k,m))^2*(M*K)); %Calculated NMSE over each channel iteration
                    end
                end
                %Calculation of Downlink NMSE Monte-carlo
            elseif DOWNLINK_NMSE==1
                for k_1=1:K
                    for k_2=1:K
                        DOWN_NMSE_SIM(ch) = DOWN_NMSE_SIM(ch) + abs(a_hat(k_1,k_2)-a(k_1,k_2))^2/(abs(a(k_1,k_2))^2*K); %Calculated NMSE over each channel iteration
                    end
                end
                %Calculation of montecarlo rate
            elseif RATE==1
                if PARTIAL_RATE==1
                    for k=1:K
                        SE_MC(ch,k)=pre_log_fact*(TAU_D/(T_C))*log2(1+SINR_MC(ch,k)); %Evaluate Ergodic SE for each user and each channel iteratioin
                    end
                elseif TOTAL_RATE==1
                    for k=1:K
                        SE_MC(ch,k)=pre_log_fact*(TAU_D/(T_C))*(log2(1+UP_SINR_MC(ch,k))+log2(1+DOWN_SINR_MC(ch,k))); %Evaluate Ergodic SE for each user and each channel iteratioin
                    end
                end
            end
        end
        
        % Averaging over Channel Iterations
        if UPLINK_NMSE==1
            UP_NMSE_MC(n)=sum(UP_NMSE_SIM)/CHAN_ITER;  % Evaluate Ergodic Uplink NMSE for each P_p/tau_p
        elseif DOWNLINK_NMSE==1
            DOWN_NMSE_MC(n)=sum(DOWN_NMSE_SIM)/CHAN_ITER;  % Evaluate Ergodic Downlik NMSE for each P_p/tau_p
        elseif RATE==1
            SUM_SE_MC(n)=sum(sum(SE_MC))/CHAN_ITER; %Evaluate Ergodic SE
        end
        
        disp(['n = ' num2str(n) '/' num2str(loop) ' is done.']);
        
        %% Figure for UPLINK/DOWNLINK NMSE vs Pilot power or Pilot length
        if UPLINK_NMSE==1
            figure(ue_iters)
            if (PILOT_POWER == 1)
                semilogy(P_P_UP_dB_PLOT(1:n),UP_NMSE_MC(1:n),'-hr');
                hold on
                semilogy(P_P_UP_dB_PLOT(1:n),UP_NMSE_TH(1:n),'-^b');
                hold on
                grid on
                drawnow
                legend('Monte Carlo','Theoritical')
                title(['UE SETUP =',num2str(ue_iters)])
            elseif (PILOT_LENGTH == 1)
                semilogy(TAU_UP_PLOT(1:n),real(UP_NMSE_MC(1:n)),'-hr');
                hold on
                semilogy(TAU_UP_PLOT(1:n),real(UP_NMSE_TH(1:n)),'-^b');
                hold on
                grid on
                drawnow
                legend('Monte Carlo','Theoritical')
                title(['UE SETUP =',num2str(ue_iters)])
            end
        elseif DOWNLINK_NMSE==1
            figure(ue_iters)
            if (PILOT_POWER == 1)
                semilogy(P_P_DOWN_dB_PLOT(1:n),DOWN_NMSE_MC(1:n),'-hr');
                hold on
                semilogy(P_P_DOWN_dB_PLOT(1:n),DOWN_NMSE_TH(1:n),'-^b');
                hold on
                grid on
                drawnow
                legend('Monte Carlo','Theoritical')
                title(['UE SETUP =',num2str(ue_iters)])
            elseif (PILOT_LENGTH == 1)
                semilogy(TAU_DOWN_PLOT(1:n),real(DOWN_NMSE_MC(1:n)),'-hr');
                hold on
                semilogy(TAU_DOWN_PLOT(1:n),real(DOWN_NMSE_TH(1:n)),'-^b');
                hold on
                grid on
                drawnow
                legend('Monte Carlo','Theoritical')
                title(['UE SETUP =',num2str(ue_iters)])
            end
        elseif RATE==1
            if AP_ANT==1
                plot(N_PLOT(1:n),SUM_SE_MC(1:n),'-hr');
                hold on
                plot(N_PLOT(1:n),SUM_SE_TH(1:n),'-^b');
                hold on
                drawnow
                legend('simulated','lwr bnd')
                title(['UE SETUP =',num2str(ue_iters)])
            elseif M_AP==1
                plot(M_PLOT(1:n),SUM_SE_MC(1:n),'-hr');
                hold on
                plot(M_PLOT(1:n),SUM_SE_TH(1:n),'-^b');
                hold on
                drawnow
                legend('simulated','lwr bnd')
                title(['UE SETUP =',num2str(ue_iters)])
            elseif K_USER==1
                plot(K_PLOT(1:n),SUM_SE_MC(1:n),'-hr');
                hold on
                plot(K_PLOT(1:n),SUM_SE_TH(1:n),'-^b');
                hold on
                drawnow
                legend('simulated','lwr bnd')
                title(['UE SETUP =',num2str(ue_iters)])
            elseif TX_SNR==1
                %plot(P_DATA_dB_PLOT(1:n),SUM_SE_MC(1:n),'-hr');
                hold on
                plot(P_DATA_dB_PLOT(1:n),SUM_SE_TH(1:n),'-^r');
                hold on
                drawnow
                legend('simulated','lwr bnd')
                title(['UE SETUP =',num2str(ue_iters)])
            elseif PILOT_POWER==1
                %plot(P_P_UP_dB_PLOT(1:n),SUM_SE_MC(1:n),'-hr');
                hold on
                plot(P_P_UP_dB_PLOT(1:n),SUM_SE_TH(1:n),'-^b');
                hold on
                drawnow
                legend('simulated','lwr bnd')
                title(['UE SETUP =',num2str(ue_iters)])               
            end
        elseif GEE_EPA==1
            if AP_ANT==1
                plot(N_PLOT(1:n),GEE_TH(1:n),'-^b');
                hold on
                drawnow
                legend('lwr bnd')
                title(['UE SETUP =',num2str(ue_iters)])
            elseif M_AP==1
                plot(M_PLOT(1:n),GEE_TH(1:n),'-^r');
                hold on
                drawnow
                legend('lwr bnd')
                title(['UE SETUP =',num2str(ue_iters)])
            elseif K_USER==1
                plot(K_PLOT(1:n),GEE_TH(1:n),'-^b');
                hold on
                drawnow
                legend('lwr bnd')
                title(['UE SETUP =',num2str(ue_iters)])
            elseif TX_SNR==1
                plot(P_DATA_dB_PLOT(1:n),GEE_TH(1:n),'-^b');
                hold on
                drawnow
                legend('lwr bnd')
                title(['UE SETUP =',num2str(ue_iters)])
            end
        elseif GEE_OPA==1
            if AP_ANT==1
                plot(N_PLOT(1:n),GEE_TH(1:n),'-^b');
                hold on
                plot(N_PLOT(1:n),GEE_OPT_TH(1:n),'-r');
                hold on
                drawnow
                legend('lwr bnd')
                title(['UE SETUP =',num2str(ue_iters)])
            elseif M_AP==1
                plot(M_PLOT(1:n),GEE_TH(1:n),'-b');
                hold on
                plot(M_PLOT(1:n),GEE_OPT_TH(1:n),'-r');
                hold on
                drawnow
                legend('lwr bnd')
                title(['UE SETUP =',num2str(ue_iters)])
            elseif K_USER==1
                plot(K_PLOT(1:n),GEE_TH(1:n),'-b');
                hold on
                plot(K_PLOT(1:n),GEE_OPT_TH(1:n),'-r');
                hold on
                drawnow
                legend('lwr bnd')
                title(['UE SETUP =',num2str(ue_iters)])
            elseif TX_SNR==1
                plot(P_DATA_dB_PLOT(1:n),GEE_TH(1:n),'-^b');
                hold on
                plot(P_DATA_dB_PLOT(1:n),GEE_OPT_TH(1:n),'-^r');
                hold on
                drawnow
                legend('lwr bnd')
                title(['UE SETUP =',num2str(ue_iters)])
            end
        end
    end
    %storing the values of SE obtained for each setup.
    if UPLINK_NMSE==1
        UP_NMSE_TH_SETUP(ue_iters,:)  = UP_NMSE_TH;
        UP_NMSE_MC_SETUP(ue_iters,:)  = UP_NMSE_MC;
    elseif DOWNLINK_NMSE==1
        DOWN_NMSE_TH_SETUP(ue_iters,:)  = DOWN_NMSE_TH;
        DOWN_NMSE_MC_SETUP(ue_iters,:)  = DOWN_NMSE_MC;
    elseif RATE==1
        SUM_SE_TH_SETUP(ue_iters,:)  = SUM_SE_TH;
        SUM_SE_MC_SETUP(ue_iters,:)  = SUM_SE_MC;
    elseif GEE_EPA==1
        GEE_TH_SETUP(ue_iters,:)  = GEE_TH;
    elseif GEE_OPA==1
        GEE_TH_SETUP(ue_iters,:)      = GEE_TH;
        GEE_OPT_TH_SETUP(ue_iters,:)  = GEE_OPT_TH;
    end
    toc;
end
if UE_SETUP>1
    if UPLINK_NMSE==1
        if (PILOT_POWER == 1)
            figure;
            semilogy(P_P_UP_dB_PLOT,mean(UP_NMSE_TH_SETUP),'-hr');
            hold on
            semilogy(P_P_UP_dB_PLOT,mean(UP_NMSE_MC_SETUP),'-^b');
            hold on
            drawnow
            legend('lwr bnd','simulated')
        elseif (PILOT_LENGTH == 1)
            figure;
            semilogy(TAU_UP_PLOT,mean(UP_NMSE_MC_SETUP),'-hr');
            hold on
            semilogy(TAU_UP_PLOT,mean(UP_NMSE_TH_SETUP),'-^b');
            hold on
            drawnow
            legend('simulated','lwr bnd')
        end
    elseif DOWNLINK_NMSE==1
        if (PILOT_POWER == 1)
            figure;
            semilogy(P_P_UP_dB_PLOT,mean(DOWN_NMSE_TH_SETUP),'-hr');
            hold on
            semilogy(P_P_UP_dB_PLOT,mean(DOWN_NMSE_MC_SETUP),'-^b');
            hold on
            drawnow
            legend('lwr bnd','simulated')
        elseif (PILOT_LENGTH == 1)
            figure;
            semilogy(TAU_UP_PLOT,mean(DOWN_NMSE_MC_SETUP),'-hr');
            hold on
            semilogy(TAU_UP_PLOT,mean(DOWN_NMSE_TH_SETUP),'-^b');
            hold on
            drawnow
            legend('simulated','lwr bnd')
        end
    elseif RATE ==1
        if AP_ANT==1
            plot(N_PLOT,mean(SUM_SE_MC_SETUP),'-hr');
            hold on
            plot(N_PLOT,mean(SUM_SE_TH_SETUP),'-^b');
            hold on
            drawnow
            legend('lwr bnd','simulated')
        elseif M_AP==1
            plot(M_PLOT,mean(SUM_SE_MC_SETUP),'-hr');
            hold on
            plot(M_PLOT,mean(SUM_SE_TH_SETUP),'-^b');
            hold on
            drawnow
            legend('lwr bnd','simulated')
        elseif K_USER==1
            plot(K_PLOT,mean(SUM_SE_MC_SETUP),'-hr');
            hold on
            plot(K_PLOT,mean(SUM_SE_TH_SETUP),'-^b');
            hold on
            drawnow
            legend('lwr bnd','simulated')
        elseif TX_SNR==1
            plot(P_DATA_dB_PLOT,mean(SUM_SE_MC_SETUP),'-hr');
            hold on
            %plot(P_DATA_dB_PLOT,mean(SUM_SE_TH_SETUP),'-^b');
            hold on
            drawnow
            legend('lwr bnd','simulated')
        elseif PILOT_POWER==1
            plot(P_P_UP_dB_PLOT,mean(SUM_SE_MC_SETUP),'-hr');
            hold on
            %plot(P_P_UP_dB_PLOT,mean(SUM_SE_TH_SETUP),'-^b');
            hold on
            drawnow
            legend('lwr bnd','simulated')
        end
    elseif GEE_EPA ==1
        if AP_ANT==1
            plot(N_PLOT,mean(GEE_TH_SETUP),'-^b');
            hold on
            drawnow
            legend('lwr bnd')
        elseif M_AP==1
            plot(M_PLOT,mean(GEE_TH_SETUP),'-^b');
            hold on
            drawnow
            legend('lwr bnd')
        elseif K_USER==1
            plot(K_PLOT,mean(GEE_TH_SETUP),'-^b');
            hold on
            drawnow
            legend('simulated')
        elseif TX_SNR==1
            plot(P_DATA_dB_PLOT,mean(GEE_TH_SETUP),'-^b');
            hold on
            drawnow
            legend('simulated')
        end
    elseif GEE_OPA ==1
        if AP_ANT==1
            plot(N_PLOT,mean(GEE_TH_SETUP),'-^b');
            hold on
            plot(N_PLOT,mean(GEE_OPT_TH_SETUP),'-r');
            hold on
            drawnow
            legend('lwr bnd')
        elseif M_AP==1
            plot(M_PLOT,mean(GEE_TH_SETUP),'-^b');
            hold on
            plot(M_PLOT,mean(GEE_OPT_TH_SETUP),'-r');
            hold on
            drawnow
            legend('lwr bnd')
        elseif K_USER==1
            plot(K_PLOT,mean(GEE_TH_SETUP),'-^b');
            hold on
            plot(K_PLOT,mean(GEE_OPT_TH_SETUP),'-r');
            hold on
            drawnow
            legend('simulated')
        elseif TX_SNR==1
            plot(P_DATA_dB_PLOT,mean(GEE_TH_SETUP),'-^b');
            hold on
            plot(P_DATA_dB_PLOT,mean(GEE_OPT_TH_SETUP),'-^r');
            hold on
            drawnow
            legend('simulated')
        end
    end
end
toc;