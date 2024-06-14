%% Copyright @2018 Ekant and Dheeraj
% UPLINK NMSE FOR CELL FREE WITH RF IMPAIRMENTS AND DYNAMIC ADC/DAC IMPAIRMENTS
clc;
clear all;
%close all;

%% CONSTANTS

UE_SETUP      = 1;     % No of random user location setups
CHAN_ITER     =150;% 100;     % No of monte-carlo realizations for each user setup.
B             = 20e6;   % Bandwidth
T_C           = 200;    % Coherence interval
pol           = 2;  % 1 for single-polarization , 2 for dual polarization.
alpha         = 0.3;
mmse          =1;
mrc           =0;
kappa         =0.0;
DA            =0;
alpha_ADC     =1;  % 0.6366 0.8825 0.96546 0.990503 0.997501
decod_err_prob=10^(-5);
tau_b         =T_C*10^(-4);
beta=0;%qfuncinv(decod_err_prob)*log2(exp(1))*sqrt(2)/sqrt(tau_b*B);
cvx          =  2;
% Correlation Paramenters
ASDdeg    = 10;%10;
ZETA      = 0.0;

% Select the value of large scale fading coefficients
Sigma_option = 0; % if, Beta_option = 1 then, beta is all ones vector, if Beta_option = 0 then beta is generated from practical system values

% For perfect/imperfect CSI
pft_csi=0; % pft_csi=1 for perfect CSI case and 0 for imperfect CSI case

% Select the correlation matrix model
CORR_option = 1; % 0 implies uncorrelated; 1 implies Local scattering model;

% Select whether all users are Rician or Rayleigh
Rician  = 0;

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
PILOT_POWER  =1;
PILOT_LENGTH =0;
AP_ANT       =0;
TX_SNR       =0;
M_AP         =0;
K_USER       =0;
alpha1       =0;
kappa_value  =0;
alpha_ADC_value  =0;
ASD_value    =0;

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
    N_PLOT        = 4;  N=N_PLOT;   % Number of antennas at AP
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
        N_PLOT          = [2:2:8]*1; %No. of AP antennas
        M_PLOT          = 128;%64;
        K_PLOT          = 5;
        P_DATA_dB_PLOT  = 0;
        loop            = length(N_PLOT);
    elseif M_AP==1
        N_PLOT         = 2*[8,4,2,1];%8;%1*[64,32,16,8,4,2,1];
        M_PLOT         =[16,32,64,128];%[4,8,16,32,64,128,256,512];%[16,32,64];%[16,32,64,128]; %No. of APs [1,4,8,16];%
        K_PLOT         = 5;
        P_DATA_dB_PLOT = 10;
        loop           = length(M_PLOT);
%         loop=1;
    elseif K_USER==1
        N_PLOT         = 2;
        M_PLOT         = 16;
        K_PLOT         = [5:10:35]; %No. of Users [1,2,5,10];%
        P_DATA_dB_PLOT = -10;
        loop           = length(K_PLOT);
    elseif TX_SNR==1
        N_PLOT         =4;%4;
        M_PLOT         =4;
        K_PLOT         =5;
        P_DATA_dB_PLOT =-10;%[-30:10:30];%[0:10:60]; %Transmit power
        loop           = length(P_DATA_dB_PLOT);
    elseif alpha1==1
         N_PLOT         =4 ;%4;
        M_PLOT         =16;
        K_PLOT         =5;
        P_DATA_dB_PLOT =0;
        alpha_PLOT =[0:0.1:1];
        loop           = length(alpha_PLOT);
    elseif kappa_value==1
        N_PLOT         =4 ;%4;
        M_PLOT         =64;
        K_PLOT         =5;
        P_DATA_dB_PLOT =10;
        kappa_PLOT =[0:0.1:0.5];
        alpha_ADC_PLOT =[1,0.997501,0.990503,0.96546,0.8825,0.6366];
        loop           = length(kappa_PLOT);
    elseif alpha_ADC_value==1
         N_PLOT         =4 ;%4;
        M_PLOT         =32;
        K_PLOT         =5;
        P_DATA_dB_PLOT =0;
        alpha_ADC_PLOT =[1,0.997501,0.990503,0.96546,0.8825,0.6366];
        loop           = length(alpha_ADC_PLOT);
    elseif ASD_value==1
         N_PLOT         =4 ;%4;
        M_PLOT         =16;
        K_PLOT         =5;
        P_DATA_dB_PLOT =20;
        ASD_PLOT =[10:10:50];
        loop           = length(ASD_PLOT);    
    elseif PILOT_POWER==1
        N_PLOT           = 2;
        M_PLOT           = 32;
        K_PLOT           = 5;
        P_DATA_dB_PLOT   = 0; 
        P_P_UP_dB_PLOT   = [-40:10:20]; %Pilot power
        P_P_DOWN_dB_PLOT = P_P_UP_dB_PLOT;
        loop             = length(P_P_UP_dB_PLOT);        
    end
TAU_UP_PLOT    = round(max(K_PLOT))-1;
TAU_DOWN_PLOT  = TAU_UP_PLOT-1;
end
% for n1=1:loop1
%  M_PLOT=M_PLOT1(n1);
%  N_PLOT=N_PLOT1(n1);
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
%     [R_G_tmp, Gmean_tmp, R_H_tmp, Hmean_tmp, K_R, SIGMA_G, R_LT_tmp, R_LR_tmp, SIGMA_LAP_tmp, SIGMA_LU_tmp] = FD_channel_cellfree_rician(Mmax,Nmax,Kmax,ASDdeg,Rician,CORR_option,Sigma_option, noiseVariancedBm, FD);
    [R_1] = channel_cellfree_GUE_Krician_non_zero_mean(max(K_PLOT),max(M_PLOT),max(N_PLOT),ASDdeg,1,1,Sigma_option,1,1);
%     channelGain_over_noise = 10^4*channelGain_over_noise;
%     R_original = [(1-alpha)*R alpha*R;  alpha*R (1-alpha)*R];
%     R_1 =cell2mat(struct2cell( load('R_2')));
%     load('R_1');
%     N_PLOT=2;
%     R_1=R_3(1:N_PLOT,1:N_PLOT,:,:);
%     R_1 = reshape(R_1,16,16,1,5);
    %This function returns all the channel correlation matrices (NxNxMxK) and channel means (NxMxK) (with respect to each user, each AP)
    %The FD loop interfernce parameters are: R_LT_tmp, R_LR_tmp (NxNxMxM)(Tx and Rx AP covariance matrices), SIGMA_LAP (MxM), SIGMA_LU (KxK) loop interference large scale fading coefficients
%     SIGMA_H=SIGMA_G;
    % The large scale fading coefficients remain same for the plink and the downlink channels
    
%     R_G_SETUP(:,:,:,:,ue_iters)   = R_G_tmp;       %Storing channel covariance matrices for each UE setup
%     Gmean_SETUP(:,:,:,ue_iters)   = Gmean_tmp;     %Storing channel mean for each UE setup
%     R_H_SETUP(:,:,:,:,ue_iters)   = R_H_tmp;       %Storing channel covariance matrices for each UE setup
%     Hmean_SETUP(:,:,:,ue_iters)   = Hmean_tmp;     %Storing channel mean for each UE setup
%     R_LR_SETUP(:,:,:,ue_iters)    = R_LT_tmp;      %Storing Rx inter-AP covariance matrices for each setup
%     R_LT_SETUP(:,:,:,ue_iters)    = R_LR_tmp;      %Storing Tx inter-AP covariance matrices for each setup
%     SIGMA_LAP_SETUP(:,:,ue_iters) = Gamma_LI*SIGMA_LAP_tmp*FD; %Storing inter-AP channel gains
%     SIGMA_LU_SETUP(:,:,ue_iters)  = SIGMA_LU_tmp*FD;  %Storing inter-UE channel gains
    
    
    for n = 1:loop % going through all values of SNR/N/K in each setup
        %Setting the plot parameter
        
        if DOWNLINK_NMSE==1
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
                N      = N_PLOT(n);                % No. of AP antennas
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
            elseif alpha1==1
                N      = N_PLOT;                   % No. of AP antennas
                M      = M_PLOT;                   % No. of APs
                K      = K_PLOT;                   % No. of users
                P_DATA = 10^(P_DATA_dB_PLOT./10); % Transmit SNR
                alpha=alpha_PLOT(n);
            elseif kappa_value==1
                N      = N_PLOT;                   % No. of AP antennas
                M      = M_PLOT;                   % No. of APs
                K      = K_PLOT;                   % No. of users
                P_DATA = 10^(P_DATA_dB_PLOT./10); % Transmit SNR
                kappa=kappa_PLOT(n);
                alpha_ADC=alpha_ADC_PLOT(n);
            elseif alpha_ADC_value==1
                N      = N_PLOT;                   % No. of AP antennas
                M      = M_PLOT;                   % No. of APs
                K      = K_PLOT;                   % No. of users
                P_DATA = 10^(P_DATA_dB_PLOT./10); % Transmit SNR
                alpha_ADC=alpha_ADC_PLOT(n);
             elseif ASD_value==1
                N      = N_PLOT;                   % No. of AP antennas
                M      = M_PLOT;                   % No. of APs
                K      = K_PLOT;                   % No. of users
                P_DATA = 10^(P_DATA_dB_PLOT./10); % Transmit SNR
                ASDdeg=ASD_PLOT(n);   
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
        
%         [R_1] = channel_cellfree_GUE_Krician_non_zero_mean(K,M,N,ASDdeg,1,1,0,1,1);
        %Select channel mean and channel covariance matrices
        Gmean     = Gmean_SETUP(1:N,1:K,1:M,ue_iters);
        R_G       = R_G_SETUP(1:N,1:N,1:K,1:M,ue_iters);
        Hmean     = Hmean_SETUP(1:N,1:K,1:M,ue_iters);
        R_H       = R_H_SETUP(1:N,1:N,1:K,1:M,ue_iters);
        R_LR      = R_LR_SETUP(1:N,1:N,1:M,ue_iters);
        R_LT      = R_LT_SETUP(1:N,1:N,1:M,ue_iters);
        SIGMA_LAP = SIGMA_LAP_SETUP(1:M,1:M,ue_iters);
        SIGMA_LUE = SIGMA_LU_SETUP(1:K,1:K,ue_iters);
%         R =reshape(R_1,16,16,1,5) ;
        
        R=R_1(1:N,1:N,1:M,1:K);
        
        
        
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
%         [R_G_HAT, C_G, CHI_G, LAMBDA_G, R_G_HAT_S, R_H_HAT, C_H, CHI_H, LAMBDA_H, R_H_HAT_S] = FD_UP_CSI_statistics_cellfree( M, K, N, TAU_UP, P_P_UP, BETA_UP_P, R_G, Gmean, R_H, Hmean, A_D_U, A_A_AP, KAPPA_TU, KAPPA_RAP, SIGMA_NAP, pft_csi);
        %This function calculates the covariance matrix of estimated channel and error covariance matrices
        %This function estimates both the uplink and downlink channels (G and H) usnig the uplink pilots
        
        % Downlink power coefficients
%         P_constraint=zeros(M,1);
%         for ap=1:M
%             for ue_1=1:K
%                 for ue_2=1:K
%                     if rem((ue_1-ue_2),TAU_UP)==0
%                         P_constraint(ap)= P_constraint(ap) + (Hmean(:,ue_1,ap)'*Hmean(:,ue_2,ap)+trace(R_H_HAT(:,:,ue_1,ap)*LAMBDA_H(:,:,ue_1,ue_2,ap)));
%                     else
%                         P_constraint(ap)= P_constraint(ap) + (Hmean(:,ue_1,ap)'*Hmean(:,ue_2,ap));
%                     end
%                 end
%             end
%         end
%         for ap=1:M
%             for ue=1:K
%                 BETA_DOWN_P(ue,ap)=1/(P_constraint(ap));
%             end
%         end
        %---------------------------------------------------
        
        
      
%                     [a_bar, R_a, C_ya, C_yy, R_a_HAT, C_a, C_D_AP, R_a_11_sum, R_a_12_sum, R_a_13_sum, R_a_14_sum] = FD_DOWN_CSI_statistics_cellfree(TAU_UP,K,M,N, A_D_AP,A_A_U,KAPPA_TAP,KAPPA_RU,SIGMA_NU, P_P_DOWN,BETA_DOWN_P, Hmean,R_H_HAT,R_H_HAT_S,C_H,LAMBDA_H, pft_csi);
%                     %This fynction calculates the variance of the estimated channel and the error variance
%                     [SINR_TH] = FD_Down_ch_Theoritical_RATE(K, P_DATA, A_A_U, A_D_U, KAPPA_RU, KAPPA_TU, SIGMA_NU, SIGMA_LUE, TAU_DOWN, a_bar, R_a_HAT, C_a, R_a_11_sum, R_a_12_sum, R_a_13_sum, R_a_14_sum, BETA_UP_DATA,FD);
%                     %This function returns the Theoritical SINR at all the K users
                               
          
        
        
        
            
%        if DOWNLINK_NMSE==1
%             for ue_1=1:K
%                 for ue_2=1:K
%                     DOWN_NMSE_TH(n) = DOWN_NMSE_TH(n) + C_a(ue_1,ue_2)/(a_bar(ue_1,ue_2)*conj(a_bar(ue_1,ue_2))+R_a(ue_1,ue_2,ue_2));%Evaluate Theoritical NMSE using the formula for each P_p/tau_p
%                 end
%             end
%         elseif RATE==1
%             
%                 for ue=1:K
%                     SE_TH(n,ue) = pre_log_fact*(TAU_D/(T_C))*log2(1+SINR_TH(ue));
%                 end
%             
%            
%             SUM_SE_TH(n) = sum(SE_TH(n,:));
%         end
        
        %-------------------------------------------------GEE------------------------------------------------
        
%         [channelGain_over_noise,R,g_LOS,K_Rician,probLOS] = channel_cellfree_GUE_Krician_non_zero_mean(UE,AP,N,ASD_VALUE,ASD_CORR,rayleigh,BETA,K_factor,cov_area);
        
        %% Monte-carlo UPLINK NMSE Calculation
        UP_NMSE_SIM   = zeros(CHAN_ITER,1);
        DOWN_NMSE_SIM = zeros(CHAN_ITER,1);
        SINR_MC       = zeros(CHAN_ITER,K);
        SE_MC         = zeros(CHAN_ITER,K);
        G_00           = zeros(N,M,K);
        G_01           = zeros(N,M,K);
        G_10           = zeros(N,M,K);
        G_11           = zeros(N,M,K);
        G_aa           = zeros(N*M,2,K);
        G_bb           = zeros(N*M,2,K);
%          p_pilot= 10^(5)*ones(M,K);
%         [SINR_TH_UatF] = THH_UATF(p,UE,AP,N ,R_G_HAT);
%        [R_G_HAT,C_G_HAT,R_G_HAT_concat,C_G_HAT_concat] = chan_estimate_1 (M,K,N,R,p_pilot);
         p_pilot=P_P_UP*ones(M,K);%10^(P_P_UP_dB_PLOT(n)/10)*ones(M,K);%10^(1)*ones(M,K);
         [R_G_HAT_concat,C_G_HAT_concat,R_G,Wd_0,Wd_1,WdC_0,WdC_1,R_G_HAT,C_G_HAT] = estimate(M,K,N,R,p_pilot,pol,alpha,kappa,alpha_ADC,DA);
         
         div=1;%[1,1/2,1/3,1/4,1/5,1/6,1/7,1/8,1/9,1/10];
         for i=1:length(div)
           p =  P_DATA/(div(i)*K);
           p1=  P_DATA*ones(M,K)/(div(i)*K);
%            p_pilot=10^(1)*ones(M,K);%10^(P_P_UP_dB_PLOT(n)/10)*ones(M,K);%10^(1)*ones(M,K);
%            tau_p=K;
%            [R_G_HAT_concat,C_G_HAT_concat,R_G] = estimate(M,K,N,R,p_pilot,pol,alpha);
          a=[1,0];
%           for v=1:2
%               DA=a(v);
           [SINR_TH_DE(i,n,ue_iters),M1,SIM_add_INT,c_mk_n0,c_mk_n1,c_mk_d0,c_mk_d1,combiner_trace,b_mk] = THH_UATF(p1,p1,K,M,N ,R_G_HAT_concat,pol,C_G_HAT_concat,kappa,alpha_ADC,WdC_0,WdC_1,p_pilot,DA,beta,R_G_HAT,C_G_HAT);%THH_UATF(p,K,M,N ,R_G_HAT_concat,pol,C_G_HAT_concat,kappa,alpha_ADC,WdC_0,WdC_1,p_pilot,DA,tau_b,B,decod_err_prob);
%            [SINR_TH_DE_opt(n,ue_iters),SINR_TH_DE(n,ue_iters),error1,p_user,p_opt(:,n),M1]=poweroptimize2(p,K,M,N ,R_G_HAT_concat,pol,C_G_HAT_concat,kappa,alpha_ADC,WdC_0,WdC_1,p_pilot,DA,beta);
%           end
%            [b_mk_opt(:,:,n),b_mk_init(:,:,n),rate_opt(n),rate0(n),GEE0(n),GEE_opt(n)]=poweroptimize2(p,K,M,beta,b_mk,c_mk_n0,c_mk_n1,c_mk_d0,c_mk_d1,T_C,combiner_trace,N,B,p_pilot,cvx);
%           gain(n)=SINR_TH_DE(i,n,ue_iters,1)/SINR_TH_DE(i,n,ue_iters,2);
           a1=sqrt(K)*(randn(N,1)+1i*randn(N,1))*1/sqrt(2);
         end
%            V_sq= 
        for ch = 1:CHAN_ITER
%             [G,H, G_LAP]    = FD_random_channel_generator(M,K,N, R_G,Gmean, R_H, Hmean,R_LR,R_LT,SIGMA_LAP); %Actual channel generated
           for ap =1:M
               for ue=1:K 
                G_00(:,ap,ue) = sqrt(1-alpha)*(sqrtm(R(:,:,ap,ue))*1/sqrt(2)*(randn(N,1)+1i*randn(N,1)));
                G_01(:,ap,ue)= sqrt(alpha)*(sqrtm(R(:,:,ap,ue))*1/sqrt(2)*(randn(N,1)+1i*randn(N,1)) ); %+g_LOS(:,ap,ue)
                G_10(:,ap,ue) = sqrt(alpha)*(sqrtm(R(:,:,ap,ue))*1/sqrt(2)*(randn(N,1)+1i*randn(N,1)) );
                G_11(:,ap,ue) = sqrt(1-alpha)*(sqrtm(R(:,:,ap,ue))*1/sqrt(2)*(randn(N,1)+1i*randn(N,1)));
               end
           end
                if pol==2
                   G=zeros(2*N,2,M,K);
                 for ue=1:K
                     for ap=1:M
                       G(:,:,ap,ue) =[G_00(:,ap,ue) G_01(:,ap,ue);G_10(:,ap,ue) G_11(:,ap,ue)];
                       
                     end
                 end
                 
                 
                else
                    G=zeros(N,M,K);
                    for ap=1:M
                     for ue=1:K
                        
                       G(:,ap,ue) = G_00(:,ap,ue)/sqrt(1-alpha);
                     end
                    end
                    G4= reshape(G,[N*M,K]);
               end
%                end
%            end
%            p =  P_DATA;
%            p1=  P_DATA*ones(M,K);
%            p_pilot= 10^(1)*ones(M,K);
%            tau_p=K;
           G_1 = G_00/sqrt(1-alpha); % chanel for 
        if ch==1
           for ch1=1:CHAN_ITER
            if pol==2
             Norm_factor=zeros(2,K,CHAN_ITER);
             [Norm_factor(:,:,ch1)] =normalization (M,K,N,R,alpha,p_pilot,pol,p1,C_G_HAT_concat,kappa,alpha_ADC,WdC_0,WdC_1);
            else
               Norm_factor=zeros(K,CHAN_ITER); 
              [Norm_factor(:,ch1)] =normalization (M,K,N,R,alpha,p_pilot,pol,p1,C_G_HAT_concat,kappa,alpha_ADC,WdC_0,WdC_1);
            end
           end
           if pol==2
            Norm_factor_true = mean(Norm_factor,3);
           else
            Norm_factor_true = mean(Norm_factor,2);
           end
        end
%         for b3=1:100
         if pft_csi==0
           [G_hat(:,:,:,:,ch),mmse_combiner,mrc_combiner,G_HAT,H] =chan_estimate (M,K,N,G,R,alpha,p_pilot,pol,p1,M1,R_G_HAT_concat,C_G_HAT_concat,a1,Norm_factor_true,kappa,alpha_ADC,Wd_0,Wd_1,WdC_0,WdC_1,DA);
        
         else
             G_hat=G;
             R_G_HAT =R;
             Q1=zeros(N,N,M);
             Q2=zeros(N,N,M);
         end
%         end

  
            if mmse==1
                combiner =mmse_combiner;
            elseif mrc ==1
                
                  combiner = mrc_combiner;
            
            else
                  combiner = G; 
           
           end
%            R_G1=R_G(:,:,:,1);
%            R_G2=R_G(:,:,:,2);
%            for ue=1:K
%             G11(:,ue)=R_G1(:,:,ue)*(randn(2*N*M,1) +1i*randn(2*N*M,1))*(1/sqrt(2));
%             G22(:,ue)=R_G2(:,:,ue)*(randn(2*N*M,1) +1i*randn(2*N*M,1))*(1/sqrt(2));
%            
%                G3(:,:,ue)=[G11(:,ue) G22(:,ue)];
%                a(ue,ch) = G11(:,ue)'*G11(:,ue) +G22(:,ue)'*G22(:,ue);%trace(R_G1(:,:,ue)) + trace(R_G2(:,:,ue));%trace(G3(:,:,ue)*G3(:,:,ue)');
%                a1(ue,ch) = H(:,1,ue)'*H(:,1,ue) +H(:,2,ue)'*H(:,2,ue);%trace(H(:,:,ue)*H(:,:,ue)');
%            end

            %This function generates original channels using the R and Gmean
           dot_prod(ch)=0;
            for ue=1:K
%                 dot_prod(ch) =(G_HAT(:,ue)'*(G4(:,ue)-G_HAT(:,ue)))/(K*N*M) +dot_prod(ch);
            end
            
            
            
                    
                      if pol==1
                       [SIM_SE(n,ch),V_sq(1:K,n,ch),Delta_sq(1:K,n,ch),V(1:K,n,ch),Delta(1:K,n,ch),add_INT(1:K,n,ch)] = SIM_SE_cal(G,p,combiner,K,M,N,pol,R,alpha,R_G,G_HAT,kappa,alpha_ADC,SIM_add_INT);
                      else
                       [SIM_SE(n,ch),V_sq(:,:,1:K,n,ch),Delta_sq(:,:,1:K,n,ch),V(:,:,1:K,n,ch),Delta(:,:,1:K,n,ch),add_INT(:,:,1:K,n,ch)] = SIM_SE_cal(G,p,combiner,K,M,N,pol,R,alpha,R_G,G_HAT,kappa,alpha_ADC,SIM_add_INT);
                      end
%            end
                    % This function returns the monte-carlo SINR of all the users
               
                
                
            
           
            
                
                %Calculation of Downlink NMSE Monte-carlo
            if DOWNLINK_NMSE==1
                for k_1=1:K
                    for k_2=1:K
                        DOWN_NMSE_SIM(ch) = DOWN_NMSE_SIM(ch) + abs(a_hat(k_1,k_2)-a(k_1,k_2))^2/(abs(a(k_1,k_2))^2*K); %Calculated NMSE over each channel iteration
                    end
                end
                %Calculation of montecarlo rate
            elseif RATE==1
               
                    for k=1:K
                        SE_MC(ch,k)=pre_log_fact*(TAU_D/(T_C))*log2(1+SINR_MC(ch,k)); %Evaluate Ergodic SE for each user and each channel iteratioin
                    end
               
            end
%         a(ch)=abs(V(1,1,ch)-V_sq(1,1,ch))/V(1,1,ch)*100;   
        end
%        m=(mean(a,2)-mean(a1,2))./(mean(a,2));
%        dot=mean(dot_prod);
%         a1=mean(a);
        for ue=1:K
            if pol==1
             SINR_UATF(n,ue) =(1-K/T_C)*log2(1 +abs(mean(V_sq(ue,n,:)))^2/(1+mean(Delta(ue,n,:))+mean(V(ue,n,:))-abs(mean(V_sq(ue,n,:)))^2+mean(add_INT(ue,n,:))));%1 +abs(mean(V_sq(ue,n,:)))^2/(1+mean(Delta(ue,n,:))+mean(V(ue,n,:))-abs(mean(V_sq(ue,n,:)))^2)
            elseif pol==2
              SIN(:,:,:,n,ue_iters) = UATF_SIM(V_sq,Delta,V,add_INT,n,K);
             SINR_UATF(n,ue,ue_iters) =(1-2*K/T_C)*log2(det(eye(2)+ SIN(:,:,ue,n,ue_iters)));%-qfuncinv(decod_err_prob)*log2(exp(1))*sqrt(2*SIN(1,1,ue,n,ue_iters)./(1+SIN(1,1,ue,n,ue_iters))./(tau_b*B))-qfuncinv(decod_err_prob)*log2(exp(1))*sqrt(2*SIN(2,2,ue,n,ue_iters)./(1+SIN(2,2,ue,n,ue_iters))./(tau_b*B));
            end
        end
        SINR_TH_UATF1(n) = sum(SINR_UATF(n,:)); %  end 
%         SINR_TH_UATF2(n1,ue_iters) = SINR_TH_UATF1(n);
%         % Averaging over Channel Iterations
%         
%            
%         if DOWNLINK_NMSE==1
%             DOWN_NMSE_MC(n)=sum(DOWN_NMSE_SIM)/CHAN_ITER;  % Evaluate Ergodic Downlik NMSE for each P_p/tau_p
%         elseif RATE==1
% %             SUM_SE_MC(n)=sum(sum(SE_MC))/CHAN_ITER; %Evaluate Ergodic SE
%             SUM_SE_MC(n)=(sum(SIM_SE(n,:)))/CHAN_ITER;
%         end
        
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
%                 plot(N_PLOT(1:n),SUM_SE_MC(1:n),'-hr');
                hold on
%                 plot(N_PLOT(1:n),SINR_TH_DE(1:n),'-^g');
                 plot(N_PLOT(1:n),gain(1:n),'-^g');
                hold on
%                 plot(N_PLOT(1:n),SINR_TH_UATF1(1:n),'-hb');
                hold on
                drawnow
                legend('simulated','lwr bnd')
                title(['UE SETUP =',num2str(ue_iters)])
            elseif M_AP==1
%                 plot(M_PLOT(1:n),SUM_SE_MC(1:n),'-hr');
                hold on
                plot(M_PLOT(1:n),SINR_TH_DE(1:n),'-^g');
%                 plot(M_PLOT(1:n),gain(1:n),'-^g');
                hold on
%                 plot(M_PLOT(1:n),SINR_TH_UATF1(1:n),'-hb');
%                 plot(M_PLOT(1:n),SUM_SE_TH(1:n),'-^b');
                hold on
                drawnow
                legend('simulated','lwr bnd')
                title(['UE SETUP =',num2str(ue_iters)])
            elseif K_USER==1
%                 plot(K_PLOT(1:n),SUM_SE_MC(1:n),'-hr');
                hold on
%                 plot(K_PLOT(1:n),SUM_SE_TH(1:n),'-^b');
                plot(K_PLOT(1:n),SINR_TH_DE(1:n),'-^g');
                hold on
                plot(K_PLOT(1:n),SINR_TH_UATF1(1:n),'-hb');
                hold on
                drawnow
                legend('simulated','lwr bnd')
                title(['UE SETUP =',num2str(ue_iters)])
            elseif TX_SNR==1
%                 plot(P_DATA_dB_PLOT(1:n),SUM_SE_MC(1:n),'-hr');
                hold on
%                 plot(P_DATA_dB_PLOT(1:n),SINR_TH_DE(1:n),'-^g');
                plot(P_DATA_dB_PLOT(1:n),GEE0(1:n),'-^g');
                hold on
                plot(P_DATA_dB_PLOT(1:n),GEE_opt(1:n),'-hb');
%                 plot(P_DATA_dB_PLOT(1:n),SINR_TH_DE_opt(1:n),'-^b');
                hold on
                drawnow
                legend('simulated','lwr bnd')
                title(['UE SETUP =',num2str(ue_iters)])
            elseif alpha1==1
%                 plot(P_DATA_dB_PLOT(1:n),SUM_SE_MC(1:n),'-hr');
                hold on
%                 plot(alpha_PLOT(1:n),SINR_TH_DE(1:n),'-^g');
                plot(alpha_PLOT(1:n),SINR_TH_UATF1(1:n),'-hb');
                hold on
                drawnow
                legend('simulated','lwr bnd')
                title(['UE SETUP =',num2str(ue_iters)]) 
            elseif kappa_value==1
%                 plot(P_DATA_dB_PLOT(1:n),SUM_SE_MC(1:n),'-hr');
                hold on
                plot(kappa_PLOT(1:n),SINR_TH_DE(1:n),'-^g');
%                 plot(kappa_PLOT(1:n),SINR_TH_UATF1(1:n),'-hb');
                hold on
                drawnow
                legend('simulated','lwr bnd')
                title(['UE SETUP =',num2str(ue_iters)])
            elseif alpha_ADC_value==1
%                 plot(P_DATA_dB_PLOT(1:n),SUM_SE_MC(1:n),'-hr');
                hold on
                plot(alpha_ADC_PLOT(1:n),SINR_TH_DE(1:n),'-^g');
%                 plot(alpha_PLOT(1:n),SINR_TH_UATF1(1:n),'-hb');
                hold on
                drawnow
                legend('simulated','lwr bnd')
                title(['UE SETUP =',num2str(ue_iters)])
            elseif ASD_value==1
%                 plot(P_DATA_dB_PLOT(1:n),SUM_SE_MC(1:n),'-hr');
                hold on
%                 plot(ASD_PLOT(1:n),SINR_TH_DE(1:n),'-^g');
                plot(ASD_PLOT(1:n),SINR_TH_UATF1(1:n),'-hb');
                hold on
                drawnow
                legend('simulated','lwr bnd')
                title(['UE SETUP =',num2str(ue_iters)])  
            elseif PILOT_POWER==1
                %plot(P_P_UP_dB_PLOT(1:n),SUM_SE_MC(1:n),'-hr');
                hold on
                plot(P_P_UP_dB_PLOT(1:n),SINR_TH_UATF1(1:n),'-^b');
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
% ecdf(real(reshape(SINR_TH_DE,1,[])));
% end
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