close all;
clear;
tStart = tic;
aID = getenv('SLURM_ARRAY_TASK_ID');

% This is for running on a cluster in parallel
% the bash script should give the aID as input
if (isempty(aID))
    warning('aID is empty. Trying SLURM ID.')
    aID = getenv('SLURM_ARRAY_TASK_ID');
end
if(isempty(aID))
    warning('aID is empty. Replacing it with 0010.')
    aID = '0022';
end
%RNG seed.
rng(str2double(aID),'twister');

%% GUE channel parameters
params.K_Factor = 9;         %dB -- %rician factor Ground UE  % if beta_gains=1
params.RAYLEIGH=0;   %1= rayleigh, % 0=rician
params.Perf_CSI =0;
params.cov_area = 1; %0.25; % 4; %km
%%
params.TAU_P_K_by_two = 0; %1;  
params.CH_estimation = 1;  % 1= have channel estimation
%%
params.LB=1;  %Lower bound
params.UB =1;  %Upper bound
params.no_of_rea = 10;     % no.of channel realizations
%%
% snr_db = -50:10:40;
params.snr_db = 30;
params.snr_db_mmw = 70;
params.ASD_VALUE = 0.25;%[0,0.25,0.5,0.75,1];  % [0,30,10]; %
params.ASD_CORR = 0;
params.Kt_Kr_vsUE  = 0; %0.175^2; %0.175^2; %0.175^2; %[1,2,3,4];  %to save 1=AP 0.1,UE=0.1;  2=AP 0.1,UE=0.3;  3=AP 0.3,UE=0.1

params.pilot_pow = 100;  % 0.1W   % UL pilot. power (W)
params.noiseFigure = 9; % gue
params.sigma_sf =4;
params.Band = 100e6;%20e6; %Communication bandwidth
params.tau_c = 200;      % coherence block length
% rng(2,'twister');
%%
% load('params.mat')
params.simTime = 2*60; %sec Total Simulation time should be more than 100.
%% Room Setup, UE placement, UE height
% We are considering an outdoor scenario where the UE is located at the
% center and gNBs are distributed around the UE. We only need to consider
% the coverageRange amount of distance from the UE.
params.coverageRange = 100;
length_area = 2*params.coverageRange;   
width_area = 2*params.coverageRange;
height_transmitter = 5;
params.areaDimensions = [width_area, length_area, height_transmitter];


params.coverageRange_sub6 = 1000;
length_area_sub6 = 2*params.coverageRange_sub6;   
width_area_sub6 = 2*params.coverageRange_sub6;
height_transmitter_sub6 = 4;
params.areaDimensions_sub6 = [width_area_sub6, length_area_sub6, height_transmitter_sub6];

params.hr = 1.4; %height receiver (UE), approximately the height a human holds the phone
params.ht = height_transmitter; %height transmitter (BS)
params.ht_sub6 = height_transmitter_sub6; %height transmitter (BS)
% params.r_min = rmin*rand(params.numUE,1);
% lambda_BS = [200,300,400,500]; %densityBS
% lambda_BS =[200,300]; %densityBS
% lambda_BS = 50:50:200;
lambda_BS = 10; %25;
% num_BS_arr = [2,5,10,20]; %densityBS
% numUE_sub6_arr = 2:2:10;
% numUE_sub6_arr = 10;
lambda_UE = 10;
lambda_UE_sub6 = 10; %:10:50;
% lambda_UE_sub6 = lambda_BS./2;
% lambda_UE_sub6 = lambda_BS.*2;
% for idxnumUEsub6 = 1:length(numUE_sub6_arr)
for idxUEDensity = 1:length(lambda_UE_sub6)
    for idxBSDensity = 1:length(lambda_BS)
        %% gNB locations
        % n = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange/1000)^2);
        % while (n==0)
        %     n = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange/1000)^2);       
        % end
        % params.numGNB = n;
        % params.RgNB = params.coverageRange * sqrt(rand(params.numGNB,1)); %location of gNBs (distance from origin)
        % % params.RgNB = (2*params.coverageRange/3) * ones(params.numGNB,1); %location of gNBs (distance from origin)
        % params.angleGNB = 2*pi*rand(params.numGNB,1);%location of gNBs (angle from x-axis)
        % params.locationsBS = [params.RgNB.*cos(params.angleGNB), params.RgNB.*sin(params.angleGNB)];
        % n = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange_sub6/1000)^2);
        % while (n==0)
        %     n = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange_sub6/1000)^2);       
        % end
        params.numGNB = 0;
        params.RgNB = params.coverageRange * sqrt(rand(params.numGNB,1)); %location of gNBs (distance from origin)
        % params.RgNB = (2*params.coverageRange/3) * ones(params.numGNB,1); %location of gNBs (distance from origin)
        params.angleGNB = 2*pi*rand(params.numGNB,1);%location of gNBs (angle from x-axis)
        params.locationsBS = [params.RgNB.*cos(params.angleGNB), params.RgNB.*sin(params.angleGNB)];

        n = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange_sub6/1000)^2);
        while (n==0)
            n = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange_sub6/1000)^2);       
        end
        params.numGNB_sub6 = n;
        % params.numGNB_sub6 = 1;
        params.RgNB_sub6 = params.coverageRange_sub6 * sqrt(rand(params.numGNB_sub6 - params.numGNB,1)); %location of gNBs (distance from origin)
        % params.RgNB_sub6 = (2*params.coverageRange_sub6/3) * ones(params.numGNB_sub6,1); %location of gNBs (distance from origin)
        % params.RgNB_sub6 = params.coverageRange_sub6 * ones(params.numGNB_sub6,1); %location of gNBs (distance from origin)
        params.angleGNB_sub6 = 2*pi*rand(params.numGNB_sub6 - params.numGNB,1);%location of gNBs (angle from x-axis)
        params.locationsBS_sub6 = [params.RgNB_sub6.*cos(params.angleGNB_sub6), params.RgNB_sub6.*sin(params.angleGNB_sub6)];  
        % params.locationsBS_sub6 = [params.locationsBS_sub6(1, :); params.locationsBS_sub6(3:10, :)];
        % params.numGNB_sub6 = 9;
        params.num_antennas_per_gNB = 10;
        params.num_antennas_per_gNB_mmW = 10;
       %%UE locations


        %%UE location
        % n = poissrnd(lambda_UE*pi*(params.coverageRange_sub6/1000)^2);
        % while (n==0)
        %     n = poissrnd(lambda_UE*pi*(params.coverageRange_sub6/1000)^2);
        % end
        % params.numUE = n;
        params.numUE = 1;
        params.RUE = 0; %params.coverageRange * sqrt(rand(params.numUE,1)); %location of UEs (distance from origin)
        params.angleUE = 2*pi*rand(params.numUE,1);%location of UEs (angle from x-axis)
        params.UE_locations = [params.RUE.*cos(params.angleUE), params.RUE.*sin(params.angleUE)];
        rmin = 1e9;
        params.r_min = rmin*ones(params.numUE,1);  %stores min rate requirement for all mmWave users
        params.bw_alloc = zeros(params.numUE,1);
        params.num_antennas_per_UE_mmW = 4;
      
        % params.numUE = 1;
        % params.RUE = 0; %params.coverageRange * sqrt(rand(params.numUE,1)); %location of UEs (distance from origin)
        % params.angleUE = 2*pi*rand(params.numUE,1);%location of UEs (angle from x-axis)
        % params.UE_locations = [params.RUE.*cos(params.angleUE), params.RUE.*sin(params.angleUE)];
        % rmin = 1e9;
        % params.r_min = rmin*ones(params.numUE,1);  %stores min rate requirement for all mmWave user
        % params.bw_alloc = zeros(params.numUE,1);

        % params.numUE_sub6 = 10;
        % params.numUE_sub6 = numUE_sub6_arr(idxnumUEsub6);
        % params.numUE_sub6 = poissrnd(lambda_UE_sub6(idxUEDensity)*pi*(params.coverageRange/1000)^2);
        % params.numUE_sub6 = poissrnd(lambda_UE_sub6(idxBSDensity)*pi*(params.coverageRange/1000)^2);
        % params.RUE_sub6 = params.coverageRange * sqrt(rand(params.numUE_sub6,1)); %location of UEs (distance from origin)
        % params.angleUE_sub6 = 2*pi*rand(params.numUE_sub6,1);%location of UEs (angle from x-axis)
        % params.UE_locations_sub6 = [params.RUE_sub6.*cos(params.angleUE_sub6), params.RUE_sub6.*sin(params.angleUE_sub6)];        
        % rmin_sub6 = 1e6;
        % params.r_min_sub6 = rmin_sub6*ones(params.numUE_sub6,1);  %stores min rate requirement for all sub-6 users

        % params.numUE_sub6 = poissrnd(lambda_UE_sub6(idxBSDensity)*pi*(params.coverageRange_sub6/1000)^2);
        % params.numUE_sub6 = poissrnd(lambda_UE_sub6(idxUEDensity)*pi*(params.coverageRange_sub6/1000)^2);
        params.numUE_sub6 = 60;
        params.RUE_sub6 = params.coverageRange_sub6*sqrt(rand(params.numUE_sub6,1)); %location of UEs (distance from origin)
        params.angleUE_sub6 = 2*pi*rand(params.numUE_sub6,1);%location of UEs (angle from x-axis)
        params.UE_locations_sub6 = [params.RUE_sub6.*cos(params.angleUE_sub6), params.RUE_sub6.*sin(params.angleUE_sub6)];        
        rmin_sub6 = 1e7;
        params.r_min_sub6 = rmin_sub6*ones(params.numUE_sub6,1);  %stores min rate requirement for all sub-6 users
        params.bw_alloc_sub6 = params.Band*ones(params.numUE_sub6,1);
        %% PHY layer params
        params.scs_mmw = 2e9;     %not using this parameter now
        params.scs_sub6 = 1e8;   %sub-6 GHz bandwidth 100 MHz
        params.num_sc_mmw = 1;    %not using this parameter now
        params.num_sc_sub6 = 1;   %sub-6 GHz considered as one full band
        
        %% UE angular coverage range (full 360 coverage for now)
        lookAngleCell{1} = [0,360];
        
        %% Blocker Properties and Simulation Duration
        params.lambdaBlockers = 0.01; %How many blockers around
        params.numBlockers = 4*(params.coverageRange)^2*params.lambdaBlockers;
        params.V = 1; %velocity of blocker m/s
        % 160-190 cm truncated gaussian with mean at 3 sigma to each sides.
        % params.hb = (175 + TruncatedGaussian(5, [-15,15], [1 params.numBlockers])) / 100;
        params.hb = 1.8*ones(1,params.numBlockers); %height blocker
        params.mu = 2; %Expected bloc dur =1/mu sec
        
        
        %% Protocol Characteristics    
        % gNB discovery time. This should include the initial beamsearch delays to
        % establish a viable channel btw UE and a recently unblocked gNB. This can
        % run in paralel, i.e., even when UE is connected to another gNB, UE can
        % discover other gNB in the background.
        % protocolParams.discovery_time = [20 50]*10^(-3);
        protocolParams.discovery_time = 50*10^(-3);
        
        % in ms, measurement report trigger time. 
        % BeamFailureMaxCount*MeasurementFrequency (10*2 = 20).
        % How long for the measurements to trigger before the UE starts the procedure to change
        % its gNB after it got blocked. 
        % protocolParams.FailureDetectionTime = [2 10 20]*10^(-3); 
        protocolParams.FailureDetectionTime = 1e-8; %20*10^(-3); 
        % Frame delivery delay, for FR2 small subframe duration around 1 ms
        protocolParams.frameDeliveryDelay = 0; %1*10^(-3); 
        % in numbers, how many hops until the MCG Failure info delivered and HO initated
        protocolParams.frameHopCount = 0; %6;
        % RACH delay, how long it takes to setup a connection with a discovered
        % gNB.
        % protocolParams.connection_time = [10 20 50]*10^(-3);
        protocolParams.connection_time = 50*10^(-3);
        
        % Signaling delay to start data transfer after rach is completed. Singaling
        % to setup UE data plane path Can be modeled as addition to the RACH in our
        % previous rotation simulations. Say Instead of RACH delay =20ms we have 
        % real RACH delay 20ms plus a extra delay coming from signaling around 10ms
        % protocolParams.signalingAfterRachTime = [5 10 20]*10^(-3);
        protocolParams.signalingAfterRachTime = 0; %20*10^(-3);
        
        
        %% protocol params from other paper
        frac = (mean(params.hb)-params.hr)/(params.ht-params.hr);
        protocolParams.theta = 2*params.V.*params.lambdaBlockers*frac/pi;
        % density_limits = [30,40,50,60,70];      % up to how many BS in coverage area, it will be ok since poisson distr.
        % K_list = [1,2,3,4];                      % Degree of Connectivity
        % protocolParams.omega_list = 1./protocolParams.connection_time;
        % self_blockage = 5/6;
       
        % N = params.num_antennas_per_gNB;  % antennas per AP
        N_mmW = params.num_antennas_per_gNB_mmW;  % antennas per AP
        N = params.num_antennas_per_gNB;  % antennas per AP
        N_UE = params.num_antennas_per_UE_mmW;
        L = params.numGNB_sub6;
        K = params.numUE + params.numUE_sub6;  % --Ground UEs
        K_mmW = params.numUE;
        snr_db = params.snr_db;
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
        % [channelGain_GUE,R_GUE,h_LOS_GUE,K_Rician,PLOS_GUE] = channel_cellfree_GUE3(K,L,N,ASD_VALUE,ASD_CORR,RAYLEIGH,0,K_Factor,cov_area,Band, [params.locationsBS; params.locationsBS_sub6], [params.UE_locations; params.UE_locations_sub6]);
        [channelGain_GUE_mmW,R_GUE_mmW,h_LOS_GUE_mmW,K_Rician_mmW,probLOS_mmW] = channel_cellfree_GUE3_mmW_only(K,K_mmW,L,N,N_mmW,N_UE,ASD_VALUE,ASD_CORR,RAYLEIGH,0,K_Factor,cov_area,Band, [params.locationsBS; params.locationsBS_sub6], [params.UE_locations; params.UE_locations_sub6]);
        [channelGain_GUE,R_GUE,h_LOS_GUE,~,~,K_Rician,probLOS] = channel_cellfree_GUE3(K,K_mmW,L,N,N_mmW,ASD_VALUE,ASD_CORR,RAYLEIGH,0,K_Factor,cov_area,Band, [params.locationsBS; params.locationsBS_sub6], [params.UE_locations; params.UE_locations_sub6]);
        params.BETA = channelGain_GUE';
        params.BETA_mmW = channelGain_GUE_mmW';
        params.ricianFactor = K_Rician';
        params.ricianFactor_mmW = K_Rician_mmW';
        params.R_GUE = R_GUE;
        params.h_LOS_GUE = h_LOS_GUE;
        params.R_GUE_mmW = R_GUE_mmW;
        params.h_LOS_GUE_mmW = h_LOS_GUE_mmW;  
        %%
        % UE states
        UE.numGNB = params.numGNB;
        numBS = params.numGNB;
        numUE = params.numUE;
        numUE_sub6 = params.numUE;

        %General timers for primary and secondary
        % UE.RACH_eff = RACH_eff;
        % UE.failureDetectionDelay = failureDetectionDelay;
        % UE.beamFailureRecoveryTimer = beamFailureRecoveryTimer;
        
        % Recording variables
        UE.primaryConnectionStarts = [];
        UE.primaryConnectionStartIndices = [];
        UE.primaryConnectionEnds = [];
        UE.primaryConnectionEndIndices = [];
        UE.primaryBSHistory = [];
        UE.primaryEventTimes = [];
        UE.primaryEventIndices = [];
        UE.primaryEventDescriptions = [];
        UE.primaryConnectionStateHistory = [];            
        
        UE.secondaryConnectionStarts = [];
        UE.secondaryConnectionStartIndices = [];
        UE.secondaryConnectionEnds = [];
        UE.secondaryConnectionEndIndices = [];
        UE.secondaryBSHistory = [];
        UE.secondaryEventTimes = [];
        UE.secondaryEventIndices = [];
        UE.secondaryEventDescriptions = [];
        UE.secondaryConnectionStateHistory = [];
        
        UE.sub6ConnectionStarts = [];
        UE.sub6ConnectionStartIndices = [];
        UE.sub6ConnectionEnds = [];
        UE.sub6ConnectionEndIndices = [];
        UE.sub6EventTimes = [];
        UE.sub6EventIndices = [];
        UE.sub6EventDescriptions = [];
        UE.sub6ConnectionStateHistory = [];
        
        %Primary State variables
        % UE.primaryConnectionState = 0;
        UE.primaryConnectionState = zeros(numUE,1);
        % UE.primaryConnectionStateHistory = [UE.primaryConnectionStateHistory; UE.primaryConnectionState];
        UE.primaryConnectionStateHistory = [UE.primaryConnectionStateHistory, UE.primaryConnectionState];
        % UE.primaryBSIdx = [];
        % UE.primaryTargetIdx = [];
        UE.primaryBSIdx = zeros(numUE,1);
        UE.primaryTargetIdx = zeros(numUE,1);
        % UE.primaryNextEventTime = -100;
        UE.primaryNextEventTime = -100*ones(numUE,1);
        
        %Secondary State variables
        % UE.secondaryConnectionState = 0;
        UE.secondaryConnectionState = zeros(numUE,1);
        % UE.secondaryConnectionStateHistory = [UE.secondaryConnectionStateHistory; UE.secondaryConnectionState];
        UE.secondaryConnectionStateHistory = [UE.secondaryConnectionStateHistory, UE.secondaryConnectionState];
        % UE.secondaryBSIdx = [];
        % UE.secondaryTargetIdx = [];
        UE.secondaryBSIdx = zeros(numUE,1);
        UE.secondaryTargetIdx = zeros(numUE,1);
        % UE.secondaryNextEventTime = -100;
        UE.secondaryNextEventTime = -100*ones(numUE,1);
        UE.tmpMCGBSIdx = zeros(numUE,1);
        
        %Sub 6 State variables
        % UE.sub6ConnectionState = 0;
        UE.sub6ConnectionState = zeros(numUE,1);
        % UE.sub6ConnectionStateHistory = [UE.sub6ConnectionStateHistory; UE.sub6ConnectionState];
        UE.sub6ConnectionStateHistory = [UE.sub6ConnectionStateHistory, UE.sub6ConnectionState];
        % UE.sub6NextEventTime = -100;
        UE.sub6NextEventTime = -100*ones(numUE,1);
        %%
        %offloading
        sub6ConnectionState = UE.sub6ConnectionState;
        ue_idx=1;
        sub6ConnectionState(ue_idx) = 1;
        r_calc_mmw = rate_analyticalv5_mmW_only(params, sub6ConnectionState); %= compute_link_rates_w_rician(params, link, ue_idx, UE.sub6ConnectionState);
        sub6ConnectionState(ue_idx) = 0;
        r_calc_sub6 = rate_analyticalv4(params, sub6ConnectionState); 
        rates_on_sub6_handoff = zeros(numUE+numUE_sub6,1);  %[r_min;r_min_sub6]; %r_min.*ones(numUE+numUE_sub6,1);
        for ue_idx_2 = 1:numUE
            if ((UE.sub6ConnectionState(ue_idx_2) == 1) || ue_idx_2 == ue_idx)
                % rates_on_sub6_handoff(ue_idx_2) = r_calc_sub6(ue_idx_2);
                rates_on_sub6_handoff(ue_idx_2) = r_calc_mmw(ue_idx_2);
            end
        end
        % for ue_idx_2 = 1+numUE:numUE+numUE_sub6
        for ue_idx_2 = 1:numUE_sub6
            % rates_on_sub6_handoff(ue_idx_2) = r_calc_sub6(ue_idx_2);
            rates_on_sub6_handoff(ue_idx_2+numUE) = r_calc_sub6(ue_idx_2);
        end
    end
end
tEnd = toc(tStart);
fprintf('Total runtime: %f seconds\n',tEnd)