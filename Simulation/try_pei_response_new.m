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
% aID = randi([0, 99]);
% rng(aID,'twister');
%% GUE channel parameters
params.K_Factor = 9;         %dB -- %rician factor Ground UE  % if beta_gains=1
params.RAYLEIGH=0;   %1= rayleigh, % 0=rician
params.Perf_CSI =1;
params.cov_area = 1; %0.25; % 4; %km
%%
params.TAU_P_K_by_two = 0; %1;  
params.CH_estimation = 0;  % 1= have channel estimation
%%
params.LB=1;  %Lower bound
params.UB =1;  %Upper bound
params.no_of_rea =2;     % no.of channel realizations
%%
% snr_db = -50:10:40;
params.snr_db = 40;
params.ASD_VALUE = 0;%[0,0.25,0.5,0.75,1];  % [0,30,10]; %
params.ASD_CORR = 0;
params.Kt_Kr_vsUE  = 0; %0.175^2; %0.175^2; %[1,2,3,4];  %to save 1=AP 0.1,UE=0.1;  2=AP 0.1,UE=0.3;  3=AP 0.3,UE=0.1

params.pilot_pow = 100;  % 0.1W   % UL pilot. power (W)
params.noiseFigure = 9; % gue
params.sigma_sf =4;
params.Band = 20e6; %Communication bandwidth


%% Define simulation setup

%Angular standard deviation in the local scattering model (in radians)
params.ASD_varphi = 0; %deg2rad(30); %azimuth angle
params.ASD_theta = 0; %deg2rad(15);  %elevation angle

%Total uplink transmit power per UE (mW)
params.p = 100;

params.num_antennas_per_gNB = 64;

%Total downlink transmit power per AP (mW)
params.rho_tot = 10^(3.6)*params.num_antennas_per_gNB; %200
% rho_tot_arr = [10:10:100, 200:100:1000, 2000:1000:10000];

%Power factor division
p_fac_arr = 10^2; %10.^(0:1:5);
% params.p_fac = 10;

%Prepare to save simulation results

% rng(2,'twister');
%%
% load('params.mat')
params.simTime = 10*60; %sec Total Simulation time should be more than 100.
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
height_transmitter_sub6 = 5;
params.areaDimensions_sub6 = [width_area_sub6, length_area_sub6, height_transmitter_sub6];
%%UE location
params.numUE = 1;
params.RUE = 0;  %params.coverageRange * sqrt(rand(params.numUE,1)); %location of UEs (distance from origin)
params.angleUE = 2*pi*rand(params.numUE,1);%location of UEs (angle from x-axis)
params.UE_locations = [params.RUE.*cos(params.angleUE), params.RUE.*sin(params.angleUE)];

params.hr = 1.4; %height receiver (UE), approximately the height a human holds the phone
params.ht = height_transmitter; %height transmitter (BS)
params.ht_sub6 = height_transmitter_sub6; %height transmitter (BS)
%Number of antennas per UE
% N_UE_mmW_arr = 2.^(0:1:5);
params.N_UE_mmW = 8;
params.N_UE_sub6 = 4;
rmin_arr = 4*10^8;
% params.r_min = rmin*rand(params.numUE,1);
% lambda_BS = 50:50:200;%densityBS
lambda_BS = 25;
% num_BS_arr = [2,5,10,20]; %densityBS
% numUE_sub6_arr = 2:2:10;
% numUE_sub6_arr = 10;
lambda_UE_sub6 = [30:20:90, 100]; %:100:2000;
% for idxnumUEsub6 = 1:length(numUE_sub6_arr)
params.loss_pc_thresh = 10;
params.Lmax=4;
for idxBSDensity = 1:length(lambda_BS)
    %% gNB locations
    % params.numGNB = 10;
    n = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange/1000)^2);
    while (n==0)
        n = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange/1000)^2);       
    end
    params.numGNB = 1;
    %     params.numGNB = n;
    % params.numGNB = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange/1000)^2);
    % params.numGNB = floor(lambda_BS(idxBSDensity)*pi*(params.coverageRange/1000)^2);
    params.RgNB = params.coverageRange * sqrt(rand(params.numGNB,1)); %location of gNBs (distance from origin)
    % params.RgNB = (2*params.coverageRange/3) * ones(params.numGNB,1); %location of gNBs (distance from origin)
    params.angleGNB = 2*pi*rand(params.numGNB,1);%location of gNBs (angle from x-axis)
    params.locationsBS = [params.RgNB.*cos(params.angleGNB), params.RgNB.*sin(params.angleGNB)];
    
    n = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange_sub6/1000)^2);
    while (n<=params.numGNB) %(n==0)
        n = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange_sub6/1000)^2);       
    end
%     params.numGNB_sub6 = 1;
    params.numGNB_sub6 = n;
    params.RgNB_sub6 = params.coverageRange_sub6 * sqrt(rand(params.numGNB_sub6 - params.numGNB,1)); %location of gNBs (distance from origin)
    % params.RgNB = (2*params.coverageRange/3) * ones(params.numGNB,1); %location of gNBs (distance from origin)
    params.angleGNB_sub6 = 2*pi*rand(params.numGNB_sub6 - params.numGNB,1);%location of gNBs (angle from x-axis)
    params.locationsBS_sub6 = [params.RgNB_sub6.*cos(params.angleGNB_sub6), params.RgNB_sub6.*sin(params.angleGNB_sub6)];  
%     params.Lmax=n;
    for idxUEDensity = 1:length(lambda_UE_sub6)
        params.ue_rearranged = [];    
        %%UE locations
        n = poissrnd(lambda_UE_sub6(idxUEDensity)*pi*(params.coverageRange_sub6/1000)^2);
        while (n==0)
            n = poissrnd(lambda_UE_sub6(idxUEDensity)*pi*(params.coverageRange_sub6/1000)^2);       
        end
%         params.numUE_sub6 = 1;
        params.numUE_sub6 = n;
        params.RUE_sub6 = params.coverageRange_sub6*sqrt(rand(params.numUE_sub6,1)); %location of UEs (distance from origin)
        params.angleUE_sub6 = 2*pi*rand(params.numUE_sub6,1);%location of UEs (angle from x-axis)
        params.UE_locations_sub6 = [params.RUE_sub6.*cos(params.angleUE_sub6), params.RUE_sub6.*sin(params.angleUE_sub6)];   
        for idxrmin = 1:length(rmin_arr)
            rmin = rmin_arr(idxrmin);
            params.r_min = rmin*ones(params.numUE,1);  %stores min rate requirement for all mmWave users
            rmin_sub6 = 1e5;
            params.r_min_sub6 = rmin_sub6*ones(params.numUE_sub6,1);  %stores min rate requirement for all sub-6 users
            params.rate_reduce_threshold = 5e7;
            %Length of the coherence block
            params.tau_c = 200;
            
            %Compute number of pilots per coherence block
            params.tau_p = params.numUE+params.numUE_sub6;
            
            %Compute the prelog factor assuming only downlink data transmission
    %         params.preLogFactor = (params.tau_c-params.tau_p)/params.tau_c;
            params.preLogFactor = 1;
    
            %Number of setups with random UE locations
            params.nbrOfSetups = 100;
                    
                  
            %Number of channel realizations per setup
            params.nbrOfRealizations = 100;
            
            %% PHY layer params
            params.scs_mmw = 2e9;     %not using this parameter now
            params.scs_sub6 = 2e7;   %sub-6 GHz bandwidth 100 MHz
            params.num_sc_mmw = 1;    %not using this parameter now
            params.num_sc_sub6 = 1;   %sub-6 GHz considered as one full band
            
            %% UE angular coverage range (full 360 coverage for now)
            lookAngleCell{1} = [0,360];
            
            %% Blocker Properties and Simulation Duration
            params.lambdaBlockers = 0.01; %How many blockers around
            params.numBlockers = 4*(params.coverageRange)^2*params.lambdaBlockers;
    %         params.numBlockers = 4*(params.coverageRange_sub6)^2*params.lambdaBlockers;
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
            protocolParams.discovery_time = 20*10^(-3);
            
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
            
            for p_idx = 1:length(p_fac_arr)
                params.p_fac = p_fac_arr(p_idx);
                params.p_fac_rearrange = 0.1*p_fac_arr(p_idx);                
                %% protocol params from other paper
                frac = (mean(params.hb)-params.hr)/(params.ht-params.hr);
                protocolParams.theta = 2*params.V.*params.lambdaBlockers*frac/pi;
                % density_limits = [30,40,50,60,70];      % up to how many BS in coverage area, it will be ok since poisson distr.
                % K_list = [1,2,3,4];                      % Degree of Connectivity
                % protocolParams.omega_list = 1./protocolParams.connection_time;
                % self_blockage = 5/6;
               
                N = params.num_antennas_per_gNB;  % antennas per AP
                L = params.numGNB_sub6;
                K = params.numUE + params.numUE_sub6;  % --Ground UEs
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
                
%                 [gainOverNoisedB,R_gNB,R_ue_mmW,R_ue_sub6,pilotIndex,D,D_small,APpositions,UEpositions,distances] = generateSetup(params.numGNB,params.numGNB_sub6,params.numUE,params.numUE+params.numUE_sub6,params.num_antennas_per_gNB,params.N_UE_mmW,params.N_UE_sub6,params.coverageRange,params.coverageRange_sub6,params.tau_p,1,0,params.ASD_varphi,params.ASD_theta);
                [gainOverNoisedB,R_gNB,R_ue_mmW,R_ue_sub6,pilotIndex,D,D_small,APpositions,UEpositions,distances] = generateSetup(params,1,str2double(aID));
                params.BETA = db2pow(gainOverNoisedB);   
                params.D = D;
                params.R_gNB = R_gNB;
                params.R_ue_mmW = R_ue_mmW;
                params.R_ue_sub6 = R_ue_sub6;
                numUE = params.numUE;
                sub6ConnectionState = zeros(numUE,1);
                ap_idxs = find(D(:,1));
                ue_idxs = 1;
                for a = 1:length(ap_idxs)
                    ap_idx = ap_idxs(a);
                    ue_idxs = union(ue_idxs,find(D(ap_idx,:)));
                end
                [channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW] = computePhysicalChannels_sub6_MIMO(params);
                ue_idx = 1;
                rate_dl_before_handoff = compute_link_rates_MIMO(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,ue_idx,sub6ConnectionState);                                              
                sub6ConnectionState(ue_idx) = 1;
                [params.D, ue_idxs_affected] = AP_reassign(params,ue_idx);
                params.ue_rearranged = ue_idxs_affected;
                [channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW] = computePhysicalChannels_sub6_MIMO(params);
                rate_dl_after_handoff = compute_link_rates_MIMO(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,ue_idx,sub6ConnectionState);                                              
%                 rate_dl_after_handoff = compute_link_rates_MIMO_fdm(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,ue_idx,sub6ConnectionState);                                              
                numUE = params.numUE;
                numUE_sub6 = params.numUE_sub6;
                numBS = size(params.locationsBS,1);
                numBlockers = params.numBlockers;
                %% Recording the Results
        
                %Taking care of folder directory creation etc
                dataFolder = 'resultData';
                impactFolder = strcat(dataFolder,'/impactResults');
                if not(isfolder(impactFolder))
                    mkdir(impactFolder)
                end
                result_string = ['/handoff_impact_',num2str(numUE),...
                    'UE_',num2str(lambda_BS(idxBSDensity)),...
                    'lambdaBS_',num2str(lambda_UE_sub6(idxUEDensity)),...
                    'lambdaUE_',num2str(numBlockers), 'Blockers_randomHeight_', num2str(aID),'Pow_factor_', num2str(params.p_fac)];    
                %Since we are mostly interested in blockage probability, we want to
                %transfer the data quickly to our local machine from server. We will save
                %the results as a txt file.
                recording_text_file_string = strcat(impactFolder,result_string,'.csv');
                fileID = fopen(recording_text_file_string,'w');
                output_categories = ['UE idx,','lambdaBS,','lambdaUE,','numBlockers,',...
                    'powFactor,','rate_before_handoff,','rate_after_handoff,','rate_before_handoff_affected,','rate_after_handoff_affected,','rate_mmW\n'];  
                fprintf(fileID,output_categories);
                min_rate_req = params.r_min(1);   
                mean_rate_before_handoff = mean(rate_dl_before_handoff(2:end));
                mean_rate_after_handoff = mean(rate_dl_after_handoff(2:end)); 
                mean_rate_before_handoff_affected = mean(rate_dl_before_handoff(ue_idxs_affected));
                mean_rate_after_handoff_affected = mean(rate_dl_after_handoff(ue_idxs_affected)); 
                rate_mmW = rate_dl_after_handoff (1);
                formatSpec = '%d,%d,%d,%d,%d,%.16f,%.16f,%.16f,%.16f,%.16f\n';
                fprintf(fileID,formatSpec,ue_idx, lambda_BS(idxBSDensity),lambda_UE_sub6(idxUEDensity),numBlockers,...
                   params.p_fac,mean_rate_before_handoff,mean_rate_after_handoff,mean_rate_before_handoff_affected,mean_rate_after_handoff_affected,rate_mmW);                      
                fclose(fileID);
            end
        end
    end
end
tEnd = toc(tStart);
fprintf('Total runtime: %f seconds\n',tEnd)
% end