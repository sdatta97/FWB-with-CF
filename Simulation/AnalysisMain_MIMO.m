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
params.Perf_CSI =1;
params.cov_area = 1; %0.25; % 4; %km
%%
params.TAU_P_K_by_two = 0; %1;  
params.CH_estimation = 0;  % 1= have channel estimation
%%
params.LB=1;  %Lower bound
params.UB =1;  %Upper bound
params.no_of_rea = 1;     % no.of channel realizations
%%
% snr_db = -50:10:40;
params.snr_db = 40;
params.ASD_VALUE = 0;%[0,0.25,0.5,0.75,1];  % [0,30,10]; %
params.ASD_CORR = 0;
params.Kt_Kr_vsUE  = 0; %0.175^2; %0.175^2; %[1,2,3,4];  %to save 1=AP 0.1,UE=0.1;  2=AP 0.1,UE=0.3;  3=AP 0.3,UE=0.1

params.pilot_pow = 100;  % 0.1W   % UL pilot. power (W)
params.noiseFigure = 9; % gue
params.sigma_sf =4;
params.Band = 100e6; %Communication bandwidth


%% Define simulation setup

%Angular standard deviation in the local scattering model (in radians)
params.ASD_varphi = deg2rad(30); %azimuth angle
params.ASD_theta = 0; %deg2rad(15);  %elevation angle

%Total uplink transmit power per UE (mW)
params.p = 100;

%Total downlink transmit power per AP (mW)
% rho_tot_arr = [10:10:100, 200:100:1000, 2000:1000:10000];

%Power factor division
p_fac_arr = 10; %.^(1:1:2);
% params.p_fac = 10;
percent_fr2_UE_arr = 20; %15:5:20;

%Prepare to save simulation results

% rng(2,'twister');
%%
% load('params.mat')
params.simTime = 60*60; %sec Total Simulation time should be more than 100.
%% Room Setup, UE placement, UE height
% We are considering an outdoor scenario where the UE is located at the
% center and gNBs are distributed around the UE. We only need to consider
% the coverageRange amount of distance from the UE.
params.deployRange = 300; %20:20:100;
params.deployRange_sub6 = 1000;

params.coverageRange_sub6 = 430;
length_area_sub6 = 2*params.coverageRange_sub6;   
width_area_sub6 = 2*params.coverageRange_sub6;
% length_area_sub6 = 2*(params.deployRange + params.coverageRange_sub6);   
% width_area_sub6 = 2*(params.deployRange + params.coverageRange_sub6);
height_transmitter = 5;
height_transmitter_sub6 = 4; % 5;
params.areaDimensions_sub6 = [width_area_sub6, length_area_sub6, height_transmitter_sub6];

params.hr = 1.4; %height receiver (UE), approximately the height a human holds the phone
params.ht = height_transmitter; %height transmitter (BS)
params.ht_sub6 = height_transmitter_sub6; %height transmitter (BS)
params.num_antennas_per_gNB = 64;
params.rho_tot = 10^(3.6)*params.num_antennas_per_gNB; %200;

% params.num_antennas_per_gNB = 8;
%Number of antennas per UE
params.N_UE_mmW = 1; %8;
params.N_UE_sub6 = 1; %4;
rmin_arr = 4*10^8;
lambda_BS = 25; %([6 7 8 9 10]).^2;
lambda_UE_sub6 = 250; %250:250:1000; %200:10:250; %150; %100:50:200; %[30:20:90, 100]; %100;
params.Lmax = 4;
% for idxnumUEsub6 = 1:length(numUE_sub6_arr)
lb_thresh = 0.1; %[0.1, 0.25, 0.5]; [0:0.05:0.1 0.5 1];
for idxnumUE = 1:length(percent_fr2_UE_arr)
    % n = poissrnd((percent_fr2_UE_arr(idxnumUE)/100)*lambda_UE_sub6*pi*(params.deployRange/1000)^2);
    % while(n==0)
    %     n = poissrnd((percent_fr2_UE_arr(idxnumUE)/100)*lambda_UE_sub6*pi*(params.deployRange/1000)^2);
    % end
    % params.numUE = n;
    params.numUE = ceil((percent_fr2_UE_arr(idxnumUE)/100)*lambda_UE_sub6*pi*(params.deployRange/1000)^2);
    % params.numUE = 20;
    %%UE location
    deployRange = params.deployRange; %(idxdeployRange);
    params.RUE =  deployRange*sqrt(rand(params.numUE,1)); %location of UEs (distance from origin)
    params.angleUE = 2*pi*rand(params.numUE,1);%location of UEs (angle from x-axis)
    params.UE_locations = [params.RUE.*cos(params.angleUE), params.RUE.*sin(params.angleUE)];
    for idxBSDensity = 1:length(lambda_BS)
        %% gNB locations
        params.numGNB = ceil(lambda_BS(idxBSDensity)*(params.deployRange_sub6/1000)^2);
        params.xGNB = (-params.deployRange_sub6/2):(params.deployRange_sub6/(sqrt(params.numGNB)-1)):(params.deployRange_sub6/2);
        params.yGNB = (-params.deployRange_sub6/2):(params.deployRange_sub6/(sqrt(params.numGNB)-1)):(params.deployRange_sub6/2);
        % params.locationsBS = [params.RgNB.*cos(params.angleGNB), params.RgNB.*sin(params.angleGNB)];
        params.locationsBS = (combvec(params.xGNB,params.yGNB)).';
        params.coverageRange =  (params.deployRange_sub6/(sqrt(params.numGNB)-1))/sqrt(2);%125*sqrt(2); %100;
        length_area = 2*params.coverageRange;   
        width_area = 2*params.coverageRange;
        % length_area = 2*(params.deployRange + params.coverageRange);   
        % width_area = 2*(params.deployRange + params.coverageRange);
        params.areaDimensions = [width_area, length_area, height_transmitter];
                %Length of the coherence block
%                 params.tau_c = 200;
            
            %Compute number of pilots per coherence block
%                 params.tau_p = params.numUE+params.numUE_sub6;
            
            %Compute the prelog factor assuming only downlink data transmission
    %         params.preLogFactor = (params.tau_c-params.tau_p)/params.tau_c;
        params.preLogFactor = 1;

        %Number of setups with random UE locations
        params.nbrOfSetups = 100;
                
              
        %Number of channel realizations per setup
        params.nbrOfRealizations = 100;
        
        %% PHY layer params
        params.scs_mmw = 2e9;     %not using this parameter now
        params.scs_sub6 = [0.5*(params.Band), 0.5*(params.Band)];   %sub-6 GHz bandwidth 100 MHz
        params.num_sc_mmw = 1;    %not using this parameter now
        params.num_sc_sub6 = 2;   %sub-6 GHz considered as one full band                    
        %% UE angular coverage range (full 360 coverage for now)
        lookAngleCell{1} = [0,360];
        
        %% Blocker Properties and Simulation Duration
        params.lambdaBlockers = 0.01; %How many blockers around
        % params.numBlockers = 4*(params.coverageRange)^2*params.lambdaBlockers;
        params.numBlockers = floor(pi*(params.coverageRange)^2*params.lambdaBlockers);
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
        % protocolParams.discovery_time = [5 20]*10^(-3);
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
        
        
        %% protocol params from other paper
        frac = (mean(params.hb)-params.hr)/(params.ht-params.hr);
        protocolParams.theta = 2*params.V.*params.lambdaBlockers*frac/pi;
        % density_limits = [30,40,50,60,70];      % up to how many BS in coverage area, it will be ok since poisson distr.
        % K_list = [1,2,3,4];                      % Degree of Connectivity
        % protocolParams.omega_list = 1./protocolParams.connection_time;
        % self_blockage = 5/6;
        
        %% Mobile blockage events
        % tic
        numBS_mobile = [];
        dataBS_mobile = [];
        for i = 1:(params.numUE)
            % dataBS_mobile = [dataBS_mobile; computeBlockageEvents(params,i)];
            [dataBS, nBS] = computeBlockageEvents(params,i);
            dataBS_mobile = [dataBS_mobile; dataBS];
            numBS_mobile = [numBS_mobile; nBS];
        end

        fprintf('Blocker generation, physical blockage and channel computation done : %f seconds\n',toc)
        %% Simulation FR1 setup

        for idxUEDensity = 1:length(lambda_UE_sub6)
            params.ue_rearranged = [];        
            %%UE locations
            % n = poissrnd(lambda_UE_sub6(idxUEDensity)*pi*(params.coverageRange_sub6/1000)^2);
            % while (n==0)
            %     n = poissrnd(lambda_UE_sub6(idxUEDensity)*pi*(params.coverageRange_sub6/1000)^2);       
            % end
            n = ceil(lambda_UE_sub6(idxUEDensity)*pi*(params.coverageRange_sub6/1000)^2);
            params.numUE_sub6 = n;
            % params.RUE_sub6 = params.coverageRange_sub6*sqrt(rand(params.numUE_sub6,1)); %location of UEs (distance from origin)
            % params.RUE_sub6 = (norm(params.UE_locations(:,1) - params.UE_locations(:,2))/2+params.coverageRange_sub6)*sqrt(rand(params.numUE_sub6,1)); %location of UEs (distance from origin)
            params.RUE_sub6 = (max(sqrt(sum(params.UE_locations.^2,2)))+params.coverageRange_sub6)*sqrt(rand(params.numUE_sub6,1)); %location of UEs (distance from origin)
            params.angleUE_sub6 = 2*pi*rand(params.numUE_sub6,1);%location of UEs (angle from x-axis)
            % params.UE_locations_sub6 =  mean(params.UE_locations,1) + [params.RUE_sub6.*cos(params.angleUE_sub6), params.RUE_sub6.*sin(params.angleUE_sub6)];   
            params.UE_locations_sub6 =  [params.RUE_sub6.*cos(params.angleUE_sub6), params.RUE_sub6.*sin(params.angleUE_sub6)];   
            
            K = params.numUE + params.numUE_sub6;  % --Ground UEs
            snr_db = params.snr_db;
            ASD_VALUE = params.ASD_VALUE;
            ASD_CORR = params.ASD_CORR;
            Kt_Kr_vsUE = params.Kt_Kr_vsUE;
            K_Factor = params.K_Factor;
            RAYLEIGH=params.RAYLEIGH;   %1= rayleigh, % 0=rician
            Perf_CSI = params.Perf_CSI;
            cov_area = params.cov_area;
            %%
            no_of_rea = params.no_of_rea;     % no.of channel realizations
            %%
            pilot_pow = params.pilot_pow; 
            noiseFigure = params.noiseFigure;
            sigma_sf = params.sigma_sf;
            numBS = params.numGNB;
            numUE = params.numUE;
            Band = params.Band; %Communication bandwidth
            num_sc_sub6 = params.num_sc_sub6;
            % params.user_sc_alloc = randi([1,num_sc_sub6],K,1);
            params.user_sc_alloc = zeros(K,num_sc_sub6);
%                 tau_c = params.tau_c;      % coherence block length 
            [gainOverNoisedB,R_gNB,R_ue_mmW,R_ue_sub6,pilotIndex,D,D_small,APpositions,UEpositions,distances] = generateSetup(params,1,str2double(aID));
            num_sc_sub6 = params.num_sc_sub6;
            params.BETA = db2pow(gainOverNoisedB);   
            params.D = D;
            % params.D = D_small;
            params.R_gNB = R_gNB;
            params.R_ue_mmW = R_ue_mmW;
            params.R_ue_sub6 = R_ue_sub6;                       
            for idxrmin = 1:length(rmin_arr)
                for idx_p = 1:length(p_fac_arr)
                    for idxlbthres = 1:length(lb_thresh)
                        lb_thres = lb_thresh(idxlbthres);
                        rmin = rmin_arr(idxrmin);
                        rmin_sub6 = 35e6;
                        params.loss_pc_thresh = lb_thres;
                        params.r_min_sub6 = rmin_sub6;  %stores min rate requirement for all sub-6 users
                        params.r_min = rmin;  %stores min rate requirement for all mmWave users
                        params.rate_reduce_threshold = 5e7;
                        params.p_fac = p_fac_arr(idx_p);
                        params.p_fac_rearrange = 1; % 0.1*p_fac_arr(idx_p);         
                        outage_probability_analysis = zeros(params.numUE,length(protocolParams.discovery_time),length(protocolParams.connection_time));
                        outage_probability_analysis_wo_cf = zeros(params.numUE,length(protocolParams.discovery_time),length(protocolParams.connection_time));
                        outage_duration_analysis = zeros(params.numUE,length(protocolParams.discovery_time),length(protocolParams.connection_time));
                        outage_duration_analysis_wo_cf = zeros(params.numUE,length(protocolParams.discovery_time),length(protocolParams.connection_time));

                        for idxDiscDelay = 1:length(protocolParams.discovery_time)
                            for idxFailureDetectionDelay = 1:length(protocolParams.FailureDetectionTime)
                                for idxConnDelay = 1:length(protocolParams.connection_time)
                                    for idxSignalingAfterRachDelay = 1:length(protocolParams.signalingAfterRachTime)
                                        theta = protocolParams.theta;
                                        omega = 1/(protocolParams.connection_time(idxConnDelay) + protocolParams.signalingAfterRachTime(idxSignalingAfterRachDelay) + protocolParams.FailureDetectionTime(idxFailureDetectionDelay));
                                        psi = 1/(protocolParams.discovery_time(idxDiscDelay) + 1/params.mu);
                                        plos2 = pLoS2(params.locationsBS, [params.UE_locations;params.UE_locations_sub6], theta,omega,psi);
                                        plos = zeros(params.numUE,1);
                                        % [~,idx_max] = mink(plos2(:,1:params.numUE),2,1);
                                        % [~,idx_max] = mink(plos2(:,1:params.numUE),1,1);
                                        [~,idx_max] = mink(plos2(:,1:params.numUE),params.numGNB,1);
                                        for k = 1:params.numUE
                                            % plos(k) = prod(plos2(:,k),1);
                                            plos(k) = prod(plos2(idx_max(:,k),k),1);
                                        end
                                        % rate_dl = rate_analytical(params);
                                        for ue_idx = 1:numUE
                                           % plos3 = pLoS3(params.locationsBS, params.UE_locations(k,:), theta,omega,psi,idx_max);
                                           % plos4 = pLoS4(params.locationsBS, params.UE_locations(k,:), theta,omega,psi,idx_max);
                                           [pos3, tos3] = pLoS3(params.locationsBS, params.UE_locations(ue_idx,:), theta,omega,psi,idx_max);
                                           sub6ConnectionState = zeros(params.numUE,1);
                                           rate_dl_before_handoff = real(rate_analytical(params,sub6ConnectionState));
                                           [~, params.ue_rearranged] = AP_reassign(params,ue_idx);
                                           ues_not_affected = setdiff((1+params.numUE):(params.numUE+params.numUE_sub6),params.ue_rearranged);
                                           lb = quantile(rate_dl_before_handoff(params.ue_rearranged),params.loss_pc_thresh);
                                           bw_alloc = Band - params.r_min_sub6/lb;
                                           params.scs_sub6(1) = bw_alloc;
                                           params.scs_sub6(2) = Band - bw_alloc;
                                           sub6ConnectionState(ue_idx) = 1;
                                           user_sc_alloc = zeros(params.numUE+params.numUE_sub6,params.num_sc_sub6);
                                           user_sc_alloc(find(sub6ConnectionState),1) = 1;
                                           user_sc_alloc(find(sub6ConnectionState),2) = 0;
                                           user_sc_alloc(ues_not_affected,1) = 1;
                                           user_sc_alloc(ues_not_affected,2) = 1;
                                           user_sc_alloc(params.ue_rearranged,1) = 0;
                                           user_sc_alloc(params.ue_rearranged,2) = 1;
                                           params.user_sc_alloc = user_sc_alloc;
                                           ues_sharing = union(((1:numUE).*sub6ConnectionState),ues_not_affected);
                                           rate_dl_after_handoff = real(rate_analyticalv2(params,sub6ConnectionState));
                                           lb = quantile(rate_dl_after_handoff(ues_not_affected),params.loss_pc_thresh);
                                           if all(rate_dl_after_handoff(ue_idx) >= params.r_min) && (lb >= params.r_min_sub6)
                                               p_out = 0;
                                           else
                                               p_out = 1;
                                           end
                                           outage_probability_analysis(ue_idx,idxDiscDelay, idxFailureDetectionDelay, idxConnDelay, idxSignalingAfterRachDelay) = outage_probability_analysis(ue_idx,idxDiscDelay, idxConnDelay) + p_out*pos3;
                                           outage_probability_analysis_wo_cf(ue_idx,idxDiscDelay, idxFailureDetectionDelay, idxConnDelay, idxSignalingAfterRachDelay) = outage_probability_analysis_wo_cf(ue_idx,idxDiscDelay, idxConnDelay) + pos3;
                                           outage_duration_analysis(ue_idx,idxDiscDelay, idxFailureDetectionDelay, idxConnDelay, idxSignalingAfterRachDelay) = outage_probability_analysis(ue_idx,idxDiscDelay, idxConnDelay) + p_out*tos3;
                                           outage_duration_analysis_wo_cf(ue_idx,idxDiscDelay, idxFailureDetectionDelay, idxConnDelay, idxSignalingAfterRachDelay) = outage_probability_analysis_wo_cf(ue_idx,idxDiscDelay, idxConnDelay) + tos3;
                                        end
                                    end
                                end
                            end
                        end                
                        %% Create Discrete Time Event Simulation input
                        simInputs.params = params;
                        simInputs.protocolParams = protocolParams;
                        simOutputs = cell(length(protocolParams.discovery_time),...
                                          length(protocolParams.FailureDetectionTime),...
                                          length(protocolParams.connection_time),...
                                          length(protocolParams.signalingAfterRachTime));                        
                          
                        %% Recording the Results
                
                        %Taking care of folder directory creation etc
                        dataFolder = 'resultData';
                        analFolder = strcat(dataFolder,'/analysisResults');
                        if not(isfolder(dataFolder))
                            mkdir(dataFolder)
                        end
                        if not(isfolder(analFolder))
                            mkdir(analFolder)
                        end
                
                
                        %Saving all results as a structure
                        dataDescription = {'simOutputs is a 4D array';...
                            ', for mesh of params ordered as follows';...
                            'First Dimension: discovery_time';...
                            'Second Dimension: FailureDetectionTime';...
                            'Third Dimension: connection_time (RACH)';...
                            'Fourth Dimension: signalingAfterRachTime';...
                            '=================================';...
                            'Each element is a struct'};
                
                        numUE = params.numUE;
                        numUE_sub6 = params.numUE_sub6;
                        numBS = size(params.locationsBS,1);
                        numBlockers = params.numBlockers;
                        result_string = strcat('/results_',num2str(percent_fr2_UE_arr(idxnumUE)),...
                            'percentfr2UE_',num2str(lambda_BS(idxBSDensity)),...
                            'lambdaBS_',num2str(lambda_UE_sub6(idxUEDensity)),...
                            'lambdaUE_',num2str(deployRange),...
                            'deployRange_',num2str(numBlockers), 'Blockers_randomHeight_', num2str(aID),'Min_rate', num2str(rmin), "Pow_fac", num2str(params.p_fac), "lb_thres", num2str(100*params.loss_pc_thresh));
                        recording_text_file_string = strcat(analFolder,result_string,'.csv');
                        fileID = fopen(recording_text_file_string,'w');
                        output_categories = ['UE idx,','lambdaBS,','lambdaUE,','numBlockers,',...
                            'deployRange,','discoveryDelay,','failureDetectionDelay,','connectionSetupDelay,',...
                            'signalingAfterRachDelay,','frameHopCount,','frameDeliveryDelay,'...
                            'minRatereq,','powerFac,','lower_bound_thresh,', 'meanOutageDuration_wo_CF,','outageProbability_wo_CF,','meanOutageDuration,','outageProbability\n'];
                
                        fprintf(fileID,output_categories);
                
                        for idxDiscDelay = 1:length(protocolParams.discovery_time)
                            for idxFailureDetectionDelay = 1:length(protocolParams.FailureDetectionTime)
                                for idxConnDelay = 1:length(protocolParams.connection_time)
                                    for idxSignalingAfterRachDelay = 1:length(protocolParams.signalingAfterRachTime)
                                        discDelay               = protocolParams.discovery_time(idxDiscDelay);
                                        failureDetectionDelay   = protocolParams.FailureDetectionTime(idxFailureDetectionDelay);
                                        connDelay               = protocolParams.connection_time(idxConnDelay);
                                        signalingAfterRachDelay = protocolParams.signalingAfterRachTime(idxSignalingAfterRachDelay);
                                        frameHopCount           = protocolParams.frameHopCount;
                                        frameDeliveryDelay      = protocolParams.frameDeliveryDelay;

                
                                        for ue_idx = 1:params.numUE   %storing outage probability and duration for each user
                                            out_prob_analysis = outage_probability_analysis(ue_idx,idxDiscDelay, idxFailureDetectionDelay, idxConnDelay, idxSignalingAfterRachDelay);
                                            out_prob_analysis_wo_cf = outage_probability_analysis_wo_cf(ue_idx,idxDiscDelay, idxFailureDetectionDelay, idxConnDelay, idxSignalingAfterRachDelay);
                                            out_dur_analysis = outage_duration_analysis(ue_idx,idxDiscDelay, idxFailureDetectionDelay, idxConnDelay, idxSignalingAfterRachDelay);
                                            out_dur_analysis_wo_cf = outage_duration_analysis_wo_cf(ue_idx,idxDiscDelay, idxFailureDetectionDelay, idxConnDelay, idxSignalingAfterRachDelay);                        
                                            min_rate_req = params.r_min;
                                            p_fac = params.p_fac;
                                            formatSpec = '%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%.16f,%.16f,%.16f,%.16f,%.16f\n';
                                            fprintf(fileID,formatSpec,ue_idx, lambda_BS(idxBSDensity),lambda_UE_sub6(idxUEDensity),numBlockers,...
                                                deployRange,discDelay,failureDetectionDelay,connDelay,...
                                                signalingAfterRachDelay,frameHopCount,frameDeliveryDelay,...
                                                min_rate_req, p_fac, params.loss_pc_thresh, out_dur_analysis_wo_cf,out_prob_analysis_wo_cf,out_dur_analysis,out_prob_analysis);
                                        end
                                    end
                                end
                            end
                        end
                        fclose(fileID);
                    end
                end
            end
        end
    end
end
tEnd = toc(tStart);
fprintf('Total runtime: %f seconds\n',tEnd)