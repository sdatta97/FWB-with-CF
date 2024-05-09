close all;
clear;
tStart = tic;
% aID = getenv('SLURM_ARRAY_TASK_ID');
% % This is for running on a cluster in parallel
% % the bash script should give the aID as input
% if (isempty(aID))
%     warning('aID is empty. Trying SLURM ID.')
%     aID = getenv('SLURM_ARRAY_TASK_ID');
% end
% if(isempty(aID))
%     warning('aID is empty. Replacing it with 0010.')
%     aID = '0022';
% end
% RNG seed.
% rng(str2double(aID),'twister');
% aID = randi([0, 99]);
for aID = 1:99
    rng(aID,'twister');
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
    params.Band = 100e6; %Communication bandwidth
    
    
    %% Define simulation setup
    
    %Angular standard deviation in the local scattering model (in radians)
    params.ASD_varphi = deg2rad(30); %azimuth angle
    params.ASD_theta = 0; %deg2rad(15);  %elevation angle
    
    %Total uplink transmit power per UE (mW)
    params.p = 100;
    
    params.num_antennas_per_gNB = 64;
    
    %Total downlink transmit power per AP (mW)
    params.rho_tot = 10^(3.6)*params.num_antennas_per_gNB; %200
    % rho_tot_arr = [10:10:100, 200:100:1000, 2000:1000:10000];
    
    %Power factor division
    % p_fac_arr = 10; %10.^(0:1:4);
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
    
    
    params.coverageRange_sub6 = 430;
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
    % lambda_BS = 50:50:200;%densityBS
    lambda_BS = 25;
    % num_BS_arr = [2,5,10,20]; %densityBS
    % numUE_sub6_arr = 2:2:10;
    % numUE_sub6_arr = 10;
    lambda_UE_sub6 = [250:250:1000, 1500, 2000]; % [30:20:90, 100]; %:100:2000;
    % for idxnumUEsub6 = 1:length(numUE_sub6_arr)
    params.loss_pc_thresh = 10;
    params.Lmax=4;
    lb_thresh = 0:0.1:1;
    for idxBSDensity = 1:length(lambda_BS)
        %% gNB locations
        % params.numGNB = 10;
        n = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange/1000)^2);
        while (n==0)
            n = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange/1000)^2);       
        end
        % params.numGNB = 1;
        params.numGNB = n;
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
            params.scs_sub6 = [0.5*(params.Band), 0.5*(params.Band)];   %sub-6 GHz bandwidth 100 MHz
            params.num_sc_mmw = 1;    %not using this parameter now
            params.num_sc_sub6 = 2;   %sub-6 GHz considered as one full band
         
           
            N = params.num_antennas_per_gNB;  % antennas per AP
            L = params.numGNB_sub6;
            K_mmW = params.numUE;
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
            
            % [gainOverNoisedB,R_gNB,R_ue_mmW,R_ue_sub6,pilotIndex,D,D_small,APpositions,UEpositions,distances] = generateSetup(params,1,str2double(aID));
            [gainOverNoisedB,R_gNB,R_ue_mmW,R_ue_sub6,pilotIndex,D,D_small,APpositions,UEpositions,distances] = generateSetup(params,1,aID);                           
            params.BETA = db2pow(gainOverNoisedB);   
            params.D = D;
            params.R_gNB = R_gNB;
            params.R_ue_mmW = R_ue_mmW;
            params.R_ue_sub6 = R_ue_sub6;
            [channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW] = computePhysicalChannels_sub6_MIMO(params); 
            % rate_dl_before_handoff = compute_link_rates_MIMO(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,1,0);                                              
            rate_dl_before_handoff = compute_link_rates_MIMO_mmse(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,1,0);                                              
            for idxlbthres = 1:length(lb_thresh)
                lb_thres = lb_thresh(idxlbthres);
                %% Recording the Results
        
                %Taking care of folder directory creation etc
                dataFolder = 'resultData';
                impactFolder = strcat(dataFolder,'/StatsResults');
                if not(isfolder(impactFolder))
                    mkdir(impactFolder)
                end
                result_string = ['/rate_stats_',num2str(K_mmW),...
                    'UE_',num2str(lambda_BS(idxBSDensity)),...
                    'lambdaBS_',num2str(lambda_UE_sub6(idxUEDensity)),...
                    'lambdaUE_',num2str(100*lb_thres),...
                    'lb_thresh_',num2str(aID),'mmse'];    
                %Since we are mostly interested in blockage probability, we want to
                %transfer the data quickly to our local machine from server. We will save
                %the results as a txt file.
                recording_text_file_string = strcat(impactFolder,result_string,'.csv');
                fileID = fopen(recording_text_file_string,'w');
                output_categories = ['lambdaBS,','lambdaUE,','loss_thresh,',...
                    'mean_rate,','threshold rate\n'];  
                fprintf(fileID,output_categories);
                mean_rate_before_handoff = mean(rate_dl_before_handoff((1+K_mmW):end));
                mean_rate_before_handoff_threshold = quantile(rate_dl_before_handoff((1+K_mmW):end),lb_thres);
                formatSpec = '%d,%d,%.16f,%.16f,%.16f\n';
                fprintf(fileID,formatSpec,lambda_BS(idxBSDensity),lambda_UE_sub6(idxUEDensity),...
                   (lb_thres*100),mean_rate_before_handoff,mean_rate_before_handoff_threshold);                      
                fclose(fileID);
            end
        end
    end
    tEnd = toc(tStart);
    fprintf('Total runtime: %f seconds\n',tEnd)
end