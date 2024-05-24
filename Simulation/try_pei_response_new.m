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
%RNG seed.
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
    
    %Total downlink transmit power per AP (mW)
    % rho_tot_arr = [10:10:100, 200:100:1000, 2000:1000:10000];
    
    %Power factor division
    p_fac_arr = 1; %10^2; %10:10:100; %10.^(0:1:5);
    % params.p_fac = 10;
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
    params.deployRange = 20:20:100;
    %% Define simulation setup
    
    %Angular standard deviation in the local scattering model (in radians)
    params.ASD_varphi = deg2rad(30); %azimuth angle
    params.ASD_theta = 0; %deg2rad(15);  %elevation angle
    params.coverageRange = 100;
    length_area = 2*params.coverageRange;   
    width_area = 2*params.coverageRange;
    % length_area = 2*(params.deployRange + params.coverageRange);   
    % width_area = 2*(params.deployRange + params.coverageRange);
    height_transmitter = 5;
    params.areaDimensions = [width_area, length_area, height_transmitter];
    
    params.coverageRange_sub6 = 430;
    length_area_sub6 = 2*params.coverageRange_sub6;   
    width_area_sub6 = 2*params.coverageRange_sub6;
    % length_area_sub6 = 2*(params.deployRange + params.coverageRange_sub6);   
    % width_area_sub6 = 2*(params.deployRange + params.coverageRange_sub6);
    height_transmitter_sub6 = 4; % 5;
    params.areaDimensions_sub6 = [width_area_sub6, length_area_sub6, height_transmitter_sub6];
    for idxdeployRange = 1:length(params.deployRange)
        %%UE location
        params.numUE = 2;
        deployRange = params.deployRange(idxdeployRange);
        params.RUE =  deployRange*rand(params.numUE,1);%params.coverageRange*sqrt(rand(params.numUE,1)); %location of UEs (distance from origin)
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
        lambda_UE_sub6 = 1000; % [30:20:90, 100]; %:100:2000;
        % for idxnumUEsub6 = 1:length(numUE_sub6_arr)
        params.loss_pc_thresh = 10;
        params.Lmax=4;
        lb_thresh = 0.1; %0:0.25:1; %[0, 0.05, 0.1, 1]; %[0.05, 0.1]; %[0.1, 0.25, 0.5];
        for idxBSDensity = 1:length(lambda_BS)
           %% gNB locations
            % params.numGNB = 10;
        %     n = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange/1000)^2);
            n = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange/1000)^2);
            while (n==0)
        %         n = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange/1000)^2);       
                n = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange/1000)^2);       
            end
            params.numGNB = n;
            % params.numGNB = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange/1000)^2);
            % params.numGNB = floor(lambda_BS(idxBSDensity)*pi*(params.coverageRange/1000)^2);
            % params.RgNB = params.coverageRange * sqrt(rand(params.numGNB,1)); %location of gNBs (distance from origin)
            % params.RgNB = (params.deployRange+params.coverageRange) * sqrt(rand(params.numGNB,1)); %location of gNBs (distance from origin)
            % params.RgNB = (2*params.coverageRange/3) * ones(params.numGNB,1); %location of gNBs (distance from origin)
            % params.angleGNB = 2*pi*rand(params.numGNB,1);%location of gNBs (angle from x-axis)
            % params.locationsBS = [params.RgNB.*cos(params.angleGNB), params.RgNB.*sin(params.angleGNB)];
            numBS = params.numGNB;
            numUE = params.numUE;
            RgNB = params.coverageRange*sqrt(rand(numBS,numUE));
            angleGNB =  2*pi*rand(numBS,numUE);
            locationsBS = zeros(numBS*numUE,2);
            for k = 1:numUE
                locationsBS(numBS*(k-1)+1:numBS*k,:) = params.UE_locations(k,:) + [RgNB(:,k).*cos(angleGNB(:,k)), RgNB(:,k).*sin(angleGNB(:,k))];
            end
            params.RgNB = RgNB;
            params.angleGNB = angleGNB;
            params.locationsBS = locationsBS;
            n = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange_sub6/1000)^2);
            while (n<=params.numGNB) %(n==0)
                n = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange_sub6/1000)^2);       
            end
            params.numGNB_sub6 = n;
        %     params.RgNB_sub6 = params.coverageRange_sub6 * sqrt(rand(params.numGNB_sub6 - params.numGNB,1)); %location of gNBs (distance from origin)
            % params.RgNB_sub6 = params.coverageRange_sub6 * sqrt(rand(params.numGNB_sub6 - params.numGNB,1)); %location of gNBs (distance from origin)
            params.RgNB_sub6 = (norm(params.UE_locations(:,1) - params.UE_locations(:,2))/2 + params.coverageRange_sub6) * sqrt(rand(params.numGNB_sub6 - params.numGNB*params.numUE,1)); %location of gNBs (distance from origin)
            % params.RgNB = (2*params.coverageRange/3) * ones(params.numGNB,1); %location of gNBs (distance from origin)
            params.angleGNB_sub6 = 2*pi*rand(params.numGNB_sub6 - params.numGNB*params.numUE,1);%location of gNBs (angle from x-axis)
            % params.locationsBS_sub6 = [params.RgNB_sub6.*cos(params.angleGNB_sub6), params.RgNB_sub6.*sin(params.angleGNB_sub6)];  
            params.locationsBS_sub6 = mean(params.UE_locations,1) + [params.RgNB_sub6.*cos(params.angleGNB_sub6), params.RgNB_sub6.*sin(params.angleGNB_sub6)];  
            %     params.Lmax=n;
            for idxUEDensity = 1:length(lambda_UE_sub6)
                params.ue_rearranged = [];    
                %%UE locations
                n = poissrnd(lambda_UE_sub6(idxUEDensity)*pi*(params.coverageRange_sub6/1000)^2);
                while (n==0)
                    n = poissrnd(lambda_UE_sub6(idxUEDensity)*pi*(params.coverageRange_sub6/1000)^2);       
                end
                params.numUE_sub6 = n;
                % params.RUE_sub6 = params.coverageRange_sub6*sqrt(rand(params.numUE_sub6,1)); %location of UEs (distance from origin)
                params.RUE_sub6 = (norm(params.UE_locations(:,1) - params.UE_locations(:,2))/2+params.coverageRange_sub6)*sqrt(rand(params.numUE_sub6,1)); %location of UEs (distance from origin)
                params.angleUE_sub6 = 2*pi*rand(params.numUE_sub6,1);%location of UEs (angle from x-axis)
                params.UE_locations_sub6 =  mean(params.UE_locations,1) + [params.RUE_sub6.*cos(params.angleUE_sub6), params.RUE_sub6.*sin(params.angleUE_sub6)];   
                  %                 [gainOverNoisedB,R_gNB,R_ue_mmW,R_ue_sub6,pilotIndex,D,D_small,APpositions,UEpositions,distances] = generateSetup(params.numGNB,params.numGNB_sub6,params.numUE,params.numUE+params.numUE_sub6,params.num_antennas_per_gNB,params.N_UE_mmW,params.N_UE_sub6,params.coverageRange,params.coverageRange_sub6,params.tau_p,1,0,params.ASD_varphi,params.ASD_theta);
                            % [gainOverNoisedB,R_gNB,R_ue_mmW,R_ue_sub6,pilotIndex,D,D_small,APpositions,UEpositions,distances] = generateSetup(params,1,str2double(aID));
                [gainOverNoisedB,R_gNB,R_ue_mmW,R_ue_sub6,pilotIndex,D,D_small,APpositions,UEpositions,distances] = generateSetup(params,1,aID);
                
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
                
                N = params.num_antennas_per_gNB;  % antennas per AP
                L = params.numGNB_sub6;
                K_mmW = params.numUE;
                K = params.numUE + params.numUE_sub6;  % --Ground UEs
                numUE_sub6 = params.numUE_sub6;
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
                num_sc_sub6 = params.num_sc_sub6;
                params.user_sc_alloc = ones(K,num_sc_sub6);                                
                params.BETA = db2pow(gainOverNoisedB);   
                params.D = D;
                params.R_gNB = R_gNB;
                params.R_ue_mmW = R_ue_mmW;
                params.R_ue_sub6 = R_ue_sub6;
                
                for idxrmin = 1:length(rmin_arr)
                    for idx_p = 1:length(p_fac_arr)
                        for idxlbthres = 1:length(lb_thresh)
                            lb_thres = lb_thresh(idxlbthres);
                            params.lb_thres = lb_thres;
                            rmin = rmin_arr(idxrmin);
                            params.r_min = rmin*ones(params.numUE,1);  %stores min rate requirement for all mmWave users
                            rmin_sub6 = 35e6;
                            params.r_min_sub6 = rmin_sub6*ones(params.numUE_sub6,1);  %stores min rate requirement for all sub-6 users
                            params.rate_reduce_threshold = 5e7;
                            params.p_fac = p_fac_arr(idx_p);
                            params.p_fac_rearrange = 1; % 0.1*p_fac_arr(idx_p);  
                            [channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW] = computePhysicalChannels_sub6_MIMO(params); 
                            % for p_idx = 1:length(p_fac_arr)
                            %     params.p_fac = p_fac_arr(p_idx);
                            %     params.p_fac_rearrange = 0.1*p_fac_arr(p_idx); 
                            numUE = params.numUE;
                            sub6ConnectionState = ones(numUE,1);
                            for ue_idx = 1:numUE 
                                [~, ue_idxs_affected] = AP_reassign(params,ue_idx);
                                params.ue_rearranged = union(ue_idxs_affected, params.ue_rearranged);
                            end
                            % D = params.D;
                            % BETA = params.BETA;
                            % params.D = D(:,[(1:numUE)'; ue_idxs_affected]);
                            % params.BETA = BETA(:,[(1:numUE)'; ue_idxs_affected]);
                            % rate_dl_before_handoff = compute_link_rates_MIMO_mmse(params,channel_dl(:,ue_idxs_affected-numUE,:,:), channel_est_dl(:,ue_idxs_affected-numUE,:,:),channel_dl_mmW, channel_est_dl_mmW,ue_idx,sub6ConnectionState);                                              
                            % params.D = D;
                            % params.BETA = BETA;
                            rate_dl_before_handoff = compute_link_rates_MIMO_mmse(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,ue_idx,sub6ConnectionState);                                              
                            % lb = quantile(rate_dl_before_handoff(union((1+numUE):end,nonzeros((1:numUE)'.*sub6ConnectionState)))./params.Band,params.lb_thres);
                            lb = quantile(rate_dl_before_handoff(union(ue_idxs_affected,nonzeros((1:numUE)'.*sub6ConnectionState)))./params.Band,params.lb_thres);
                            bw_alloc = Band - rmin_sub6/lb;
                            params.scs_sub6(1) = bw_alloc;
                            params.scs_sub6(2) = Band - bw_alloc;
                            % params.ue_rearranged = ue_idxs_affected;
                            ues_not_affected = setdiff((1+numUE):(numUE+numUE_sub6),params.ue_rearranged);
                            % user_sc_alloc = ones(numUE+numUE_sub6,params.num_sc_sub6);                               
                            user_sc_alloc = params.user_sc_alloc; %zeros(numUE+numUE_sub6,1);                               
                            user_sc_alloc(ue_idx,1) = 1;
                            user_sc_alloc(ue_idx,2) = 0;
                            user_sc_alloc(ues_not_affected,1) = 1;
                            user_sc_alloc(ues_not_affected,2) = 0;
                            user_sc_alloc(ue_idxs_affected,1) = 0;
                            user_sc_alloc(ue_idxs_affected,2) = 1;
                            params.user_sc_alloc = user_sc_alloc;
                            sub6ConnectionState(ue_idx) = 1;
        %                         rate_dl_after_handoff = compute_link_rates_MIMO(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,ue_idx,sub6ConnectionState);                                              
                            rate_dl_after_handoff = compute_link_rates_MIMOv4(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,ue_idx,sub6ConnectionState);                                                         
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
                            result_string = strcat('/handoff_impact_',num2str(numUE),...
                                'UE_',num2str(lambda_BS(idxBSDensity)),...
                                'lambdaBS_',num2str(lambda_UE_sub6(idxUEDensity)),...
                                'lambdaUE_',num2str(deployRange),...
                                'deployRange_',num2str(numBlockers), 'Blockers_randomHeight_', num2str(aID),'Min_rate', num2str(rmin), "Pow_fac", num2str(params.p_fac), "lb_thres", num2str(100*params.lb_thres));
                    
                            %Since we are mostly interested in blockage probability, we want to
                            %transfer the data quickly to our local machine from server. We will save
                            %the results as a txt file.
                            recording_text_file_string = strcat(impactFolder,result_string,'.csv');
                            fileID = fopen(recording_text_file_string,'w');
                            output_categories = ['UE idx,','lambdaBS,','lambdaUE,','numBlockers,',...
                                'deployRange,','minRatereq,','powerFac,','lower_bound_thresh,', 'mean_rate_before_handoff,','mean_rate_after_handoff\n'];
                    
                            fprintf(fileID,output_categories);
                        
                            p_fac = params.p_fac;
                            formatSpec = '%d,%d,%d,%d,%f,%f,%f,%f,%.16f,%.16f\n';
                            fprintf(fileID,formatSpec,ue_idx, lambda_BS(idxBSDensity),lambda_UE_sub6(idxUEDensity),numBlockers,...
                                deployRange,rmin, p_fac, lb_thres, mean(rate_dl_before_handoff),mean(rate_dl_after_handoff((1+K_mmW):end)));
                            fclose(fileID);
                        end
                    end
                end
            end
        end
    end
end