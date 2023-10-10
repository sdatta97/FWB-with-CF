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


%%
% load('params.mat')
params.simTime = 10*60; %sec Total Simulation time should be more than 100.
%% Room Setup, UE placement, UE height
% We are considering an outdoor scenario where the UE is located at the
% center and gNBs are distributed around the UE. We only need to consider
% the coverageRange amount of distance from the UE.
params.coverageRange = 100;
length_area = 2*params.coverageRange;   width_area = 2*params.coverageRange;
height_transmitter = 5;
params.areaDimensions = [width_area, length_area, height_transmitter];
%%UE location
params.numUE = 1;
params.RUE = 0;  %params.coverageRange * sqrt(rand(params.numUE,1)); %location of UEs (distance from origin)
params.angleUE = 2*pi*rand(params.numUE,1);%location of UEs (angle from x-axis)
params.UE_locations = [params.RUE.*cos(params.angleUE), params.RUE.*sin(params.angleUE)];

params.hr = 1.4; %height receiver (UE), approximately the height a human holds the phone
params.ht = height_transmitter; %height transmitter (BS)
height_transmitter_sub6 = 4;
params.ht_sub6 = height_transmitter_sub6; %height transmitter (BS)
rmin = 1e8;
params.r_min = rmin*ones(params.numUE,1);  %stores min rate requirement for all mmWave users
% params.r_min = rmin*rand(params.numUE,1);

%% gNB locations
lambda_BS = [200,300,400,500]; %densityBS
% lambda_BS =[200,300]; %densityBS
% lambda_BS = 200;
% num_BS_arr = [2,5,10,20]; %densityBS
numUE_sub6_arr = 2:2:10;
% numUE_sub6_arr = 10;
for idxnumUEsub6 = 1:length(numUE_sub6_arr)
    for idxBSDensity = 1:length(lambda_BS)
        % params.numGNB = 10;
        params.numGNB = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange/1000)^2);
        % params.numGNB = floor(lambda_BS(idxBSDensity)*pi*(params.coverageRange/1000)^2);
        params.RgNB = params.coverageRange * sqrt(rand(params.numGNB,1)); %location of gNBs (distance from origin)
        % params.RgNB = (2*params.coverageRange/3) * ones(params.numGNB,1); %location of gNBs (distance from origin)
        params.angleGNB = 2*pi*rand(params.numGNB,1);%location of gNBs (angle from x-axis)
        params.locationsBS = [params.RgNB.*cos(params.angleGNB), params.RgNB.*sin(params.angleGNB)];
        %%UE locations
        % params.numUE_sub6 = 10;
        params.numUE_sub6 = numUE_sub6_arr(idxnumUEsub6);
        params.RUE_sub6 = params.coverageRange * sqrt(rand(params.numUE_sub6,1)); %location of UEs (distance from origin)
        params.angleUE_sub6 = 2*pi*rand(params.numUE_sub6,1);%location of UEs (angle from x-axis)
        params.UE_locations_sub6 = [params.RUE_sub6.*cos(params.angleUE_sub6), params.RUE_sub6.*sin(params.angleUE_sub6)];
        rmin_sub6 = 1e6;
        params.r_min_sub6 = rmin_sub6*ones(params.numUE_sub6,1);  %stores min rate requirement for all sub-6 users

        %% PHY layer params
        params.scs_mmw = 2e9;     %not using this parameter now
        params.scs_sub6 = 20e6;   %sub-6 GHz bandwidth
        params.num_sc_mmw = 1;    %not using this parameter now
        params.num_sc_sub6 = 1;   %sub-6 GHz considered as one full band
        
        %% UE angular coverage range (full 360 coverage for now)
        lookAngleCell{1} = [0,360];
        
        %% Blocker Properties and Simulation Duration
        params.lambdaBlockers = 0.1; %How many blockers around
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
        protocolParams.discovery_time = [20 50]*10^(-3);
        % protocolParams.discovery_time = 50*10^(-3);
        
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
        protocolParams.connection_time = [10 20 50]*10^(-3);
        % protocolParams.connection_time = 50*10^(-3);
        
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
        % a_list = theta_list.*2*R/3;                       %Blocker Arrivals
        % self_blockage = 5/6;
        % protocolParams.psi_list = 1000./(protocolParams.discovery_time + params.mu);
        
        [BETA, ricianFactor] = channel_cellfree_2v2(params.numUE+params.numUE_sub6,params.numGNB, params.num_sc_sub6*params.scs_sub6, params.locationsBS, [params.UE_locations;params.UE_locations_sub6]);
        params.BETA = BETA;
        params.ricianFactor = ricianFactor;
        outage_probability_analysis = zeros(params.numUE,length(protocolParams.discovery_time),length(protocolParams.connection_time));
        outage_probability_analysis_wo_cf = zeros(params.numUE,length(protocolParams.discovery_time),length(protocolParams.connection_time));
        
        for idxDiscDelay = 1:length(protocolParams.discovery_time)
            for idxFailureDetectionDelay = 1:length(protocolParams.FailureDetectionTime)
                for idxConnDelay = 1:length(protocolParams.connection_time)
                    for idxSignalingAfterRachDelay = 1:length(protocolParams.signalingAfterRachTime)
                        theta = protocolParams.theta;
                        % omega = protocolParams.omega_list(idxConnDelay);
                        omega = 1/(protocolParams.connection_time(idxConnDelay) + protocolParams.signalingAfterRachTime(idxSignalingAfterRachDelay) + protocolParams.FailureDetectionTime(idxFailureDetectionDelay));
                        % omega = 1/protocolParams.connection_time(idxConnDelay); 
                        % psi = protocolParams.psi_list(idxDiscDelay);
                        psi = 1/(protocolParams.discovery_time(idxDiscDelay) + 1/params.mu);
                        plos2 = pLoS2(params.locationsBS, [params.UE_locations;params.UE_locations_sub6], theta,omega,psi);
                        plos = zeros(params.numUE,1);
                        % [~,idx_max] = mink(plos2(:,1:params.numUE),2,1);
                        % [~,idx_max] = mink(plos2(:,1:params.numUE),1,1);
                        [~,idx_max] = mink(plos2(:,1:params.numUE),params.numGNB,1);
                        for k = 1:params.numUE
                            % plos(k) = prod(plos2(:,k),1);
                            plos(k) = prod(plos2(idx_max(:,k),k),1);
                            % plos(k) = mean(plos2(idx_max(:,k),k),1);
                        end
                        % rate_dl = rate_analytical(params, plos2, plos);
                        rate_dl = rate_analytical(params, plos2, plos, BETA, ricianFactor);
                        for k = 1:params.numUE
                           plos3 = pLoS3(params.locationsBS, params.UE_locations(k,:), theta,omega,psi,idx_max);
                           % plos3 = pLoS3(theta,omega,params.coverageRange);
                           plos4 = pLoS4(params.locationsBS, params.UE_locations(k,:), theta,omega,psi,idx_max);
                           if (rate_dl(k) >= params.r_min(k) && all(rate_dl(1+params.numUE:params.numUE+params.numUE_sub6)' >= params.r_min_sub6(:)))
                               p1 = 1;
                           else
                               % p1 = 1-plos4;
                               % p1 = 1-plos(k);
                               p1 = 1-plos3;
                           end
                           outage_probability_analysis(k,idxDiscDelay, idxFailureDetectionDelay, idxConnDelay, idxSignalingAfterRachDelay) = outage_probability_analysis(k,idxDiscDelay, idxConnDelay) + (1-p1);
                           % outage_probability_analysis_wo_cf(k,idxDiscDelay, idxFailureDetectionDelay, idxConnDelay, idxSignalingAfterRachDelay) = outage_probability_analysis_wo_cf(k,idxDiscDelay, idxConnDelay) + plos4;
                           outage_probability_analysis_wo_cf(k,idxDiscDelay, idxFailureDetectionDelay, idxConnDelay, idxSignalingAfterRachDelay) = outage_probability_analysis_wo_cf(k,idxDiscDelay, idxConnDelay) + plos3;
                           % outage_probability_analysis_wo_cf(k,idxDiscDelay, idxFailureDetectionDelay, idxConnDelay, idxSignalingAfterRachDelay) = outage_probability_analysis_wo_cf(k,idxDiscDelay, idxConnDelay) + plos(k);
                        end
                    end
                end
            end
        end
        %% Recording the Results

        %Taking care of folder directory creation etc
        dataFolder = 'resultDataAnalysis';
        outageFolder = strcat(dataFolder,'/outageResults');
        eventFolder = strcat(dataFolder,'/allResults');
        if not(isfolder(dataFolder))
            mkdir(dataFolder)
        end
        if not(isfolder(outageFolder))
            mkdir(outageFolder)
        end
        if not(isfolder(eventFolder))
            mkdir(eventFolder)
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
        result_string = ['/results_',num2str(numUE),...
            'UE_',num2str(numUE_sub6),...
            'UE_sub6_',num2str(lambda_BS(idxBSDensity)),...
            'lambdaBS_',num2str(numBlockers), 'Blockers_randomHeight_', num2str(aID),'Min_rate', num2str(rmin)];
        % results_save_string = strcat(eventFolder,result_string,'.mat');
        % save(results_save_string,'simOutputs','protocolParams','dataDescription')

        %Since we are mostly interested in blockage probability, we want to
        %transfer the data quickly to our local machine from server. We will save
        %the results as a txt file.
        recording_text_file_string = strcat(outageFolder,result_string,'.csv');
        fileID = fopen(recording_text_file_string,'w');
        output_categories = ['UE idx,','numUEsub6,','lambdaBS,','numBlockers,',...
            'discoveryDelay,','failureDetectionDelay,','connectionSetupDelay,',...
            'signalingAfterRachDelay,','frameHopCount,','frameDeliveryDelay,'...
            'minRatereq,','outageProbability_wo_CF_Analysis,','outageProbabilityAnalysis,','meanOutageDuration_wo_CF,','outageProbability_wo_CF\n'];

        fprintf(fileID,output_categories);

        for idxDiscDelay = 1:length(protocolParams.discovery_time)
            for idxFailureDetectionDelay = 1:length(protocolParams.FailureDetectionTime)
                for idxConnDelay = 1:length(protocolParams.connection_time)
                    for idxSignalingAfterRachDelay = 1:length(protocolParams.signalingAfterRachTime)
                        % thisOutputs = simOutputs{idxDiscDelay,idxFailureDetectionDelay,idxConnDelay,idxSignalingAfterRachDelay};                    
                        % numBS                   = size(thisOutputs.params.locationsBS,1);                
                        % numBlockers             = thisOutputs.params.numBlockers;
                        % discDelay               = thisOutputs.discovery_delay;
                        % failureDetectionDelay   = thisOutputs.failureDetectionDelay;
                        % connDelay               = thisOutputs.connection_setup_delay;
                        % signalingAfterRachDelay = thisOutputs.signalingAfterRachDelay;
                        % frameHopCount           = thisOutputs.frameHopCount;
                        % frameDeliveryDelay      = thisOutputs.frameDeliveryDelay;
                        numBS                   = params.numGNB;                
                        numBlockers             = params.numBlockers;
                        discDelay               = protocolParams.discovery_time(idxDiscDelay);
                        failureDetectionDelay   = protocolParams.FailureDetectionTime(idxFailureDetectionDelay);
                        connDelay               = protocolParams.connection_time(idxConnDelay);
                        signalingAfterRachDelay = protocolParams.signalingAfterRachTime(idxSignalingAfterRachDelay);
                        frameHopCount           = protocolParams.frameHopCount;
                        frameDeliveryDelay      = protocolParams.frameDeliveryDelay;                   
                        for ue_idx = 1:params.numUE   %storing outage probability and duration for each user
                            out_prob_analysis = outage_probability_analysis(ue_idx,idxDiscDelay, idxFailureDetectionDelay, idxConnDelay, idxSignalingAfterRachDelay);
                            out_prob_analysis_wo_cf = outage_probability_analysis_wo_cf(ue_idx,idxDiscDelay, idxFailureDetectionDelay, idxConnDelay, idxSignalingAfterRachDelay);
                            min_rate_req = params.r_min(ue_idx);
                            formatSpec = '%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%.16f,%.16f\n';
                            fprintf(fileID,formatSpec,ue_idx,numUE_sub6, lambda_BS(idxBSDensity),numBlockers,...
                                discDelay,failureDetectionDelay,connDelay,...
                                signalingAfterRachDelay,frameHopCount,frameDeliveryDelay,...
                                min_rate_req, out_prob_analysis_wo_cf, out_prob_analysis);
                        end
                    end
                end
            end
        end
        fclose(fileID);
    end
end
tEnd = toc(tStart);
fprintf('Total runtime: %f seconds\n',tEnd)