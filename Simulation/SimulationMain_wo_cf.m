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
height_transmitter_sub6 = 4;
params.areaDimensions_sub6 = [width_area_sub6, length_area_sub6, height_transmitter_sub6];
%%UE location
params.numUE = 1;
params.RUE = 0;  %params.coverageRange * sqrt(rand(params.numUE,1)); %location of UEs (distance from origin)
params.angleUE = 2*pi*rand(params.numUE,1);%location of UEs (angle from x-axis)
params.UE_locations = [params.RUE.*cos(params.angleUE), params.RUE.*sin(params.angleUE)];

params.hr = 1.4; %height receiver (UE), approximately the height a human holds the phone
params.ht = height_transmitter; %height transmitter (BS)
lambda_BS = 200:200:5000; %densityBS
for idxBSDensity = 1:length(lambda_BS)
    %% gNB locations
    % params.numGNB = 10;
    n = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange/1000)^2);
    while (n==0)
        n = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange/1000)^2);       
    end
    params.numGNB = n;
    % params.numGNB = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange/1000)^2);
    % params.numGNB = floor(lambda_BS(idxBSDensity)*pi*(params.coverageRange/1000)^2);
    params.RgNB = params.coverageRange * sqrt(rand(params.numGNB,1)); %location of gNBs (distance from origin)
    % params.RgNB = (2*params.coverageRange/3) * ones(params.numGNB,1); %location of gNBs (distance from origin)
    params.angleGNB = 2*pi*rand(params.numGNB,1);%location of gNBs (angle from x-axis)
    params.locationsBS = [params.RgNB.*cos(params.angleGNB), params.RgNB.*sin(params.angleGNB)];
    
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
    % a_list = theta_list.*2*R/3;                       %Blocker Arrivals
    % self_blockage = 5/6;
    % protocolParams.psi_list = 1000./(protocolParams.discovery_time + params.mu);
    
    outage_probability_analysis_wo_cf = zeros(params.numUE,length(protocolParams.discovery_time),length(protocolParams.connection_time));
    
    %% Simulation
    for idxDiscDelay = 1:length(protocolParams.discovery_time)
        for idxFailureDetectionDelay = 1:length(protocolParams.FailureDetectionTime)
            for idxConnDelay = 1:length(protocolParams.connection_time)
                for idxSignalingAfterRachDelay = 1:length(protocolParams.signalingAfterRachTime)
                    theta = protocolParams.theta;
                    omega = 1/(protocolParams.connection_time(idxConnDelay) + protocolParams.signalingAfterRachTime(idxSignalingAfterRachDelay) + protocolParams.FailureDetectionTime(idxFailureDetectionDelay));
                    psi = 1/(protocolParams.discovery_time(idxDiscDelay) + 1/params.mu);
                    for k = 1:params.numUE
                       % plos3 = pLoS3(params.locationsBS, params.UE_locations(k,:), theta,omega,psi,1:params.numGNB);
                       plos4 = pLoS4(params.locationsBS, params.UE_locations(k,:), theta,omega,psi,1:params.numGNB);
                       % outage_probability_analysis_wo_cf(k,idxDiscDelay, idxFailureDetectionDelay, idxConnDelay, idxSignalingAfterRachDelay) = outage_probability_analysis_wo_cf(k,idxDiscDelay, idxConnDelay) + plos3;
                       outage_probability_analysis_wo_cf(k,idxDiscDelay, idxFailureDetectionDelay, idxConnDelay, idxSignalingAfterRachDelay) = outage_probability_analysis_wo_cf(k,idxDiscDelay, idxConnDelay) + plos4;
                    end
                end
            end
        end
    end
    
    %% Mobile blockage events
    tic
    dataBS_mobile = [];
    for i = 1:(params.numUE)
        dataBS_mobile = [dataBS_mobile; computeBlockageEvents(params,i)];
    end
    % [phy_channel_mmw, phy_channel_sub6] = computePhysicalChannels(params);

    fprintf('Blocker generation, physical blockage and channel computation done : %f seconds\n',toc)

    %% Create Discrete Time Event Simulation input
    simInputs.params = params;
    simInputs.dataBS_mobile = dataBS_mobile;
    simInputs.protocolParams = protocolParams;
    simOutputs = cell(length(protocolParams.discovery_time),...
                      length(protocolParams.FailureDetectionTime),...
                      length(protocolParams.connection_time),...
                      length(protocolParams.signalingAfterRachTime));              
    simOutputs_wo_cf = cell(length(protocolParams.discovery_time),...
                      length(protocolParams.FailureDetectionTime),...
                      length(protocolParams.connection_time),...
                      length(protocolParams.signalingAfterRachTime));   
    % for ue_idx = 1:params.numUE
    for idxDiscDelay = 1:length(protocolParams.discovery_time)
        for idxFailureDetectionDelay = 1:length(protocolParams.FailureDetectionTime)
            for idxConnDelay = 1:length(protocolParams.connection_time)
                for idxSignalingAfterRachDelay = 1:length(protocolParams.signalingAfterRachTime)
                    simInputs.discovery_delay = protocolParams.discovery_time(idxDiscDelay);
                    simInputs.failureDetectionDelay = protocolParams.FailureDetectionTime(idxFailureDetectionDelay);
                    simInputs.connection_setup_delay = protocolParams.connection_time(idxConnDelay);
                    simInputs.signalingAfterRachDelay = protocolParams.signalingAfterRachTime(idxSignalingAfterRachDelay);
                    simOutputs_wo_cf{idxDiscDelay,idxFailureDetectionDelay,...
                        idxConnDelay,idxSignalingAfterRachDelay} =  discreteSimulatorv2(simInputs); 
                end
            end
        end
    end
    %% Recording the Results

    %Taking care of folder directory creation etc
    dataFolder = 'resultData';
    outageFolder = strcat(dataFolder,'/outageResults_new');
    eventFolder = strcat(dataFolder,'/allResults_new');
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
    numBS = size(params.locationsBS,1);
    numBlockers = params.numBlockers;
    result_string = ['/results_',num2str(numUE),...
        'UE_',num2str(lambda_BS(idxBSDensity)),...
        'lambdaBS_',num2str(numBlockers), 'Blockers_randomHeight_', num2str(aID)];
    results_save_string = strcat(eventFolder,result_string,'.mat');
    save(results_save_string,'simOutputs','protocolParams','dataDescription')

    %Since we are mostly interested in blockage probability, we want to
    %transfer the data quickly to our local machine from server. We will save
    %the results as a txt file.
    recording_text_file_string = strcat(outageFolder,result_string,'.csv');
    fileID = fopen(recording_text_file_string,'w');
    output_categories = ['UE idx,','lambdaBS,','numBlockers,',...
        'discoveryDelay,','failureDetectionDelay,','connectionSetupDelay,',...
        'signalingAfterRachDelay,','frameHopCount,','frameDeliveryDelay,'...
        ,'outageProbability_wo_CF_Analysis,','meanOutageDuration_wo_CF,','outageProbability_wo_CF\n'];

    fprintf(fileID,output_categories);

    for idxDiscDelay = 1:length(protocolParams.discovery_time)
        for idxFailureDetectionDelay = 1:length(protocolParams.FailureDetectionTime)
            for idxConnDelay = 1:length(protocolParams.connection_time)
                for idxSignalingAfterRachDelay = 1:length(protocolParams.signalingAfterRachTime)
                    thisOutputs_wo_cf = simOutputs_wo_cf{idxDiscDelay,idxFailureDetectionDelay,idxConnDelay,idxSignalingAfterRachDelay};                    
                    numBS                   = size(thisOutputs_wo_cf.params.locationsBS,1); 
                    numBlockers             = thisOutputs_wo_cf.params.numBlockers;
                    discDelay               = thisOutputs_wo_cf.discovery_delay;
                    failureDetectionDelay   = thisOutputs_wo_cf.failureDetectionDelay;
                    connDelay               = thisOutputs_wo_cf.connection_setup_delay;
                    signalingAfterRachDelay = thisOutputs_wo_cf.signalingAfterRachDelay;
                    frameHopCount           = thisOutputs_wo_cf.frameHopCount;
                    frameDeliveryDelay      = thisOutputs_wo_cf.frameDeliveryDelay;
                    for ue_idx = 1:params.numUE   %storing outage probability and duration for each user
                        mean_outage_duration_wo_cf    = thisOutputs_wo_cf.mean_outage_duration_wo_cf(ue_idx);
                        outage_probability_wo_cf      = thisOutputs_wo_cf.outage_probability_wo_cf(ue_idx);
                        out_prob_analysis_wo_cf = outage_probability_analysis_wo_cf(ue_idx,idxDiscDelay, idxFailureDetectionDelay, idxConnDelay, idxSignalingAfterRachDelay);
                        formatSpec = '%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%.16f,%.16f\n';
                        fprintf(fileID,formatSpec,ue_idx, lambda_BS(idxBSDensity),numBlockers,...
                            discDelay,failureDetectionDelay,connDelay,...
                            signalingAfterRachDelay,frameHopCount,frameDeliveryDelay,...
                            out_prob_analysis_wo_cf, mean_outage_duration_wo_cf, outage_probability_wo_cf);
                    end
                end
            end
        end
    end
    fclose(fileID);
end
tEnd = toc(tStart);
fprintf('Total runtime: %f seconds\n',tEnd)