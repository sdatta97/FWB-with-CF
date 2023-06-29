close all;
clear;
params.simTime = 10*60; %sec Total Simulation time should be more than 100.
%% Room Setup, UE placement, UE height
% We are considering an outdoor scenario where the UE is located at the
% center and gNBs are distributed around the UE. We only need to consider
% the coverageRange amount of distance from the UE.
params.coverageRange = 100;
params.coverageRange_sub6 = 1000;
length_area = 2*params.coverageRange;   width_area = 2*params.coverageRange;
height_transmitter = 5;
params.areaDimensions = [width_area, length_area, height_transmitter];
%%UE locations 
params.numUE = 1;
params.RUE = 0;  %params.coverageRange * sqrt(rand(params.numUE,1)); %location of UEs (distance from origin)
params.angleUE = 2*pi*rand(params.numUE,1);%location of UEs (angle from x-axis)
params.UE_locations = [params.RUE.*cos(params.angleUE), params.RUE.*sin(params.angleUE)];
params.hr = 1.4; %height receiver (UE), approximately the height a human holds the phone
params.ht = height_transmitter; %height transmitter (BS)
height_transmitter_sub6 = 4;
params.ht_sub6 = height_transmitter_sub6; %height transmitter (BS)
rmin = 1e9;
params.r_min = rmin*ones(params.numUE,1);  %stores min rate requirement for all mmWave users

%% gNB locations
% params.numGNB = 10;
% lambda_BS =[200,300,400,500]; %densityBS
lambda_BS = 50:50:200; %densityBS
lambda_UE_sub6 = lambda_BS./4;
rate_anal_mmW_user = zeros(size(lambda_BS));
rate_sim_mmW_user = zeros(size(lambda_BS));    
for idxBSDensity = 1:length(lambda_BS)
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
    n = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange_sub6/1000)^2);
    while (n==0)
        n = poissrnd(lambda_BS(idxBSDensity)*pi*(params.coverageRange_sub6/1000)^2);       
    end
    params.numGNB_sub6 = n;
    params.RgNB_sub6 = params.coverageRange_sub6 * sqrt(rand(params.numGNB_sub6 - params.numGNB,1)); %location of gNBs (distance from origin)
    % params.RgNB = (2*params.coverageRange/3) * ones(params.numGNB,1); %location of gNBs (distance from origin)
    params.angleGNB_sub6 = 2*pi*rand(params.numGNB_sub6 - params.numGNB,1);%location of gNBs (angle from x-axis)
    params.locationsBS_sub6 = [params.RgNB_sub6.*cos(params.angleGNB_sub6), params.RgNB_sub6.*sin(params.angleGNB_sub6)];  
    params.num_antennas_per_gNB = 8;
    % numGNB = 10;
    % RgNB = params.coverageRange * sqrt(rand(numGNB-params.numGNB,1)); %location of gNBs (distance from origin)
    % angleGNB = 2*pi*rand(numGNB-params.numGNB,1);%location of gNBs (angle from x-axis)
    % locationsBS = [params.RgNB.*cos(params.angleGNB), params.RgNB.*sin(params.angleGNB)];
    % % params.numGNB = numGNB;
    % params.RgNB = [params.RgNB;RgNB]; %location of gNBs (distance from origin)
    % params.angleGNB = [params.angleGNB;angleGNB];%location of gNBs (angle from x-axis)
    % params.locationsBS = [params.RgNB.*cos(params.angleGNB), params.RgNB.*sin(params.angleGNB)];
    
    
    params.numUE_sub6 = poissrnd(lambda_UE_sub6(idxBSDensity)*pi*(params.coverageRange_sub6/1000)^2);
    params.RUE_sub6 = params.coverageRange_sub6*sqrt(rand(params.numUE_sub6,1)); %location of UEs (distance from origin)
    params.angleUE_sub6 = 2*pi*rand(params.numUE_sub6,1);%location of UEs (angle from x-axis)
    params.UE_locations_sub6 = [params.RUE_sub6.*cos(params.angleUE_sub6), params.RUE_sub6.*sin(params.angleUE_sub6)];        
    rmin_sub6 = 1e6;
    params.r_min_sub6 = rmin_sub6*ones(params.numUE_sub6,1);  %stores min rate requirement for all sub-6 users
    % params.numUE_sub6_max = 40;
    params.numUE_sub6_max = params.numUE_sub6;
    %% PHY layer params
    params.scs_mmw = 2e9;     %not using this parameter now
    params.scs_sub6 = 1e8;   %sub-6 GHz bandwidth 100 MHz
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
    
    
    %% Sequence of events
    % https://nyu0-my.sharepoint.com/:o:/g/personal/mfo254_nyu_edu/EjQnC6mtkpNKt7F_2dbwjrEBzPBXZmDauUCn6EjiUwhl4Q?e=rpnrVT
    
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
    tic
    dataBS_mobile = [];
    for i = 1:params.numUE
        dataBS_mobile = [dataBS_mobile; computeBlockageEvents(params,i)];
    end
    % [phy_channel_mmw, phy_channel_sub6] = computePhysicalChannels(params);
    
    [BETA, ricianFactor] = channel_cellfree_2v2(params.numUE+params.numUE_sub6,params.numGNB_sub6, params.num_sc_sub6*params.scs_sub6, [params.locationsBS;params.locationsBS_sub6], [params.UE_locations;params.UE_locations_sub6]);
    params.BETA = BETA;
    params.ricianFactor = ricianFactor;
    fprintf('Blocker generation, physical blockage and channel computation done : %f seconds\n',toc)
   
    %% protocol params from other paper
    frac = (mean(params.hb)-params.hr)/(params.ht-params.hr);
    protocolParams.theta = 2*params.V.*params.lambdaBlockers*frac/pi;
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
                    % plos2 = pLoS2(params.locationsBS, [params.UE_locations;params.UE_locations_sub6], theta,omega,psi);
                    % plos = zeros(params.numUE,1);
                    % [~,idx_max] = maxk(plos2(:,1:params.numUE),2,1);
                    % for k = 1:params.numUE
                    %     plos(k) = prod(plos2(idx_max(:,k),k),1);
                    %     % plos(k) = mean(plos2(idx_max(:,k),k),1);
                    % end
                    plos2 = zeros(params.numGNB_sub6,params.numUE+params.numUE_sub6);
                    plos = ones(params.numUE+params.numUE_sub6,1);
                    % rate_anal = rate_analytical(params, plos2, plos);
                    rate_anal = rate_analytical(params, plos2, plos,BETA,ricianFactor);
                    % rate_anal = rate_analyticalv2(params, plos2, plos,BETA,ricianFactor);
                    rate_anal_mmW_user(idxBSDensity) = rate_anal(1);
                    discovery_delay = protocolParams.discovery_time(idxDiscDelay);
                    failureDetectionDelay = protocolParams.FailureDetectionTime(idxFailureDetectionDelay);
                    connection_setup_delay = protocolParams.connection_time(idxConnDelay);
                    signalingAfterRachDelay = protocolParams.signalingAfterRachTime(idxSignalingAfterRachDelay);
                    [discoveredTimes, bsBlockageTimes] = computeDiscoveredTimes(dataBS_mobile,params,discovery_delay,failureDetectionDelay);
                    %% Simulation Initialization
                    disp('=====================================');
                    disp('Starting Simulation:')
                    tic
                    % numBS = size(dataBS_mobile,1);
                    numBS = params.numGNB;
                    % numBS_mmW = params.numGNB;
                    % numBS = params.numGNB_sub6;
                    numUE = params.numUE;
                    numUE_sub6 = params.numUE_sub6;
                    % fprintf('numBS mmW: %d. \n',numBS_mmW)
                    fprintf('numBS: %d. \n',numBS)
                    fprintf('DiscD: %.3f, FailDet: %.3f, ConnD: %.3f, SingAfterRach: %.3f\n',...
                    discovery_delay,failureDetectionDelay,connection_setup_delay,signalingAfterRachDelay);
                    
                    %initial setup
                    currentTime=0;
                    % Protocol related things
                    % bsPriorities = 1:numBS; %Priority of bs in terms of starting a connection decision
                    bsPriorities = repmat(1:1:numBS,[numUE,1]); %Priority of bs in terms of starting a connection decision
                    % bsPriorities = repmat(1:1:numBS_mmW,[numUE,1]); %Priority of bs in terms of starting a connection decision
                    % [~,bsPriorities] = sort(link_rates_mmw,'descend'); %Priority of bs in terms of starting a connection decision
                    % bsLastConnectionTimes = -100*ones(1,numBS); %When was the last time this bs was in connected state
                    bsLastConnectionTimes = -100*ones(numUE,numBS); %When was the last time this bs was in connected state
                    % bsLastConnectionTimes = -100*ones(numUE,numBS_mmW); %When was the last time this bs was in connected state
                    
                    % Discovery computations
                    % link = cell(numBS,1);
                    link = cell(numUE,numBS);
                    for ue_idx = 1:numUE
                        for idxBS = 1:numBS
                            % link{idxBS}.discoveredTimes = discoveredTimes{idxBS};
                            % link{idxBS}.discovery_state = discoveryStatus(discoveredTimes,idxBS,currentTime);
                            % link{idxBS}.nonBlockedTimes =  bsBlockageTimes{idxBS};
                            % link{idxBS}.blockageStatus = blockageStatus(bsBlockageTimes,idxBS,currentTime);
                            link{ue_idx,idxBS}.discoveredTimes = discoveredTimes{(ue_idx-1)*numBS + idxBS};
                            link{ue_idx,idxBS}.discovery_state = discoveryStatus(discoveredTimes,idxBS,currentTime, ue_idx, numBS);
                            link{ue_idx,idxBS}.nonBlockedTimes = bsBlockageTimes{(ue_idx-1)*numBS + idxBS};
                            link{ue_idx,idxBS}.blockageStatus = blockageStatus(bsBlockageTimes,idxBS,currentTime, ue_idx, numBS);    % if link is discovered the next event time will be next blockage
                            % arrival. If it is not discovered, next event time will be discovery
                            % completion time.
                            if link{ue_idx,idxBS}.discovery_state %link is discovered
                                next_blockage = link{ue_idx,idxBS}.discoveredTimes(3,find(link{ue_idx,idxBS}.discoveredTimes(3,:)> currentTime,1));
                                link{ue_idx,idxBS}.nextEventTime = next_blockage;
                            else
                                next_discovery = link{ue_idx,idxBS}.discoveredTimes(1,find(link{ue_idx,idxBS}.discoveredTimes(1,:)> currentTime,1));
                                link{ue_idx,idxBS}.nextEventTime = next_discovery;
                            end
                        end
                    end
                    % % UE Primary/Secondary STATES:
                    % % Not-Connected(idle) = 0
                    % % Connected = 1
                    % % Establishing Connection = 2
                    % % Blockage Detection State = 3
                    % % Recovery Attempt state (Beam Failure Declared) = 4
                    % % MCG Failure converting secondary to primary = 5
                    
                    
                    
                    % UE states
                    UE.numGNB = params.numGNB;
                    %General timers for primary and secondary
                    %UE.RACH_eff = RACH_eff;
                    UE.failureDetectionDelay = failureDetectionDelay;
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
                    
                    % Try RACH if BS available
                    % if UE.primaryConnectionState == 0
                    %     UE = tryPrimaryConnecting(UE,currentTime,link,bsPriorities,bsLastConnectionTimes);
                    % end
                    for ue_idx = 1:numUE
                        if UE.primaryConnectionState(ue_idx) == 0
                            UE = tryPrimaryConnecting(UE,currentTime,link,bsPriorities,bsLastConnectionTimes,ue_idx);
                        end
                    end
                    % Try secondary RACH if BS available
                    % if UE.primaryConnectionState == 0
                    %     UE = trySecondaryConnecting(UE,currentTime,link,bsPriorities,bsLastConnectionTimes);
                    % end
                    for ue_idx = 1:numUE
                        if UE.secondaryConnectionState(ue_idx) == 0
                            UE = trySecondaryConnecting(UE,currentTime,link,bsPriorities,bsLastConnectionTimes,ue_idx);
                        end
                    end
                    currentTime=0;
                    %% Simulation next step
                    nextEventTimes = getNextEvents(link,params.simTime);
                    nextEventTime = min(nextEventTimes(:));
                    %need to check Ue event times as well
                    UEEventTimes = [UE.primaryNextEventTime;UE.secondaryNextEventTime];
                    if min(UEEventTimes) > currentTime
                        nextEventTime = min(UEEventTimes);
                    end
                    prevTime = currentTime;
                    currentTime = nextEventTime;
                    if currentTime <= prevTime
                        disp('Error infinite loop use ctrl c')
                    end
                    %Physical + Discovery Updates of all links,
                    link = updatePhysicalDiscovery(currentTime,link,discoveredTimes,bsBlockageTimes);
    
                    % [BETA, phy_channel_sub6, phy_channel_sub6_est] = computePhysicalChannels_sub6v2(params,link);
                    % [phy_channel_sub6, phy_channel_sub6_est] = computePhysicalChannels_sub6v2(params,link,BETA,ricianFactor);
                    [phy_channel_sub6, phy_channel_sub6_est] = computePhysicalChannels_sub6v2(params,link);
                    rate_sim = compute_link_ratesv3(BETA,phy_channel_sub6, phy_channel_sub6_est, params.scs_sub6,1,ones(params.numUE,1));
                    % rate_sim = compute_link_ratesv4(BETA,phy_channel_sub6, phy_channel_sub6_est, params.scs_sub6,1,ones(params.numUE,1));
                    % rate_sim = compute_link_ratesv4(params,phy_channel_sub6, phy_channel_sub6_est,1,ones(params.numUE,1));
                    rate_sim_mmW_user(idxBSDensity) = rate_sim(1);
                end
            end
        end
    end
end