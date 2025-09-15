% function simOutputs = discreteSimulator(simInputs, ue_idx)
function simOutputs = discreteSimulator_MIMO(simInputs) 
%% Inputs
params = simInputs.params;
protocolParams = simInputs.protocolParams;
dataBS_mobile = simInputs.dataBS_mobile; %{(ue_idx-1)*params.numGNB+1:ue_idx*params.numGNB,1};
numBS_mobile = simInputs.numBS_mobile;
r_min = params.r_min;
r_min_sub6 = params.r_min_sub6;
rate_reduce_threshold = params.rate_reduce_threshold;
Band = params.Band;
D = params.D;
ap_idxs = find(D(:,1));
ue_idxs = 1;
for a = 1:length(ap_idxs)
    ap_idx = ap_idxs(a);
    ue_idxs = union(ue_idxs,find(D(ap_idx,:)));
end
% BETA = params.BETA;
% ricianFactor = params.ricianFactor;
params.no_of_rea = 1;     % no.of channel realizations

discovery_delay         = simInputs.discovery_delay;
failureDetectionDelay   = simInputs.failureDetectionDelay;
connection_setup_delay  = simInputs.connection_setup_delay;
signalingAfterRachDelay = simInputs.signalingAfterRachDelay;

beamFailureRecoveryTimer = discovery_delay;

% RACH_eff = frameHopCount * frameDeliverDelay + RACH + signalingAfterRach
RACH_eff = protocolParams.frameHopCount * protocolParams.frameDeliveryDelay + ...
           connection_setup_delay + signalingAfterRachDelay;

%  discovery_delay        =   simInputs.discovery_delay;
%  connection_setup_delay =   simInputs.connection_setup_delay;
%  TTT_delay              =   simInputs.TTT_delay;
%  pkt_delay              =   simInputs.Pkt_delay;
% 


%% Aply discovery
% [discoveredTimes, bsBlockageTimes] = computeDiscoveredTimes(dataBS_mobile,params,discovery_delay,failureDetectionDelay);
[discoveredTimes, bsBlockageTimes] = computeDiscoveredTimes(dataBS_mobile,numBS_mobile,params,discovery_delay,failureDetectionDelay);
%% Simulation Initialization
disp('=====================================');
disp('Starting Simulation:')
tic
% numBS = size(dataBS_mobile,1);
numBS = params.numGNB;
numUE = params.numUE;
numUE_sub6 = params.numUE_sub6;
fprintf('numBS: %d. \n',numBS)
fprintf('DiscD: %.3f, FailDet: %.3f, ConnD: %.3f, SingAfterRach: %.3f\n',...
discovery_delay,failureDetectionDelay,connection_setup_delay,signalingAfterRachDelay);

%initial setup
currentTime=0;
% Protocol related things
% bsPriorities = 1:numBS; %Priority of bs in terms of starting a connection decision
bsPriorities = repmat(1:1:numBS,[numUE,1]); %Priority of bs in terms of starting a connection decision
% bsPriorities = zeros(numUE,numBS);
% for ue_idx = 1:numUE
%     r = zeros(numBS,1);
%     for ap_idx = 1:numBS
%         r(ap_idx) = sqrt(sum((params.UE_locations(ue_idx,:)-params.locationsBS(ap_idx,:)).^2)); 
%     end
%     [~, bsPriorities(ue_idx,:)] = sort(r,'ascend');
% end
% bsPriorities = repmat(1:1:numBS_mmW,[numUE,1]); %Priority of bs in terms of starting a connection decision
% [~,bsPriorities] = sort(link_rates_mmw,'descend'); %Priority of bs in terms of starting a connection decision
bsLastConnectionTimes = -100*ones(numUE,numBS); %When was the last time this bs was in connected state

% Discovery computations
link = cell(numUE,numBS);
% link = cell((numUE+numUE_sub6),numBS);

for ue_idx = 1:numUE
% for ue_idx = 1:(numUE+numUE_sub6)
    for idxBS = 1:numBS_mobile(ue_idx) %1:numBS        
        % link{ue_idx,idxBS}.discoveredTimes = discoveredTimes{(ue_idx-1)*numBS + idxBS};
        link{ue_idx,idxBS}.discoveredTimes = discoveredTimes{sum(numBS_mobile(1:(ue_idx-1)))+idxBS};
        % link{ue_idx,idxBS}.discovery_state = discoveryStatus(discoveredTimes,idxBS,currentTime, ue_idx, numBS);
        link{ue_idx,idxBS}.discovery_state = discoveryStatus(discoveredTimes,idxBS,currentTime, ue_idx, numBS_mobile);
        % link{ue_idx,idxBS}.nonBlockedTimes = bsBlockageTimes{(ue_idx-1)*numBS + idxBS};
        link{ue_idx,idxBS}.nonBlockedTimes = bsBlockageTimes{sum(numBS_mobile(1:(ue_idx-1))) + idxBS};
        % link{ue_idx,idxBS}.blockageStatus = blockageStatus(bsBlockageTimes,idxBS,currentTime, ue_idx, numBS);    % if link is discovered the next event time will be next blockage
        link{ue_idx,idxBS}.blockageStatus = blockageStatus(bsBlockageTimes,idxBS,currentTime, ue_idx, numBS_mobile);    % if link is discovered the next event time will be next blockage
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
UE.RACH_eff = RACH_eff;
UE.failureDetectionDelay = failureDetectionDelay;
UE.beamFailureRecoveryTimer = beamFailureRecoveryTimer;

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
UE.primaryConnectionState = zeros(numUE,1);
UE.primaryConnectionStateHistory = [UE.primaryConnectionStateHistory, UE.primaryConnectionState];
UE.primaryBSIdx = zeros(numUE,1);
UE.primaryTargetIdx = zeros(numUE,1);
UE.primaryNextEventTime = -100*ones(numUE,1);

%Secondary State variables
UE.secondaryConnectionState = zeros(numUE,1);
UE.secondaryConnectionStateHistory = [UE.secondaryConnectionStateHistory, UE.secondaryConnectionState];
UE.secondaryBSIdx = zeros(numUE,1);
UE.secondaryTargetIdx = zeros(numUE,1);
UE.secondaryNextEventTime = -100*ones(numUE,1);
UE.tmpMCGBSIdx = zeros(numUE,1);

%Sub 6 State variables
UE.sub6ConnectionState = zeros(numUE,1);
UE.sub6ConnectionStateHistory = [UE.sub6ConnectionStateHistory, UE.sub6ConnectionState];
UE.sub6NextEventTime = -100*ones(numUE,1);

% Try RACH if BS available
for ue_idx = 1:numUE
    if (numBS_mobile(ue_idx) > 0)
        if UE.primaryConnectionState(ue_idx) == 0
            % UE = tryPrimaryConnecting(UE,currentTime,link,bsPriorities,bsLastConnectionTimes,ue_idx);
            UE = tryPrimaryConnecting(UE,currentTime,link,bsPriorities,bsLastConnectionTimes,ue_idx,numBS_mobile);
        end
    end
end
% Try secondary RACH if BS available
for ue_idx = 1:numUE
    if (numBS_mobile(ue_idx) > 0)
        if UE.secondaryConnectionState(ue_idx) == 0
            % UE = trySecondaryConnecting(UE,currentTime,link,bsPriorities,bsLastConnectionTimes,ue_idx);
            UE = trySecondaryConnecting(UE,currentTime,link,bsPriorities,bsLastConnectionTimes,ue_idx,numBS_mobile);
        end
    end
end
%% Simulation next step
% nextEventTimes = getNextEvents(link,params.simTime);
nextEventTimes = getNextEvents(link,params.simTime,numBS_mobile);
nextEventTimes = nextEventTimes(nextEventTimes>0);
% link_mat = cell2mat(link);
% nextEventTimes = getNextEvents(mat2cell(link_mat(1:numUE,1:params.numGNB),[1,1]),params.simTime);
nextEventTime = min(nextEventTimes,[],"all"); %min(nextEventTimes(1:numUE,1:params.numGNB));
%need to check Ue event times as well
UEEventTimes = [UE.primaryNextEventTime;UE.secondaryNextEventTime];
if min(UEEventTimes) > currentTime
    nextEventTime = min(UEEventTimes);
end
while nextEventTime < params.simTime
    prevTime = currentTime;
    currentTime = nextEventTime;
    if currentTime <= prevTime
        disp('Error infinite loop use ctrl c')
    end
    %Physical + Discovery Updates of all links,
    % link = updatePhysicalDiscovery(currentTime,link,discoveredTimes,bsBlockageTimes);
    link = updatePhysicalDiscovery(currentTime,link,discoveredTimes,bsBlockageTimes,numBS_mobile);

    %Protocol updates
    for ue_idx = 1:numUE
        if (numBS_mobile(ue_idx) > 0)
            if UE.secondaryConnectionState(ue_idx) == 1
                if UE.secondaryNextEventTime(ue_idx) == currentTime
                    UE.secondaryConnectionState(ue_idx) = 3;
                    UE.secondaryConnectionStateHistory = [UE.secondaryConnectionStateHistory, UE.secondaryConnectionState];
                    UE.secondaryNextEventTime(ue_idx) = currentTime + UE.failureDetectionDelay;
                    UE.secondaryConnectionEnds = [UE.secondaryConnectionEnds, currentTime];
                    UE.secondaryConnectionEndIndices = [UE.secondaryConnectionEndIndices, ue_idx];
                    UE.secondaryEventTimes = [UE.secondaryEventTimes; currentTime];
                    UE.secondaryEventIndices = [UE.secondaryEventIndices; ue_idx];
                    UE.secondaryEventDescriptions = [UE.secondaryEventDescriptions; {'UE '+string(ue_idx)+': secondaryGotBlocked-NotDetectedYet-MeasurementsGoingOn'}];
                end
            elseif UE.secondaryConnectionState(ue_idx) == 3
                % a blockage occured but we still didnt detected it yet. We are
                % in measurement state but at this point our measurements are
                % done, if current time is our nextevent time
                if UE.secondaryNextEventTime(ue_idx) == currentTime
                    idxSecondaryBS = UE.secondaryBSIdx(ue_idx);
                    isBSRecovered = link{ue_idx,idxSecondaryBS}.blockageStatus;
                    if isBSRecovered
                        UE.secondaryConnectionState(ue_idx) = 1;
                        UE.secondaryConnectionStateHistory = [UE.secondaryConnectionStateHistory, UE.secondaryConnectionState];
                        UE.secondaryConnectionStarts = [UE.secondaryConnectionStarts, currentTime];
                        UE.secondaryConnectionStartIndices = [UE.secondaryConnectionStartIndices, ue_idx];
                        UE.secondaryBSIdx(ue_idx) = idxSecondaryBS;
                        UE.secondaryBSHistory = [UE.secondaryBSHistory, UE.secondaryBSIdx];
                        next_blockage = link{ue_idx,idxSecondaryBS}.discoveredTimes(3,find(link{ue_idx,idxSecondaryBS}.discoveredTimes(3,:)> currentTime,1));
                        UE.secondaryNextEventTime(ue_idx) = next_blockage;
                        UE.secondaryEventTimes = [UE.secondaryEventTimes; currentTime];
                        UE.secondaryEventIndices = [UE.secondaryEventIndices; ue_idx];
                        UE.secondaryEventDescriptions =  [UE.secondaryEventDescriptions; {'UE '+string(ue_idx)+': MeasurementDone-Recovered-b4-BFDT-expired-No-Failure-ReconnectedWithoutDelays'}];
                    else
                        %BS didnt recover now we declare beam faiure and go to
                        %beam recovery procedure.
                        % Now that we measured the problem, and declared beam
                        % failure we go to recovery states with our own BS. 
                        % Now LA will wait until it determines if any of the
                        % beams from the currentBS is suitable for connection
                        % re-establishment. 
                        UE.secondaryConnectionState(ue_idx) = 4;
                        UE.secondaryConnectionStateHistory = [UE.secondaryConnectionStateHistory, UE.secondaryConnectionState];
                        UE.secondaryNextEventTime(ue_idx) = currentTime + UE.beamFailureRecoveryTimer;
                        % In order to give this BS a priority in the future
                            lastConnectedBS = UE.secondaryBSHistory(ue_idx,end);
                            bsLastConnectionTimes(ue_idx,lastConnectedBS) = currentTime;
                            [~,bsPriorities(ue_idx,:)] = sort(bsLastConnectionTimes(ue_idx,:),'descend');
                            % [~,bsPriorities(ue_idx,:)] = sort(bsLastConnectionTimes(ue_idx,1:numBS_mmW),'descend');
                       % now event recording
                        UE.secondaryEventTimes = [UE.secondaryEventTimes; currentTime];
                        UE.secondaryEventIndices = [UE.secondaryEventIndices; ue_idx];
                        UE.secondaryEventDescriptions =  [UE.secondaryEventDescriptions; {'UE '+string(ue_idx)+': MeasurementDone-Declared-Beam-Failure-going2-BeamRecovery'}];
                    end
                end
            elseif UE.secondaryConnectionState(ue_idx) == 4
                if UE.secondaryNextEventTime(ue_idx) == currentTime
                    %Check if BS recovered
                    idxSecondaryBS = UE.secondaryBSIdx(ue_idx);
                    isBSRecovered = link{ue_idx,idxSecondaryBS}.blockageStatus;
                    if isBSRecovered
                        if link{ue_idx,idxSecondaryBS}.blockageStatus %&& link{idxSecondaryBS}.discovery_state
                            UE.secondaryConnectionState(ue_idx) = 2;
                            UE.secondaryConnectionStateHistory = [UE.secondaryConnectionStateHistory, UE.secondaryConnectionState];
                            UE.secondaryTargetIdx(ue_idx) = idxSecondaryBS;
                            UE.secondaryBSIdx(ue_idx) = 0;
                            UE.secondaryNextEventTime(ue_idx) = currentTime + UE.RACH_eff;
                            link{ue_idx,idxSecondaryBS}.nextEventTime = currentTime + UE.RACH_eff;
                            % now event recording
                            UE.secondaryEventTimes = [UE.secondaryEventTimes; currentTime];
                            UE.secondaryEventIndices = [UE.secondaryEventIndices; ue_idx];
                            UE.secondaryEventDescriptions =  [UE.secondaryEventDescriptions; {'UE '+string(ue_idx)+': Recovered-b4-RLF-expired-(Re)Connecting-to-BS'}];
                        else
                            disp('Arrived this line something might be wrong. DEBUG')
                        end
                    else
                        %Declare RLF to go idle state
                        UE.secondaryConnectionState(ue_idx) = 0;
                        UE.secondaryConnectionStateHistory = [UE.secondaryConnectionStateHistory, UE.secondaryConnectionState];
                        UE.secondaryBSIdx(ue_idx) = 0;
                        % now event recording
                        UE.secondaryEventTimes = [UE.secondaryEventTimes; currentTime];
                        UE.secondaryEventIndices = [UE.secondaryEventIndices; ue_idx];
                        UE.secondaryEventDescriptions =  [UE.secondaryEventDescriptions; {'UE '+string(ue_idx)+': Declared-RLF-Going2-NonConnectedState-CanConnect-any-available-BS'}];
                    end
                end
            end
        end
    end
    for ue_idx = 1:numUE
        if (numBS_mobile > 0)
            if UE.primaryConnectionState(ue_idx) == 1
                if UE.primaryNextEventTime(ue_idx) == currentTime
                    UE.primaryConnectionState(ue_idx) = 3;
                    UE.primaryConnectionStateHistory = [UE.primaryConnectionStateHistory, UE.primaryConnectionState];
                    UE.primaryNextEventTime(ue_idx) = currentTime + UE.failureDetectionDelay;
                    UE.primaryConnectionEnds = [UE.primaryConnectionEnds, currentTime];
                    UE.primaryConnectionEndIndices = [UE.primaryConnectionEndIndices, ue_idx];
                    UE.primaryEventTimes = [UE.primaryEventTimes; currentTime];
                    UE.primaryEventIndices = [UE.primaryEventIndices; ue_idx];
                    UE.primaryEventDescriptions = [UE.primaryEventDescriptions; {'UE '+string(ue_idx)+': primaryGotBlocked-NotDetectedYet-MeasurementsGoingOn-DataplaneInterruption'}];
                end
            elseif UE.primaryConnectionState(ue_idx) == 3
                if UE.primaryNextEventTime(ue_idx) == currentTime
                    idxPrimaryBS = UE.primaryBSIdx(ue_idx);
                    isBSRecovered = link{ue_idx,idxPrimaryBS}.blockageStatus;
                    if isBSRecovered
                        UE.primaryConnectionState(ue_idx) = 1;
                        UE.primaryConnectionStateHistory = [UE.primaryConnectionStateHistory, UE.primaryConnectionState];
                        UE.primaryConnectionStarts = [UE.primaryConnectionStarts, currentTime];
                        UE.primaryConnectionStartIndices = [UE.primaryConnectionStartIndices, ue_idx];
                        if UE.sub6ConnectionState(ue_idx) == 1
                            UE.sub6ConnectionEnds = [UE.sub6ConnectionEnds, currentTime];
                            UE.sub6ConnectionEndIndices = [UE.sub6ConnectionEndIndices, ue_idx];
                            UE.sub6ConnectionState(ue_idx) = 0;
                            UE.sub6ConnectionStateHistory = [UE.sub6ConnectionStateHistory, UE.sub6ConnectionState];
                        end
                        UE.primaryBSIdx(ue_idx) = idxPrimaryBS;
                        UE.primaryBSHistory = [UE.primaryBSHistory, UE.primaryBSIdx];
                        next_blockage = link{ue_idx,idxPrimaryBS}.discoveredTimes(3,find(link{ue_idx,idxPrimaryBS}.discoveredTimes(3,:)> currentTime,1));
                        UE.primaryNextEventTime(ue_idx) = next_blockage;
                        UE.primaryEventTimes = [UE.primaryEventTimes; currentTime];
                        UE.primaryEventIndices = [UE.primaryEventIndices; ue_idx];
                        UE.primaryEventDescriptions =  [UE.primaryEventDescriptions; {'UE '+string(ue_idx)+': MeasurementDone-Recovered-b4-BFDT-expired-No-Failure-ReconnectedWithoutDelays'}];
                    else
                        if UE.sub6ConnectionState(ue_idx) == 1
                            disp("Some problem")
                        else                        
                            sub6ConnectionState = UE.sub6ConnectionState;
                            sub6ConnectionState(ue_idx) = 1;
                            [channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW] = computePhysicalChannels_sub6_MIMO(params);
                            % rate_dl_before_handoff = compute_link_rates_MIMO(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,ue_idx,sub6ConnectionState);                                              
    %                         rate_dl_before_handoff = compute_link_rates_MIMOv2(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,ue_idx,sub6ConnectionState);                                         
    %                         rate_dl_before_handoff = compute_link_rates_MIMOv3(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,ue_idx,sub6ConnectionState);                                              
                            % rate_dl_before_handoff = compute_link_rates_MIMO_mmse(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,ue_idx,sub6ConnectionState);                                              
                            % rate_dl_before_handoff_old = compute_link_rates_MIMO_mmse(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,ue_idx,sub6ConnectionState);                                              
    %                         if (sub6ConnectionState == zeros(params.numUE,1))
    %                             p_fac = params.p_fac;
    %                             params.p_fac = 0;
    %                         end
    %                         rate_dl_before_handoff = compute_link_rates_MIMO_quadriga(params,link,ue_idx,sub6ConnectionState);    
    %                         if (sub6ConnectionState == zeros(params.numUE,1))
    %                             params.p_fac = p_fac;
    %                         end
    %                         D_old = params.D;
    %                         [params.D, ue_idxs_affected] = AP_reassign(params,ue_idx);
                            [~, ue_idxs_affected] = AP_reassign(params,ue_idx);
                            ue_rearranged = union(ue_idxs_affected, params.ue_rearranged);
                            rate_dl_before_handoff = compute_link_rates_MIMO_mmse(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,zeros(numUE,1));                                              
                            lb = quantile(rate_dl_before_handoff(ue_rearranged)./params.Band,params.loss_pc_thresh);
                            bw_alloc = Band - r_min_sub6/lb;
                            if (bw_alloc < 0) %|| isnan(bw_alloc)
                                bw_alloc = 0;
                                params.p_fac = 1;
                                rate_dl_after_handoff = rate_dl_before_handoff;
                            elseif isnan(bw_alloc)
                                rate_dl_after_handoff = compute_link_rates_MIMO_mmse(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,sub6ConnectionState);                                              
                                lb = quantile(rate_dl_after_handoff((1+numUE):(numUE+numUE_sub6)),params.loss_pc_thresh);                                
                            else 
                                params.ue_rearranged = ue_rearranged;
                                ues_not_affected = setdiff((1+numUE):(numUE+numUE_sub6),params.ue_rearranged);
                                params.scs_sub6(1) = bw_alloc;
                                params.scs_sub6(2) = Band - bw_alloc;
                                % user_sc_alloc = ones(numUE+numUE_sub6,params.num_sc_sub6);                               
                                user_sc_alloc = params.user_sc_alloc; %zeros(numUE+numUE_sub6,1);     
                                user_sc_alloc(find(sub6ConnectionState),1) = 1;
                                user_sc_alloc(find(sub6ConnectionState),2) = 0;
                                user_sc_alloc(ues_not_affected,1) = 1;
                                user_sc_alloc(ues_not_affected,2) = 1;
                                user_sc_alloc(params.ue_rearranged,1) = 0;
                                user_sc_alloc(params.ue_rearranged,2) = 1;
                                params.user_sc_alloc = user_sc_alloc;
                                ues_sharing = union(((1:numUE).*sub6ConnectionState),ues_not_affected);
        %                         rate_dl_after_handoff = compute_link_rates_MIMO(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,sub6ConnectionState);                                              
                                rate_dl_after_handoff = compute_link_rates_MIMOv4(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,sub6ConnectionState);    
                                ues_not_affected = setdiff((1+numUE):(numUE+numUE_sub6),params.ue_rearranged);
                                lb = quantile(rate_dl_after_handoff(ues_not_affected),params.loss_pc_thresh);
                            end
                            if (all(rate_dl_after_handoff(find((1:numUE)'.*sub6ConnectionState)) >= r_min) && (lb >= r_min_sub6))
                                UE.sub6ConnectionStarts = [UE.sub6ConnectionStarts, currentTime];
                                UE.sub6ConnectionStartIndices = [UE.sub6ConnectionStartIndices, ue_idx];
                                UE.sub6ConnectionState(ue_idx) = 1;
                                UE.sub6ConnectionStateHistory = [UE.sub6ConnectionStateHistory, UE.sub6ConnectionState];
                            end
                        end
                        %BS didnt recover now we declare beam faiure and go to
                        % either MCG recovery route or go to Beam Failure route
                        % If we have secondary BS in connected mode
                        if UE.secondaryConnectionState(ue_idx) == 1
                            UE.primaryConnectionState(ue_idx) = 5;
                            UE.primaryConnectionStateHistory = [UE.primaryConnectionStateHistory, UE.primaryConnectionState];
                            UE.primaryNextEventTime(ue_idx) = currentTime + UE.RACH_eff;
                            UE.tmpMCGBSIdx(ue_idx) = UE.secondaryBSIdx(ue_idx);
                            % In order to give this BS a priority in the future
                            lastConnectedBS = UE.primaryBSHistory(ue_idx,end);
                            bsLastConnectionTimes(ue_idx,lastConnectedBS) = currentTime;
                            [~,bsPriorities(ue_idx,:)] = sort(bsLastConnectionTimes(ue_idx,:),'descend');
                            % [~,bsPriorities(ue_idx,:)] = sort(bsLastConnectionTimes(ue_idx,1:numBS_mmW),'descend');
                            % now event recording
                            UE.primaryEventTimes = [UE.primaryEventTimes; currentTime];
                            UE.primaryEventIndices = [UE.primaryEventIndices; ue_idx];
                            UE.primaryEventDescriptions =  [UE.primaryEventDescriptions; {'UE '+string(ue_idx)+': MeasurementDone,SecondaryBS-Available-Initate-MCG-Fail-Recovery'}];
                        elseif UE.secondaryConnectionState(ue_idx) == 2
                            % Cancel the RACH of secondary as now we need RACH for
                            % primary. Secondary RACH already failed without
                            % assitance from primary.
                            UE.secondaryConnectionState(ue_idx) = 0;
                            UE.secondaryConnectionStateHistory = [UE.secondaryConnectionStateHistory, UE.secondaryConnectionState];
                            
                            UE.primaryConnectionState(ue_idx) = 2;
                            UE.primaryConnectionStateHistory = [UE.primaryConnectionStateHistory, UE.primaryConnectionState];
                            
                            UE.primaryTargetIdx(ue_idx) = UE.secondaryTargetIdx(ue_idx);
                            UE.secondaryTargetIdx(ue_idx) = 0;
                            UE.primaryBSIdx(ue_idx) = 0;
                            UE.primaryNextEventTime(ue_idx) = (UE.RACH_eff + currentTime);
                            UE.secondaryNextEventTime(ue_idx) = UE.primaryNextEventTime(ue_idx) + 1e-8;
                            UE.primaryEventTimes = [UE.primaryEventTimes; currentTime];
                            UE.primaryEventIndices = [UE.primaryEventIndices; ue_idx];
                            UE.primaryEventDescriptions = [UE.primaryEventDescriptions; {'UE '+string(ue_idx)+': MeasurementsDone,SecondaryRACH-Aborted,now-the-sameBS used as primary RACH.'}];
                        else
                            UE.primaryConnectionState(ue_idx) = 4;
                            UE.primaryConnectionStateHistory = [UE.primaryConnectionStateHistory, UE.primaryConnectionState];
                            UE.primaryNextEventTime(ue_idx) = currentTime + UE.beamFailureRecoveryTimer;
                            % In order to give this BS a priority in the future
                            lastConnectedBS = UE.primaryBSHistory(ue_idx,end);
                            bsLastConnectionTimes(ue_idx,lastConnectedBS) = currentTime;
                            [~,bsPriorities(ue_idx,:)] = sort(bsLastConnectionTimes(ue_idx,:),'descend');
                            % [~,bsPriorities(ue_idx,:)] = sort(bsLastConnectionTimes(ue_idx,1:numBS_mmW),'descend');
                            % now event recording
                            UE.primaryEventTimes = [UE.primaryEventTimes; currentTime];
                            UE.primaryEventIndices = [UE.primaryEventIndices; ue_idx];
                            UE.primaryEventDescriptions =  [UE.primaryEventDescriptions; {'UE '+string(ue_idx)+': MeasurementDone,-No-Secondary-BS-Available-4-MCG-Failure-Recovery,-Declared-Beam-Failure-going2-BeamRecovery'}];
                        end
                    end            
                end
            elseif UE.primaryConnectionState(ue_idx) == 4
                if UE.primaryNextEventTime(ue_idx) == currentTime
                    idxPrimaryBS = UE.primaryBSIdx(ue_idx);
                    % if isempty( idxPrimaryBS )
                    if idxPrimaryBS==0
                        isBSRecovered = false;
                    else
                        isBSRecovered = link{ue_idx,idxPrimaryBS}.blockageStatus;
                    end
                    if isBSRecovered
                        if link{ue_idx,idxPrimaryBS}.blockageStatus %&& link{idxPrimaryBS}.discovery_state
                            UE.primaryConnectionState(ue_idx) = 2;
                            UE.primaryConnectionStateHistory = [UE.primaryConnectionStateHistory, UE.primaryConnectionState];
                            UE.primaryTargetIdx(ue_idx) = idxPrimaryBS;
                            UE.primaryNextEventTime(ue_idx) = currentTime + UE.RACH_eff;
                            link{ue_idx,idxPrimaryBS}.nextEventTime = currentTime + UE.RACH_eff;
                            % now event recording
                            UE.primaryEventTimes = [UE.primaryEventTimes; currentTime];
                            UE.primaryEventIndices = [UE.primaryEventIndices; ue_idx];
                            UE.primaryEventDescriptions =  [UE.primaryEventDescriptions; {'UE '+string(ue_idx)+': Recovered-b4-RLF-expired-(Re)Connecting-to-BS'}];
                        else
                            disp('Arrived this line something might be wrong. DEBUG')
                        end
                    else
                        %Declare RLF to go idle state
                        UE.primaryConnectionState(ue_idx) = 0;
                        UE.primaryConnectionStateHistory = [UE.primaryConnectionStateHistory, UE.primaryConnectionState];
                        UE.primaryBSIdx(ue_idx) = 0;
                        % now event recording
                        UE.primaryEventTimes = [UE.primaryEventTimes; currentTime];
                        UE.primaryEventIndices = [UE.primaryEventIndices; ue_idx];
                        UE.primaryEventDescriptions =  [UE.primaryEventDescriptions; {'Declared-RLF-Going2-NonConnectedState-CanConnect-any-available-BS'}];
                    end
                end
            elseif UE.primaryConnectionState(ue_idx) == 5
                if UE.primaryNextEventTime(ue_idx) == currentTime
                    tmpMCGBSIdx = UE.tmpMCGBSIdx(ue_idx);
                    if link{ue_idx,tmpMCGBSIdx}.blockageStatus && link{ue_idx,tmpMCGBSIdx}.discovery_state
                        %Remove from secondary
                        UE.secondaryConnectionState(ue_idx) = 0;
                        UE.secondaryConnectionStateHistory = [UE.secondaryConnectionStateHistory, UE.secondaryConnectionState];
                        UE.secondaryBSIdx(ue_idx) = 0;
                        UE.secondaryConnectionEnds = [UE.secondaryConnectionEnds, currentTime];
                        UE.secondaryConnectionEndIndices = [UE.secondaryConnectionEndIndices, ue_idx];
                        UE.secondaryEventTimes = [UE.secondaryEventTimes; currentTime];
                        UE.secondaryEventIndices = [UE.secondaryEventIndices; ue_idx];
                        UE.secondaryEventDescriptions = [UE.secondaryEventDescriptions; {'UE '+string(ue_idx)+': primaryBlocked-and-Secondary-Changed-to-Primary'}];
                        % Add to primary
                        UE.primaryConnectionState(ue_idx) = 1;
                        UE.primaryConnectionStateHistory = [UE.primaryConnectionStateHistory, UE.primaryConnectionState];
                        UE.primaryConnectionStarts = [UE.primaryConnectionStarts, currentTime];
                        UE.primaryConnectionStartIndices = [UE.primaryConnectionStartIndices, ue_idx];
                        if UE.sub6ConnectionState(ue_idx) == 1
                            UE.sub6ConnectionEnds = [UE.sub6ConnectionEnds, currentTime];
                            UE.sub6ConnectionEndIndices = [UE.sub6ConnectionEndIndices, ue_idx];
                            UE.sub6ConnectionState(ue_idx) = 0;
                            UE.sub6ConnectionStateHistory = [UE.sub6ConnectionStateHistory, UE.sub6ConnectionState];
                        end
                        UE.primaryBSIdx(ue_idx) = tmpMCGBSIdx;
                        UE.primaryBSHistory = [UE.primaryBSHistory, UE.primaryBSIdx];
                        next_blockage = link{ue_idx,tmpMCGBSIdx}.discoveredTimes(3,find(link{ue_idx,tmpMCGBSIdx}.discoveredTimes(3,:)> currentTime,1));
                        UE.primaryNextEventTime(ue_idx) = next_blockage;
                        UE.primaryEventTimes = [UE.primaryEventTimes; currentTime];
                        UE.primaryEventIndices = [UE.primaryEventIndices; ue_idx];
                        UE.primaryEventDescriptions =  [UE.primaryEventDescriptions; {'UE '+string(ue_idx)+': primaryGotBlocked-Secondary-Become-Primary-BS-Completed'}];
                    else
                        %Secondary BS also have problems, so try recovery with the
                        %primary
                        UE.primaryConnectionState(ue_idx) = 4;
                        UE.primaryConnectionStateHistory = [UE.primaryConnectionStateHistory, UE.primaryConnectionState];
                        UE.primaryNextEventTime(ue_idx) = currentTime + UE.beamFailureRecoveryTimer;
                        % In order to give this BS a priority in the future
                        lastConnectedBS = UE.primaryBSHistory(ue_idx,end);
                        bsLastConnectionTimes(ue_idx,lastConnectedBS) = currentTime;
                        [~,bsPriorities(ue_idx,:)] = sort(bsLastConnectionTimes(ue_idx,:),'descend');
                        % [~,bsPriorities(ue_idx,:)] = sort(bsLastConnectionTimes(ue_idx,1:numBS_mmW),'descend');
                        % now event recording
                        UE.primaryEventTimes = [UE.primaryEventTimes; currentTime];
                        UE.primaryEventIndices = [UE.primaryEventIndices; ue_idx];
                        UE.primaryEventDescriptions =  [UE.primaryEventDescriptions; {'UE '+string(ue_idx)+': Secondary BS also have problems, couldnt transfer secondary to primary. Trying BFR with primary.'}];
                    end
                end
            end
        end
    end
    %Now conection establishment start or end procedures
    %Let's start with primary
    for ue_idx = 1:numUE
        if (numBS_mobile(ue_idx) > 0)
            if UE.primaryConnectionState(ue_idx) == 0
                % UE = tryPrimaryConnecting(UE,currentTime,link,bsPriorities,bsLastConnectionTimes,ue_idx);
                UE = tryPrimaryConnecting(UE,currentTime,link,bsPriorities,bsLastConnectionTimes,ue_idx,numBS_mobile);
            elseif UE.primaryConnectionState(ue_idx) == 2
                if UE.primaryNextEventTime(ue_idx) == currentTime
                    idxBS = UE.primaryTargetIdx(ue_idx);
                    if link{ue_idx,idxBS}.discovery_state
                        UE.primaryConnectionState(ue_idx) = 1;
                        UE.primaryConnectionStateHistory = [UE.primaryConnectionStateHistory, UE.primaryConnectionState];
                        UE.primaryTargetIdx(ue_idx) = 0;
                        UE.primaryConnectionStarts = [UE.primaryConnectionStarts, currentTime];
                        UE.primaryConnectionStartIndices = [UE.primaryConnectionStartIndices, ue_idx];
                        if UE.sub6ConnectionState(ue_idx) == 1
                            UE.sub6ConnectionEnds = [UE.sub6ConnectionEnds, currentTime];
                            UE.sub6ConnectionEndIndices = [UE.sub6ConnectionEndIndices, ue_idx];
                            UE.sub6ConnectionState(ue_idx) = 0;
                            UE.sub6ConnectionStateHistory = [UE.sub6ConnectionStateHistory, UE.sub6ConnectionState];
                        end
                        UE.primaryBSIdx(ue_idx) = idxBS;
                        UE.primaryBSHistory = [UE.primaryBSHistory, UE.primaryBSIdx];
                        next_blockage = link{ue_idx,idxBS}.discoveredTimes(3,find(link{ue_idx,idxBS}.discoveredTimes(3,:)> currentTime,1));
                        UE.primaryNextEventTime(ue_idx) = next_blockage;
                        UE.primaryEventTimes = [UE.primaryEventTimes; currentTime];
                        UE.primaryEventIndices = [UE.primaryEventIndices; ue_idx];
                        UE.primaryEventDescriptions =  [UE.primaryEventDescriptions; {'UE '+string(ue_idx)+': primaryRachCompleted-Connected'}];
                    else
                        UE.primaryConnectionState(ue_idx) = 0;
                        UE.primaryConnectionStateHistory = [UE.primaryConnectionStateHistory, UE.primaryConnectionState];
                        UE.primaryEventTimes = [UE.primaryEventTimes; currentTime];
                        UE.primaryEventIndices = [UE.primaryEventIndices; ue_idx];
                        UE.primaryEventDescriptions =  [UE.primaryEventDescriptions; {'UE '+string(ue_idx)+': primaryRachFailed-NotConnected-go2-Idle-State-Can-Try-another-BS'}];
                        % UE = tryPrimaryConnecting(UE,currentTime,link,bsPriorities,bsLastConnectionTimes,ue_idx);
                        UE = tryPrimaryConnecting(UE,currentTime,link,bsPriorities,bsLastConnectionTimes,ue_idx,numBS_mobile);
                    end
                end
            end
        end
    end

        %Secondary Connection Establishment Procedures
    for ue_idx = 1:numUE
        if (numBS_mobile(ue_idx) > 0)
            if UE.secondaryConnectionState(ue_idx) == 0
                % UE = trySecondaryConnecting(UE,currentTime,link,bsPriorities,bsLastConnectionTimes,ue_idx);
                UE = trySecondaryConnecting(UE,currentTime,link,bsPriorities,bsLastConnectionTimes,ue_idx,numBS_mobile);
            elseif UE.secondaryConnectionState(ue_idx) == 2
                if UE.secondaryNextEventTime(ue_idx) == currentTime
                    idxBS = UE.secondaryTargetIdx(ue_idx);
                    if link{ue_idx,idxBS}.discovery_state
                        UE.secondaryConnectionState(ue_idx) = 1;
                        UE.secondaryConnectionStateHistory = [UE.secondaryConnectionStateHistory, UE.secondaryConnectionState];
                        UE.secondaryTargetIdx(ue_idx) = 0;
                        UE.secondaryConnectionStarts = [UE.secondaryConnectionStarts, currentTime];
                        UE.secondaryConnectionStartIndices = [UE.secondaryConnectionStartIndices, ue_idx];
                        UE.secondaryBSIdx(ue_idx) = idxBS;
                        UE.secondaryBSHistory = [UE.secondaryBSHistory, UE.secondaryBSIdx];
                        next_blockage = link{ue_idx,idxBS}.discoveredTimes(3,find(link{ue_idx,idxBS}.discoveredTimes(3,:)> currentTime,1));
                        UE.secondaryNextEventTime(ue_idx) = next_blockage;
                        UE.secondaryEventTimes = [UE.secondaryEventTimes; currentTime];
                        UE.secondaryEventIndices = [UE.secondaryEventIndices; ue_idx];
                        UE.secondaryEventDescriptions =  [UE.secondaryEventDescriptions; {'UE '+string(ue_idx)+': secondaryRachCompleted-Connected'}];
                    else
                        UE.secondaryConnectionState(ue_idx) = 0;
                        UE.secondaryConnectionStateHistory = [UE.secondaryConnectionStateHistory, UE.secondaryConnectionState];
                        UE.secondaryEventTimes = [UE.secondaryEventTimes; currentTime];
                        UE.secondaryEventIndices = [UE.secondaryEventIndices; ue_idx];
                        UE.secondaryEventDescriptions =  [UE.secondaryEventDescriptions; {'UE '+string(ue_idx)+': secondaryRachFailed-NotConnected-go2-Idle-State-Can-Try-another-BS'}];
                        % UE = trySecondaryConnecting(UE,currentTime,link,bsPriorities,bsLastConnectionTimes,ue_idx);
                        UE = trySecondaryConnecting(UE,currentTime,link,bsPriorities,bsLastConnectionTimes,ue_idx,numBS_mobile);
                    end
                end
            end
        end
    end
    
    % nextEventTimes = getNextEvents(link,params.simTime);
    nextEventTimes = getNextEvents(link,params.simTime,numBS_mobile);
%     link_mat = cell2mat(link);
%     nextEventTimes = getNextEvents(mat2cell(link_mat(1:numUE,1:params.numGNB),[1,1]),params.simTime);
    nextEventTime = params.simTime;
    for ue_idx = 1:numUE
        % if (min(nextEventTimes(ue_idx,1:params.numGNB)) < nextEventTime)
        if (min(nextEventTimes(ue_idx,1:numBS_mobile(ue_idx))) < nextEventTime)
            % nextEventTime = min(nextEventTimes(ue_idx,1:params.numGNB));
            nextEventTime = min(nextEventTimes(ue_idx,1:numBS_mobile(ue_idx)));
        end
    end
    %need to check Ue event times as well
    UEEventTimes = [UE.primaryNextEventTime;UE.secondaryNextEventTime];
    if min(UEEventTimes(UEEventTimes > currentTime)) > currentTime
        UEnextEventTime = min(UEEventTimes(UEEventTimes > currentTime));
    else
        UEnextEventTime = [];
    end
    
    if UEnextEventTime < nextEventTime
        nextEventTime = UEnextEventTime;
    end
end

% Just to have complete events, if the UE is connected at the end of 
% simulationTime it will not have an connection lost event. So manually add
% that recording event.
for ue_idx = 1:numUE
    if length(UE.primaryConnectionStarts(UE.primaryConnectionStartIndices==ue_idx)) > length(UE.primaryConnectionEnds(UE.primaryConnectionEndIndices==ue_idx))
        UE.primaryConnectionEnds = [UE.primaryConnectionEnds, params.simTime];
        UE.primaryConnectionEndIndices = [UE.primaryConnectionEndIndices, ue_idx];
    end
end

for ue_idx = 1:numUE
    if length(UE.secondaryConnectionStarts(UE.secondaryConnectionStartIndices==ue_idx)) > length(UE.secondaryConnectionEnds(UE.secondaryConnectionEndIndices==ue_idx))
        UE.secondaryConnectionEnds = [UE.secondaryConnectionEnds, params.simTime];
        UE.secondaryConnectionEndIndices = [UE.secondaryConnectionEndIndices, ue_idx];
    end
end

for ue_idx = 1:numUE
    if length(UE.sub6ConnectionStarts(UE.sub6ConnectionStartIndices==ue_idx)) > length(UE.sub6ConnectionEnds(UE.sub6ConnectionEndIndices==ue_idx))
        UE.sub6ConnectionEnds = [UE.sub6ConnectionEnds, params.simTime];
        UE.sub6ConnectionEndIndices = [UE.sub6ConnectionEndIndices, ue_idx];
    end
end

%% Event structuring
connectionStarts = [];
connectionStartIndices = [];
connectionEnds = [];
connectionEndIndices = [];
sub6connectionStarts = [];
sub6connectionStartIndices = [];
sub6connectionEnds = [];
sub6connectionEndIndices = [];
connectionStarts = [connectionStarts, UE.primaryConnectionStarts];
connectionStartIndices = [connectionStartIndices, UE.primaryConnectionStartIndices];
connectionEnds = [connectionEnds, UE.primaryConnectionEnds];
connectionEndIndices = [connectionEndIndices, UE.primaryConnectionEndIndices];
sub6connectionStarts = [sub6connectionStarts, UE.sub6ConnectionStarts];
sub6connectionStartIndices = [sub6connectionStartIndices, UE.sub6ConnectionStartIndices];
sub6connectionEnds = [sub6connectionEnds, UE.sub6ConnectionEnds];
sub6connectionEndIndices = [sub6connectionEndIndices, UE.sub6ConnectionEndIndices];
outage_probability_wo_cf = ones(numUE,1);
mean_outage_duration_wo_cf = zeros(numUE,1);
outage_probability = ones(numUE,1);
mean_outage_duration = zeros(numUE,1);
for ue_idx = 1:numUE
    connectionEvents = [connectionStarts(connectionStartIndices==ue_idx);connectionEnds(connectionEndIndices==ue_idx)-connectionStarts(connectionStartIndices==ue_idx);connectionEnds(connectionEndIndices==ue_idx)];
    sub6connectionEvents = [sub6connectionStarts(sub6connectionStartIndices==ue_idx);sub6connectionEnds(sub6connectionEndIndices==ue_idx)-sub6connectionStarts(sub6connectionStartIndices==ue_idx);sub6connectionEnds(sub6connectionEndIndices==ue_idx)];
    %Merge connection events so that means track the times at
    %least one of the LAs is connected. If there is not even
    %one that means outage.
    connectionEvents = mergeConnectionEvents(connectionEvents);
    sub6connectionEvents = mergeConnectionEvents(sub6connectionEvents);
    outageEvents = getOutageEvents(connectionEvents,params);
    outage_durations_wo_cf = outageEvents(2,:);
    outage_duration_wo_cf = sum(outageEvents(2,:));
    connected_duration_wo_cf = sum(connectionEvents(2,:));
    try
%       outage_not_mitigated_by_cf = (setdiff(outageEvents',sub6connectionEvents','rows'))';
        [~,ia] = setdiff(outageEvents(1,:),sub6connectionEvents(1,:)-1e-8);
        outage_durations_wi_cf = outageEvents(2,ia);
        outage_duration = sum(outageEvents(2,:)) - sum(sub6connectionEvents(2,:));
        connected_duration = sum(connectionEvents(2,:)) + sum(sub6connectionEvents(2,:));
    catch 
        outage_durations_wi_cf = outage_durations_wo_cf;
        outage_duration = sum(outageEvents(2,:));
        connected_duration = sum(connectionEvents(2,:));
    end
    total_duration =   outage_duration + connected_duration;
    total_duration_wo_cf =   outage_duration_wo_cf + connected_duration_wo_cf;
    outage_probability_wo_cf(ue_idx) = outage_duration_wo_cf / total_duration;
    mean_outage_duration_wo_cf(ue_idx) = outage_duration_wo_cf / size(outageEvents,2);
    outage_probability(ue_idx) = outage_duration / total_duration;
    mean_outage_duration(ue_idx) = outage_duration / size(outageEvents,2);
    if abs(total_duration - params.simTime) > 1e-10
        warning('TotalTime and simTime doesnt match check here.')
    end
    if abs(total_duration_wo_cf - params.simTime) > 1e-10
        warning('TotalTime and simTime doesnt match check here.')
    end
end

%% OUTPUT preparation
%Creating the output struct

% Main results are hold in lookAngle structs.
simOutputs.UE = UE;

% General Event Metrics
simOutputs.connectionEvents = connectionEvents;
simOutputs.outageEvents = outageEvents;
simOutputs.outage_duration = outage_duration;
simOutputs.connected_duration = connected_duration;
simOutputs.total_duration =   total_duration;
simOutputs.outage_probability_wo_cf = outage_probability_wo_cf;
simOutputs.mean_outage_duration_wo_cf = mean_outage_duration_wo_cf;
simOutputs.outage_probability = outage_probability;
simOutputs.mean_outage_duration = mean_outage_duration;
simOutputs.outage_durations_wo_cf = outage_durations_wo_cf;
simOutputs.outage_durations_wi_cf = outage_durations_wi_cf;
% General parameters
simOutputs.discovery_delay = discovery_delay;
simOutputs.failureDetectionDelay = failureDetectionDelay;
simOutputs.connection_setup_delay = connection_setup_delay;
simOutputs.signalingAfterRachDelay = signalingAfterRachDelay;
simOutputs.beamFailureRecoveryTimer = beamFailureRecoveryTimer;
simOutputs.RACH_eff = RACH_eff;
simOutputs.params = params;
simOutputs.protocolParams = protocolParams;
simOutputs.frameHopCount = protocolParams.frameHopCount;
simOutputs.frameDeliveryDelay = protocolParams.frameDeliveryDelay;

% fprintf('numBS: %d. \n',numBS)


fprintf('Simulation done : %f seconds\n',toc)
end

