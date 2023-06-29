% function UE = tryPrimaryConnecting(UE,currentTime,link,bsPriorities,bsLastConnectionTimes)
function UE = tryPrimaryConnecting(UE,currentTime,link,bsPriorities,bsLastConnectionTimes,ue_idx)
%TRYCONNECTING If a UE is not connected we try all the BS to find
%a target bs to establish connection

% % 
% % %Primary State variables
% % UE.primaryConnectionState = 0;
% % UE.primaryBSIdx = [];
% % UE.primaryTargetIdx = [];
% % UE.primaryNextEventTime = -100;

% % %Secondary State variables
% % UE.secondaryConnectionState = 0;
% % UE.secondaryBSIdx = [];
% % UE.secondaryTargetIdx = [];
% % UE.secondaryNextEventTime = -100;


% % % UE.secondaryConnectionStarts = [];
% % % UE.secondaryConnectionEnds = [];
% % % UE.secondaryBSHistory = [];
% % % UE.secondaryEventTimes = [];
% % % UE.secondaryEventDescriptions = [];

% Connection establishing for primary.
% Check if secondary connection is available
% if secondary connection available
%     start MCG failure and make the secondary as the primary
% else
%     check any available BS.
%
numBS = UE.numGNB;
% for idxBS = bsPriorities
for i = 1:numBS
   % idxBS = bsPriorities(i);
   idxBS = bsPriorities(ue_idx,i);
   % if link{idxBS}.discovery_state
    if link{ue_idx,idxBS}.discovery_state
       % if (ismember(idxBS, UE.secondaryBSIdx)) 
       if (ismember(idxBS, UE.secondaryBSIdx(ue_idx))) 
            %Takeover Secondary:
            % UE.primaryConnectionState = 5;
            % UE.tmpMCGBSIdx = UE.secondaryBSIdx;
            % UE.primaryConnectionStateHistory = [UE.primaryConnectionStateHistory; UE.primaryConnectionState];
            % UE.primaryNextEventTime = currentTime + UE.RACH_eff;
            % % now event recording
            % UE.primaryEventTimes = [UE.primaryEventTimes; currentTime];
            % UE.primaryEventDescriptions =  [UE.primaryEventDescriptions; {'PrimaryWasIdle-FoundBS-,SecondaryBS-Available-Initate-MCG-Fail-Recovery'}];
            UE.primaryConnectionState (ue_idx) = 5;
            UE.tmpMCGBSIdx(ue_idx) = UE.secondaryBSIdx(ue_idx);
            UE.primaryConnectionStateHistory = [UE.primaryConnectionStateHistory, UE.primaryConnectionState];
            UE.primaryNextEventTime(ue_idx) = currentTime + UE.RACH_eff;
            % now event recording
            UE.primaryEventTimes = [UE.primaryEventTimes; currentTime];
            UE.primaryEventIndices = [UE.primaryEventIndices; ue_idx];
            UE.primaryEventDescriptions =  [UE.primaryEventDescriptions; {'UE '+string(ue_idx)+': PrimaryWasIdle-FoundBS-,SecondaryBS-Available-Initate-MCG-Fail-Recovery'}];
            return
        end
        % if (ismember(idxBS, UE.secondaryTargetIdx))
        if (ismember(idxBS, UE.secondaryTargetIdx(ue_idx)))
            %We skip this BS since its establishing connection with the UE as a secondary.
            % Cancel the RACH of secondary as now we need RACH for
            % primary. Secondary RACH already failed without
            % assitance from primary.
            % UE.secondaryConnectionState = 0;
            % UE.secondaryConnectionStateHistory = [UE.secondaryConnectionStateHistory; UE.secondaryConnectionState];
            % 
            % UE.primaryConnectionState = 2;
            % UE.primaryConnectionStateHistory = [UE.primaryConnectionStateHistory; UE.primaryConnectionState];
            % 
            % UE.primaryTargetIdx = UE.secondaryTargetIdx;
            % UE.secondaryTargetIdx = [];
            % UE.primaryBSIdx = [];
            % UE.primaryNextEventTime = (UE.RACH_eff + currentTime);
            % UE.secondaryNextEventTime = UE.primaryNextEventTime + 1e-8;
            % UE.primaryEventTimes = [UE.primaryEventTimes; currentTime];
            % UE.primaryEventDescriptions = [UE.primaryEventDescriptions; {'MeasurementsDone,SecondaryRACHFAiled,now-the-sameBS used as primary RACH.'}];
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
            UE.primaryEventDescriptions = [UE.primaryEventDescriptions; {'UE '+string(ue_idx)+': MeasurementsDone,SecondaryRACHFAiled,now-the-sameBS used as primary RACH.'}];
            return    
        end
        % UE.primaryConnectionState = 2;
        % UE.primaryConnectionStateHistory = [UE.primaryConnectionStateHistory; UE.primaryConnectionState];
        % UE.primaryTargetIdx = idxBS;
        % UE.primaryBSIdx = [];
        % UE.primaryNextEventTime = (UE.RACH_eff + currentTime);
        % UE.primaryEventTimes = [UE.primaryEventTimes; currentTime];
        UE.primaryConnectionState(ue_idx) = 2;
        UE.primaryConnectionStateHistory = [UE.primaryConnectionStateHistory, UE.primaryConnectionState];
        UE.primaryTargetIdx(ue_idx) = idxBS;
        UE.primaryBSIdx(ue_idx) = 0;
        UE.primaryNextEventTime(ue_idx) = (UE.RACH_eff + currentTime);
        UE.primaryEventTimes = [UE.primaryEventTimes; currentTime];
        UE.primaryEventIndices = [UE.primaryEventIndices; ue_idx];
        UE.primaryEventDescriptions = [UE.primaryEventDescriptions; {'UE '+string(ue_idx)+': InitiatingRach-Primary(RACH)'}];
        return %no need to visit other BSs, try to connect most fresh one
    end
end

UE.primaryEventTimes = [UE.primaryEventTimes; currentTime];
UE.primaryEventIndices = [UE.primaryEventIndices; ue_idx];
UE.primaryEventDescriptions = [UE.primaryEventDescriptions; ...
                    {'UE '+string(ue_idx)+': Idle-(No-Unblocked&Discovered-BS-Available)'}];
end

