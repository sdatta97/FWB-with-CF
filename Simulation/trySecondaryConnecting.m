% function UE = trySecondaryConnecting(UE,currentTime,link,bsPriorities,bsLastConnectionTimes)
function UE = trySecondaryConnecting(UE,currentTime,link,bsPriorities,bsLastConnectionTimes,ue_idx)
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

% find discovered BSs for this UE,
% if there is enough time to initialize connection:
% initialize connection to the highest priority BS
% else
% check next priority bs etc...
% start connection initialization, set nextEventTime
% for idxBS = bsPriorities
numBS = UE.numGNB;
for i = 1:numBS
    % idxBS = bsPriorities(i);
    idxBS = bsPriorities(ue_idx,i);
    % if (ismember(idxBS, UE.primaryBSIdx)) || (ismember(idxBS, UE.primaryTargetIdx))
    if (ismember(idxBS, UE.primaryBSIdx(ue_idx))) || (ismember(idxBS, UE.primaryTargetIdx(ue_idx)))
        continue
        %We skip this BS since its either connected with the UE as primary
        %or is establishing connection with the UE as a primary.
    end
    % if UE.primaryConnectionState ~= 1
    %     UE.secondaryNextEventTime = UE.primaryNextEventTime + 1e-8;
    %     UE.secondaryEventTimes = [UE.secondaryEventTimes; currentTime];
    %     UE.secondaryEventDescriptions = [UE.secondaryEventDescriptions; ...
    %                 {'Idle-(No-Active-Primary-Connection-Available-Wait-Until-Primary-Established)'}];
    %     return
    % end
    if UE.primaryConnectionState(ue_idx) ~= 1
        UE.secondaryNextEventTime(ue_idx) = UE.primaryNextEventTime(ue_idx) + 1e-8;
        UE.secondaryEventTimes = [UE.secondaryEventTimes; currentTime];
        UE.secondaryEventIndices = [UE.secondaryEventIndices; ue_idx];
        UE.secondaryEventDescriptions = [UE.secondaryEventDescriptions; ...
                    {'UE '+string(ue_idx)+': Idle-(No-Active-Primary-Connection-Available-Wait-Until-Primary-Established)'}];
        return
    end
    % if link{idxBS}.discovery_state
    %     UE.secondaryConnectionState = 2;
    %     UE.secondaryConnectionStateHistory = [UE.secondaryConnectionStateHistory; UE.secondaryConnectionState];
    %     UE.secondaryTargetIdx = idxBS;
    %     UE.secondaryBSIdx = [];
    %     UE.secondaryNextEventTime = (UE.RACH_eff + currentTime);
    %     UE.secondaryEventTimes = [UE.secondaryEventTimes; currentTime];
    %     UE.secondaryEventDescriptions = [UE.secondaryEventDescriptions; {'InitiatingRach-Secondary(RACH)'}];
    %     return %no need to visit other BSs, try to connect most fresh one
    % end
    if link{ue_idx,idxBS}.discovery_state
        UE.secondaryConnectionState(ue_idx) = 2;
        UE.secondaryConnectionStateHistory = [UE.secondaryConnectionStateHistory, UE.secondaryConnectionState];
        UE.secondaryTargetIdx(ue_idx) = idxBS;
        UE.secondaryBSIdx(ue_idx) = 0;
        UE.secondaryNextEventTime(ue_idx) = (UE.RACH_eff + currentTime);
        UE.secondaryEventTimes = [UE.secondaryEventTimes; currentTime];
        UE.secondaryEventIndices = [UE.secondaryEventIndices; ue_idx];
        UE.secondaryEventDescriptions = [UE.secondaryEventDescriptions; {'UE '+string(ue_idx)+': InitiatingRach-Secondary(RACH)'}];
        return %no need to visit other BSs, try to connect most fresh one
    end
end
UE.secondaryEventTimes = [UE.secondaryEventTimes; currentTime];
UE.secondaryEventIndices = [UE.secondaryEventIndices; currentTime];
UE.secondaryEventDescriptions = [UE.secondaryEventDescriptions; ...
                    {'UE '+string(ue_idx)+': Idle-(No-Unblocked&Discovered-BS-Available)(Maybe one avail. is used as primary)'}];
end

