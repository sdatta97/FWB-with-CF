% function link = updatePhysicalDiscovery(currentTime,link,discoveredTimes,bsBlockageTimes,ue_idx)
function link = updatePhysicalDiscovery(currentTime,link,discoveredTimes,bsBlockageTimes)
%Physical + Discovery Updates of all links, if this event time is realted
%to that link. (Multiple links can be updated at once, for example, if the
%event is a rotation event.)
% numBS = size(link,1);
numUE = size(link,1);
numBS = size(link,2);
for ue_idx = 1:numUE
    for idxBS = 1:numBS
        % link{idxBS}.discovery_state = discoveryStatus(discoveredTimes,idxBS,currentTime);
        % link{idxBS}.blockageStatus = blockageStatus(bsBlockageTimes,idxBS,currentTime);
        link{ue_idx,idxBS}.discovery_state = discoveryStatus(discoveredTimes,idxBS,currentTime,ue_idx,numBS);
        link{ue_idx,idxBS}.blockageStatus = blockageStatus(bsBlockageTimes,idxBS,currentTime,ue_idx,numBS);
        if link{ue_idx,idxBS}.nextEventTime ~= currentTime
            continue
        end
        if link{ue_idx,idxBS}.discovery_state %link isdiscovered
            %check discovery status (initially non-discovered)
            next_blockage = link{ue_idx,idxBS}.discoveredTimes(3,find(link{ue_idx,idxBS}.discoveredTimes(3,:)> currentTime,1));
            link{ue_idx,idxBS}.nextEventTime = next_blockage;
        else
            %link is not discovered find the next discovery event, thats the next
            %event time for this link.
            next_discovery = link{ue_idx,idxBS}.discoveredTimes(1,find(link{ue_idx,idxBS}.discoveredTimes(1,:)> currentTime,1));
            link{ue_idx,idxBS}.nextEventTime = next_discovery;
        end
    end
end
end

