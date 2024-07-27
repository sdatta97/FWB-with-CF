% function status = discoveryStatus(discoveredTimes,idxBS,currentTime, ue_idx, numBS)
function status = discoveryStatus(discoveredTimes,idxBS,currentTime, ue_idx, numBS_mobile)

% if sum((discoveredTimes{idxBS}(1,:) <= currentTime) & (discoveredTimes{idxBS}(3,:) > currentTime))
% if sum((discoveredTimes{(ue_idx-1)*numBS + idxBS}(1,:) <= currentTime) & (discoveredTimes{(ue_idx-1)*numBS + idxBS}(3,:) > currentTime))
if sum((discoveredTimes{sum(numBS_mobile(1:(ue_idx-1)))+ idxBS}(1,:) <= currentTime) & (discoveredTimes{sum(numBS_mobile(1:(ue_idx-1))) + idxBS}(3,:) > currentTime))
    status = 1; %discovered
else
    status = 0; %not discovered
end

end


