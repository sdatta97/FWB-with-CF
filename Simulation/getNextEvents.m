function nextEvents = getNextEvents(link,simTime)
% nextEvents = size(link);
% for idxBS = 1:size(link,1)
%     nextEvents(idxBS) = min([simTime,link{idxBS}.nextEventTime]);
% end
% nextEvents = size(link);
nextEvents = zeros(size(link));
for ue_idx = 1:size(link,1)
    for idxBS = 1:size(link,2)
        nextEvents(ue_idx,idxBS) = min([simTime,link{ue_idx,idxBS}.nextEventTime]);
    end
end
end
