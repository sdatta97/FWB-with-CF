function status = blockageStatus(bsBlockageTimes,idxBS,currentTime, ue_idx, numBS)

% if sum((bsBlockageTimes{idxBS}(1,:) <= currentTime) & (bsBlockageTimes{idxBS}(3,:) > currentTime))
if sum((bsBlockageTimes{(ue_idx-1)*numBS + idxBS}(1,:) <= currentTime) & (bsBlockageTimes{(ue_idx-1)*numBS + idxBS}(3,:) > currentTime))
    status = 1; %not blocked
else
    status = 0; %blocked
end

end


