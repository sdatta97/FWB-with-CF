function connectionEvents = mergeConnectionEvents(connectionEvents)
try 
    [~,sort_idx]=sort(connectionEvents(1,:));
catch ME
    return
end
connectionEvents = connectionEvents(:,sort_idx);
%Merging
isCleared = 0;
while ~isCleared
    isCleared=1;
    len = size(connectionEvents,2);
    for jj=len:-1:2
        if connectionEvents(3,jj-1) >= connectionEvents(1,jj)
            isCleared=0;
            connectionEvents(3,jj-1) = max(connectionEvents(3,jj),connectionEvents(3,jj-1));
            connectionEvents(:,jj) = [];
            connectionEvents(2,jj-1) = connectionEvents(3,jj-1) - connectionEvents(1,jj-1);
        end
    end
end
connectionEvents(2,:) = connectionEvents(3,:)-connectionEvents(1,:);
end

