function outageEvents = getOutageEvents(connectionEvents,params)
try 
    outage_starts = [0,connectionEvents(3,:)];
    outage_ends = [connectionEvents(1,:),params.simTime];
catch ME
    return
end
outageEvents = [outage_starts; outage_ends-outage_starts ;outage_ends];
if outageEvents(3,1) == 0
    outageEvents(:,1) = [];
end
if outageEvents(1,end) == params.simTime
    outageEvents(:,end)=[];
end

end

