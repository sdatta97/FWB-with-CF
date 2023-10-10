function [discoveredTimes,bsBlockageTimes] = computeDiscoveredTimes(dataBS_mobile,params,discovery_time,failureDetectionDelay)
discoveredTimes = cell(size(dataBS_mobile));
bsBlockageTimes = cell(size(dataBS_mobile));
for ue_idx = 1:params.numUE
% for ue_idx = 1:(params.numUE+params.numUE_sub6)
%     for bs_idx = 1:(params.numGNB+params.numGNB_sub6)
    for bs_idx = 1:params.numGNB  %size(dataBS_mobile,1)
    % for bs_idx = 1:params.numGNB_sub6  %size(dataBS_mobile,1)
        % ins = [0,dataBS_mobile{bs_idx}(3,1:end)];
        % outs = [dataBS_mobile{bs_idx}(1,1:end),params.simTime];
        % ins = [0,dataBS_mobile{(ue_idx-1)*params.numGNB_sub6+bs_idx}(3,1:end)];
        % outs = [dataBS_mobile{(ue_idx-1)*params.numGNB_sub6+bs_idx}(1,1:end),params.simTime];
        ins = [0,dataBS_mobile{(ue_idx-1)*params.numGNB+bs_idx}(3,1:end)];
        outs = [dataBS_mobile{(ue_idx-1)*params.numGNB+bs_idx}(1,1:end),params.simTime];    
        % ins = [0,dataBS_mobile{(ue_idx-1)*(params.numGNB+params.numGNB_sub6)+bs_idx}(3,1:end)];
        % outs = [dataBS_mobile{(ue_idx-1)*(params.numGNB+params.numGNB_sub6)+bs_idx}(1,1:end),params.simTime];  
        durs = outs - ins;
        BS_discovered_times = [ins;durs;outs];
        [~,sort_idx] = sort(BS_discovered_times(1,:));
        BS_discovered_times = BS_discovered_times(:,sort_idx);
        %First clear overlapping non-blocked-times
        isCleared = 0;
        while ~isCleared
            isCleared=1;
            len = size(BS_discovered_times,2);
            for jj=len:-1:2
                if BS_discovered_times(3,jj-1)  >= BS_discovered_times(1,jj)
                    isCleared=0;
                    BS_discovered_times(3,jj-1) = max(BS_discovered_times(3,jj),BS_discovered_times(3,jj-1));
                    BS_discovered_times(:,jj) = [];
                    BS_discovered_times(2,jj-1) = BS_discovered_times(3,jj-1) - BS_discovered_times(1,jj-1);
                end
            end  
        end
        if BS_discovered_times(3,1) ==0
            BS_discovered_times(:,1)=[];
        end
        if BS_discovered_times(1,end) == params.simTime
            BS_discovered_times(:,end)=[];
        end
        %Now impose the discovery times when necessary, by shifting t_in times by
        %discovery duration
        % bsBlockageTimes{bs_idx} = BS_discovered_times;
        bsBlockageTimes{(ue_idx-1)*params.numGNB+bs_idx} = BS_discovered_times;        
        % bsBlockageTimes{(ue_idx-1)*params.numGNB_sub6+bs_idx} = BS_discovered_times;    
        % bsBlockageTimes{(ue_idx-1)*(params.numGNB+params.numGNB_sub6)+bs_idx} = BS_discovered_times;        
        len = size(BS_discovered_times,2);
        for jj=len:-1:1
            if (BS_discovered_times(1,jj) + discovery_time) < BS_discovered_times(3,jj) 
                BS_discovered_times(1,jj) = BS_discovered_times(1,jj) + discovery_time;
                BS_discovered_times(2,jj) = BS_discovered_times(3,jj) - BS_discovered_times(1,jj);
            else
                BS_discovered_times(:,jj) = [];
            end
        end  
%         discoveredTimes{(ue_idx-1)*(params.numGNB+params.numGNB_sub6)+bs_idx} = BS_discovered_times;    
        % discoveredTimes{bs_idx} = BS_discovered_times;    
        discoveredTimes{(ue_idx-1)*params.numGNB+bs_idx} = BS_discovered_times;    
        % discoveredTimes{(ue_idx-1)*params.numGNB_sub6+bs_idx} = BS_discovered_times;    
    end
end
end

