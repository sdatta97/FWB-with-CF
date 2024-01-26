function dataBS = computeBlockageEvents(params,ue_idx)
%COMPUTEBLOCKAGEEVENTS Computes blockage arrival and blockage departure
%times for every BS. Physical characteristics of the system independent of
%the protocol.
%   Finds dataBS cell where each elemnt of cell represents one bs
%   for each cell 3 by numBlock matrix:
%       first row: blockage arrival time
%       second row: blockage duration (exponential rv)
%       third row: blockage departure time

numBlockers = params.numBlockers;
width_area  = params.areaDimensions(1);
length_area  = params.areaDimensions(2);
% width_area  = params.areaDimensions_sub6(1);
% length_area  = params.areaDimensions_sub6(2);
numUE = params.numUE;
numUE_sub6 = params.numUE_sub6;
if (ue_idx<=numUE)
    UE_location = params.UE_locations(ue_idx,:);
else
    UE_location = params.UE_locations_sub6(ue_idx-numUE,:);
end
hr = params.hr;
ht = params.ht;
% locationsBS = [params.locationsBS; params.locationsBS_sub6];
% locationsBS = params.locationsBS_sub6;
locationsBS = params.locationsBS;
V=params.V;
hb = params.hb;
simTime = params.simTime;
mu = params.mu;


%% Generating Blocker Mobility
if isfile('all_mobility_precompute.mat')
    % File exists.
    load('all_mobility_precompute.mat','all_mobility')
else
     % File does not exist.
    mobility_input = cell(1,2);
    all_mobility = cell(1,2);
    for indB=1:length(numBlockers)
        mobility_input{indB} = struct('V_POSITION_X_INTERVAL',[-width_area/2 width_area/2],...%(m)
            'V_POSITION_Y_INTERVAL',[-length_area/2 length_area/2],...%(m)
            'V_SPEED_INTERVAL',[V V],...%(m/s)
            'V_PAUSE_INTERVAL',[0 0],...%pause time (s)
            'V_WALK_INTERVAL',[1.00 60.00],...%walk time (s)
            'V_DIRECTION_INTERVAL',[-180 180],...%(degrees)
            'SIMULATION_TIME',simTime,...%(s)
            'NB_NODES',numBlockers(indB));

        % Generate_Mobility function is Copyright (c) 2011, Mathieu Boutin
        all_mobility{indB} = Generate_Mobility(mobility_input{indB});
    end

end

%% Computing Physical interaction btw Blockers and Base Stations
% Dependent Variables

nBS = size(locationsBS,1);
% Lets create nBS by 4 matrix representing line segments for UE-BS blockage
% area:
%   Each row is of the form [x1 y1 x2 y2] where (x1,y1) is the start point and
%   (x2,y2) is the end point of a line segment:
%
%                  Line Segment
%       o--------------------------------o
%       ^                                ^
%    (x1,y1)                          (x2,y2)

%BS blocker interactions - Physical Characteristics of the system
%Consider the triangle from tip of the BS to the UE antenna
% BS---------------
% --NB--------BL--
%----------------UE




for indBC=1:length(numBlockers) %For every blocker count seperate systems
    nB = numBlockers(indBC);
    s_mobility = all_mobility{indBC};
    
    dataBS = cell(nBS,1);
    for indBS = 1:nBS
        dataBS{indBS}=cell(1,nB);
    end
    for indB = 1:nB 
        blocker_height = hb(indB);
        frac = (blocker_height-hr)/(ht-hr);
        BS_blockage_coordinates = UE_location + frac*(locationsBS-UE_location);
        UE_BS_blockage_line_segments = [repmat(UE_location,[nBS,1]), BS_blockage_coordinates];
        blockerMovementStartTimes = s_mobility.VS_NODE(indB).V_TIME(1:end-1);
        blockerMovementVelocities = sqrt((s_mobility.VS_NODE(indB).V_SPEED_X(1:end-1)).^2+ ...
            (s_mobility.VS_NODE(indB).V_SPEED_Y(1:end-1)).^2);
        blockerMovementSegments= [s_mobility.VS_NODE(indB).V_POSITION_X(1:end-1),...
            s_mobility.VS_NODE(indB).V_POSITION_Y(1:end-1),...
            s_mobility.VS_NODE(indB).V_POSITION_X(2:end),...
            s_mobility.VS_NODE(indB).V_POSITION_Y(2:end)];
        intersectionResults = lineSegmentIntersect(blockerMovementSegments,UE_BS_blockage_line_segments);
        for indBS = 1:nBS
            blockages = find(intersectionResults.intAdjacencyMatrix(:,indBS));
            intersectionX = intersectionResults.intMatrixX(blockages,indBS);
            intersectionY = intersectionResults.intMatrixY(blockages,indBS);
            vectorDistances = [intersectionX,intersectionY] - blockerMovementSegments(blockages,1:2);
            scalDistances = sqrt(sum(vectorDistances.^2,2));
            intersectionTime = blockerMovementStartTimes(blockages) + scalDistances./ blockerMovementVelocities(blockages);
            dataBS{indBS}{indB} = sort(intersectionTime(intersectionTime<simTime))';
        end
    end
    
    for indBS = 1:nBS
        dataBS{indBS}=cell2mat(dataBS{indBS});
    end
    
    
    for indBS = 1:nBS
        % len =length(dataBS{(ue_idx-1)*nBS+indBS});
        len =length(dataBS{indBS});
        dataBS{indBS}(2,:) =  exprnd(1/mu,1,len); % block duration
        dataBS{indBS}(3,:) = dataBS{indBS}(2,:) + dataBS{indBS}(1,:); % end of physical blockages\
        %if a blocker arrives before the previous blocker served then that is a
        %one long blockage, for programming purposes we delete the second
        %arrival and make one long combined blockage
        [~,sort_idx]=sort(dataBS{indBS}(1,:));
        dataBS{indBS} = dataBS{indBS}(:,sort_idx);
        isCleared = 0;
        while ~isCleared
            isCleared=1;
            len = size(dataBS{indBS},2);
            for jj=len:-1:2
                if dataBS{indBS}(3,jj-1) >= dataBS{indBS}(1,jj)
                    isCleared=0;
                    dataBS{indBS}(3,jj-1) = max(dataBS{indBS}(3,jj),dataBS{indBS}(3,jj-1));
                    dataBS{indBS}(:,jj) = [];
                    dataBS{indBS}(2,jj-1) = dataBS{indBS}(3,jj-1) - dataBS{indBS}(1,jj-1);
                end
            end  
        end
        dataBS{indBS}= dataBS{indBS}(:,dataBS{indBS}(3,:)<simTime);
    end
    
end
end

