%Taken from 
% https://in.mathworks.com/matlabcentral/answers/
%   538119-how-to-import-to-matlab-many-text-files-as-table-type
projectdir = '/Users/sdatta/FWB-with-CF/Simulation/resultData/impactResults';
% projectdir = '/Users/sdatta/Desktop/data/analysisResults';
dinfo = dir(fullfile(projectdir, 'handoff_impact_incl_small_10UE_25lambdaBS_250lambdaUE_100deployRange_400Blockers_randomHeight_*Min_rate50000000Pow_fac10lb_thres10.csv'));   %use appropriate extension
filenames = fullfile({dinfo.folder}, {dinfo.name});
nfiles = length(filenames);
tables = cell(nfiles,1);
for K = 1 : nfiles
    tables{K} = readtable(filenames{K});
end

combinedTable = vertcat(tables{:});

colNames = combinedTable.Properties.VariableNames;
%% group the same parameter iterations and get stats mean and std, can add other stuff
% changingVars = cell(1,length(colNames)-2);
% for i=1:(length(colNames)-2)
%     changingVars{i} = colNames{i};
% end
% 
changingVars = cell(1,length(colNames)-6);
for i=1:(length(colNames)-6)
    changingVars{i} = colNames{i};
end
% changingVars = cell(1,length(colNames)-4);
% for i=1:(length(colNames)-4)
%     changingVars{i} = colNames{i};
% end
% changingVars = cell(1,length(colNames)-6);
% for i=1:(length(colNames)-6)
%     changingVars{i} = colNames{i};
% end
% changingVars = cell(1,length(colNames)-2);
% for i=1:(length(colNames)-2)
%     changingVars{i} = colNames{i};
% end

summaryTable  = groupsummary(combinedTable,changingVars,{'mean','std','median'});

writetable(summaryTable,'./impact_single_antenna_multi_UEv2_comp_cf_small_50Mbps.txt')
writetable(summaryTable,'./impact_single_antenna_multi_UEv2_comp_cf_small_50Mbps.csv')