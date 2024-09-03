%Taken from 
% https://in.mathworks.com/matlabcentral/answers/
%   538119-how-to-import-to-matlab-many-text-files-as-table-type
% projectdir = '/Users/sdatta/FWB-with-CF/Simulation/resultData/impactResults';
projectdir = '/Users/sdatta/Desktop/data/multiUEsqgridBS_vary_l_lambda_UE_250/outageResults';
dinfo = dir(fullfile(projectdir, '*s100.csv'));   %use appropriate extension
filenames = fullfile({dinfo.folder}, {dinfo.name});
nfiles = length(filenames);
tables = cell(nfiles,1);
for K = 1 : nfiles
    tables{K} = readtable(filenames{K});
end

combinedTable = vertcat(tables{:});

colNames = combinedTable.Properties.VariableNames;
%% group the same parameter iterations and get stats mean and std, can add other stuff
% changingVars = cell(1,length(colNames)-6);
% for i=1:(length(colNames)-6)
%     changingVars{i} = colNames{i};
% end
changingVars = cell(1,length(colNames)-4);
for i=1:(length(colNames)-4)
    changingVars{i} = colNames{i};
end
% changingVars = cell(1,length(colNames)-2);
% for i=1:(length(colNames)-2)
%     changingVars{i} = colNames{i};
% end

summaryTable  = groupsummary(combinedTable,changingVars,{'mean','std','median'});

writetable(summaryTable,'./outage_multi_UE_sqgridBS_lambda_UE_250_l_100.txt')
writetable(summaryTable,'./outage_multi_UE_sqgridBS_lambda_UE_250_l_100.csv')