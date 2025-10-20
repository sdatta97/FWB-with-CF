%Taken from 
% https://in.mathworks.com/matlabcentral/answers/
%   538119-how-to-import-to-matlab-many-text-files-as-table-type
% projectdir = '//Users/sdatta/FWB-with-CF/Simulation/resultData/ppResults';
projectdir = '/Users/sdatta/Desktop/data/multiUEsscatterBS_final/lambda_250_UE_TN/outageResults';
% projectdir = '/Users/sdatta/FWB-with-CF/Simulation/resultData/algocompResults';
dinfo = dir(fullfile(projectdir, '*.csv'));   %use appropriate extension
filenames = fullfile({dinfo.folder}, {dinfo.name});
nfiles = length(filenames);
tables = cell(nfiles,1);
for K = 1 : nfiles
    tables{K} = readtable(filenames{K});
end

combinedTable = vertcat(tables{:});

colNames = combinedTable.Properties.VariableNames;
%% group the same parameter iterations and get stats mean and std, can add other stuff
% changingVars = cell(1,length(colNames)-8);
% for i=1:(length(colNames)-8)
%     changingVars{i} = colNames{i};
% end
% changingVars = cell(1,length(colNames)-6);
% for i=1:(length(colNames)-6)
%     changingVars{i} = colNames{i};
% end
changingVars = cell(1,length(colNames)-4);
for i=1:(length(colNames)-4)
    changingVars{i} = colNames{i};
end
% changingVars = cell(1,length(colNames)-5);
% for i=1:(length(colNames)-5)
%     changingVars{i} = colNames{i};
% end
% changingVars = cell(1,length(colNames)-2);
% for i=1:(length(colNames)-2)
%     changingVars{i} = colNames{i};
% end
% changingVars = cell(1,length(colNames)-1);
% for i=1:(length(colNames)-1)
%     changingVars{i} = colNames{i};
% end

summaryTable  = groupsummary(combinedTable,changingVars,{'mean','std','median'});
% summaryTable  = groupsummary(combinedTable,changingVars,{'mean'});

writetable(summaryTable,'./outage_data_lambda_UE_250_lambda_BS_vary_modify_scatter_fin_UE_TN.txt')
writetable(summaryTable,'./outage_data_lambda_UE_250_lambda_BS_vary_modify_scatter_fin_UE_TN.csv')