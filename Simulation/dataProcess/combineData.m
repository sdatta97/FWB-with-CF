%Taken from 
% https://in.mathworks.com/matlabcentral/answers/
%   538119-how-to-import-to-matlab-many-text-files-as-table-type
projectdir = 'C:\Users\dsoum\Desktop\data\outageData\outageResultsvaryinglbthres_mmse_updated\lambda_0.01';
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
% changingVars = cell(1,length(colNames)-2);
% for i=1:(length(colNames)-2)
%     changingVars{i} = colNames{i};
% end
% 
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
summaryTable  = groupsummary(combinedTable,changingVars,{'mean','std'});

writetable(summaryTable,'./outage_mmse_updated_vary_lb_thres_low_bl.txt')
writetable(summaryTable,'./outage_mmse_updated_vary_lb_thres_low_bl.csv')