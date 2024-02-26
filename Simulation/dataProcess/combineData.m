%Taken from 
% https://in.mathworks.com/matlabcentral/answers/
%   538119-how-to-import-to-matlab-many-text-files-as-table-type
projectdir = 'C:/Users/dsoum/Desktop/outageData/outageResults0.1_Pow_fac_1000_rmin_4_high_thres';
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
% changingVars = cell(1,length(colNames)-6);
% for i=1:(length(colNames)-6)
%     changingVars{i} = colNames{i};
% end
changingVars = cell(1,length(colNames)-4);
for i=1:(length(colNames)-4)
    changingVars{i} = colNames{i};
end
% changingVars = cell(1,length(colNames)-3);
% for i=1:(length(colNames)-3)
%     changingVars{i} = colNames{i};
% end
summaryTable  = groupsummary(combinedTable,changingVars,{'mean','std'});

writetable(summaryTable,'./outage_pfac_1000_rmin_4_high_thres.txt')
writetable(summaryTable,'./outage_pfac_1000_rmin_4_high_thres.csv')
