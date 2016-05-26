clear variables;
close all;
masterlist = dlmread('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/HC list.txt');
regslocation = '/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/surveys/';

cd(regslocation)
files = dir('*survey*.mat');
all_ids = arrayfun(@(x) (x.name(isstrprop(x.name,'digit'))), files, 'UniformOutput', false);
all_ids = arrayfun(@(x) str2double(x), all_ids, 'UniformOutput', false);
all_ids = cell2mat(all_ids);

regressors_exist = intersect(masterlist, all_ids);
regressors_NTM = setxor(masterlist, regressors_exist);

 