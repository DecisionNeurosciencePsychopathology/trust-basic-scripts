data_dir_str= '/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior'; 
cd(data_dir_str)
files1 = dir('*.mat');

data_dir_str= '/Users/polinavanyukov/Box Sync/Project Trust Game/data/scan unprocessed/'; 
cd(data_dir_str)
files2 = dir('*');

ids1 = char(files1.name);
ids2 = char(files2.name);
ids1 = cellstr(ids1);
ids2 = cellstr(ids2);
ids1 = cellfun(@(x) x(isstrprop(x, 'digit')), ids1, 'UniformOutput', false);
ids2 = cellfun(@(x) x(isstrprop(x, 'digit')), ids2, 'UniformOutput', false);
ids1 = cellfun(@(x) str2num(x), ids1, 'UniformOutput', false);
ids2 = cellfun(@(x) str2num(x), ids2, 'UniformOutput', false);
ids1 = cell2mat(ids1);
ids2 = cell2mat(ids2);

unprocessed = setdiff(ids2, ids1);
unlisted = setdiff(ids1, ids2);
%errors = [110215 115438 210100 219863 220130];%behavioral subjects that still need to be processed, not included in the report
%
to_process=unprocessed([5 6]);
for c=1:length(unprocessed)
    trustbehavior(unprocessed(c));
end

data_dir_str = '/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior';
files = dir(data_dir_str);
cd(data_dir_str);
files=dir('trust*.mat');
ids = char(files.name);
ids = cellstr(ids);
ids = cellfun(@(x) x(isstrprop(x, 'digit')), ids, 'UniformOutput', false);
ids = cellfun(@(x) str2num(x), ids, 'UniformOutput', false);
ids = cell2mat(ids);

data_dir_str = '/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/beha_behavior';
files = dir(data_dir_str);
cd(data_dir_str);
files=dir('*.mat');
ids = char(files.name);
ids = cellstr(ids);
ids = cellfun(@(x) x(isstrprop(x, 'digit')), ids, 'UniformOutput', false);
ids = cellfun(@(x) str2num(x), ids, 'UniformOutput', false);
ids = cell2mat(ids);

