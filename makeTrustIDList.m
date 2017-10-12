%Quick and dirty way to grab all trust ids from Box Sync folder
%dir('C:\Users\wilsonj3\Box Sync\Suicide studies\data')
%data_dir=dir('C:\Users\wilsonj3\Box Sync\Suicide studies\data');
data_dir=dir('E:\Box Sync\Project Trust Game\data\processed\scan_behavior');
data_dir=glob('E:\Box Sync\Project Trust Game\data\processed\scan_behavior\trust*.mat');

len = length(data_dir);

trust_ids = [];
for  i = 1:length(data_dir)
    str = data_dir{i};
    expression = '[0-9]{4,6}';
    match_case = regexp(str,expression,'match');
    if ~isempty(match_case)
        match_case=str2double(cell2mat(match_case));
        trust_ids = [trust_ids match_case];
    end
end

trust_ids = trust_ids';

save trust_ids trust_ids