data_dir_str= '/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior';
data_dump_str = '/Users/polinavanyukov/Box Sync/Project Trust Game/regs/';
path_of_motion_table = '/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/group_data/fd_max.mat';
load(path_of_motion_table);

if ~exist(data_dump_str,'file')
    mkdir(data_dump_str)
    fprintf('Creating id specific reg folder in: %s\n\n',data_dump_str);
end

cd(data_dir_str)
files = dir('trust*.mat');
num_of_subjects = length(files);

for k=18
%for k=1:58
    filename = files(k).name;
    fprintf('File processing: %s\n', filename);
    id = filename(isstrprop(filename,'digit'));
    if not(str2double(id)==219956||str2double(id)==220017)
        load(filename);
        if ischar(id) 
            id = str2double(id);
        end
        %First get index of id or kick out
        id_idx=find(ismember(fd_max.Subjects,id), 1);
        
        if isempty(id_idx)
            fprintf('Subject is not on censor block list: %d\n' , id);
        end
        censor_limit = 5;
        bad_block = strcmp(b.identity, 'good');
        trials_bad = find(bad_block,1);
        if trials_bad < 49
            block = 1;
        elseif trials_bad < 49+48
            block = 2;
        elseif trials_bad < 49 + 2*48
            block = 3;
        else
            block = 4;
        end
        row_name = ['Max' num2str(block)];
        
        if fd_max.(row_name)(id_idx)>=censor_limit
             fprintf('Remove file %d\n', id);
        end
        
    end
end