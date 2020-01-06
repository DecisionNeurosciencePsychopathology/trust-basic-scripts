function s = censorMovement(id,s,block)
%Censor entire blocks if max movement is greater than 5 per block
%Load in the data

%% to set location for the input/output files
os = computer;
if strcmp(os(1:end-2),'PCWIN')
    path_of_motion_table = 'C:\kod\trust_basic_scripts\fd_max.mat';
else
    [~, me] = system('whoami');
    me = strtrim(me);    
    if strcmp(me,'polinavanyukov')==1
        path_of_motion_table = '/Users/polinavanyukov/Scripts/GitHub/trust_clinical/fd_max.mat';
    else
        path_of_motion_table = glob('?');
    end
end

%% function code begins here
load(path_of_motion_table);
%load('fd_max.mat') %Change path as needed

if ischar(id)
    id = str2double(id);
end

%First get index of id or kick out
id_idx=find(ismember(fd_max.Subjects,id), 1);
if isempty(id_idx)
    warning(sprintf('Subject %d is not on censor block list\n', id))
    return
end

censor_limit = 5;
row_name = ['Max' num2str(block)];

%Censor the entire block if greater than the censor limit
if fd_max.(row_name)(id_idx)>=censor_limit
    s.(['regressors' num2str(block)]).to_censor = ...
        logical(s.(['regressors' num2str(block)]).to_censor .* 0);
    
end

