%% this script is for checking number of L/R responses for each subject...
%% and for counting the number of missed trials
clear all;
data_dir_str= '/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/beha_behavior';

cd(data_dir_str)
files = dir('trust*.mat');
num_of_subjects = length(files);

responses = double(zeros(num_of_subjects,4));
%index=19
for index=1:num_of_subjects
    filename = files(index).name;
    fprintf('File processing: %s\n', filename);
    load(filename);
    filename = files(index).name;
    a = unique(b.partnerchoice_RESP, 'stable');
    if ~iscell(a)
        a = num2cell(a);
    end
    nonempty_cells = find(~cellfun(@isempty,a));
    [rows,columns] = size(a);
    a = a(nonempty_cells);
    if rows==3
        a(1,:) = [];
    end
    c = cellfun(@(x) sum(ismember(b.partnerchoice_RESP, x)),a, 'un',0);
    c = cell2mat(c);
    missed_trials = sum(b.decisions==0);       
    c = c'; %number of p and q responses
    id = filename(isstrprop(filename,'digit')); 
    responses(index, :) = [str2num(id) c missed_trials];   
end
ID = responses(:,1);
P_2 = responses(:,2);
Q_7 = responses(:,3);
Missed = responses(:,4);
coarse_beha = table(ID,P_2,Q_7,Missed);%$p = 2, q = 7
coarse_beha.diff = abs(coarse_beha.P_2 - coarse_beha.Q_7);
data_dir_str= '/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/beha_behavior/group_data/scan_responses';
save(data_dir_str,'coarse_beha');

clear all;

data_dir_str= '/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior';

cd(data_dir_str)
subject_files = dir('trust*.mat');
num_of_subjects = length(subject_files);

responses = double(zeros(num_of_subjects,4));
for ct=1:num_of_subjects
    file_name = subject_files(ct).name;
    fprintf('File processing: %s\n', file_name);
    load(file_name);
    id = str2num(file_name(isstrprop(file_name,'digit')));
    if iscell(b.partnerchoice_RESP)
        a = unique(b.partnerchoice_RESP, 'stable');
        nonempty_cells = find(~cellfun(@isempty,a));        
        a = a(nonempty_cells);
        c = cellfun(@(x) sum(ismember(b.partnerchoice_RESP, x)),a, 'un',0);
        c = cell2mat(c);
        c = c'; %number of p and q responses
        %c = [c missed_trials];
        missed_trials = sum(b.decisions==0);
    else
        c = [sum(b.partnerchoice_RESP==2) sum(b.partnerchoice_RESP==7)];
        missed_trials = sum(b.partnerchoice_RESP == -999);
    end    
    responses(ct, :) = [id c missed_trials];
    clear c;
end
ID = responses(:,1);
P_2 = responses(:,2);
Q_7 = responses(:,3);
Missed = responses(:,4);
coarse_behavior = table(ID,P_2,Q_7,Missed);
coarse_behavior.diff = abs(coarse_behavior.P_2 - coarse_behavior.Q_7);

data_dir_str= '/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/group_data/scan_responses';
save(data_dir_str,'coarse_behavior');