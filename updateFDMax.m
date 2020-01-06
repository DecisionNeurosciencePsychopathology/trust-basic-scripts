%To update the fdmax.mat file

%% to set location for the input/output files
os = computer;
if strcmp(os(1:end-2),'PCWIN')
    path_to_fd_output = 'C:\kod\trust_basic_scripts\';
else
    [~, me] = system('whoami');
    me = strtrim(me);    
    if strcmp(me,'polinavanyukov')==1
        path_to_fd_output = '/Users/polinavanyukov/Scripts/GitHub/trust_clinical/';
    else
        path_to_fd_output = glob('?');
    end
end

%% function code begins here

fd_mean = readtable(strcat(path_to_fd_output,'fd_mean_output_trust.csv'));
fd_max = readtable(strcat(path_to_fd_output,'fd_max_output_trust.csv'));

fd_mean = sortrows(fd_mean,'Subjects','ascend');
fd_max = sortrows(fd_max,'Subjects','ascend');

% load('C:\kod\Neuropsych_preproc\matlab\db\subjIDlistDB.mat')
% idNumbers = fd_max.Subjects;
% id_idx=ismember(subjectIDlistDB.id_number,idNumbers);


% row_name = ['Max' num2str(block)];
% fd_max.(row_name)(id_idx)

save fd_max fd_max