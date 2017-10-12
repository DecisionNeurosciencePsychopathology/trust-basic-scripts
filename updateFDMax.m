%To update the fdmax.mat file

fd_mean = readtable('C:\kod\trust_basic_scripts\fd_trust_mean_output.csv');
fd_max = readtable('C:\kod\trust_basic_scripts\fd_trust_max_output.csv');

fd_mean = sortrows(fd_mean,'Subjects','ascend');
fd_max = sortrows(fd_max,'Subjects','ascend');

% load('C:\kod\Neuropsych_preproc\matlab\db\subjIDlistDB.mat')
% idNumbers = fd_max.Subjects;
% id_idx=ismember(subjectIDlistDB.id_number,idNumbers);


% row_name = ['Max' num2str(block)];
% fd_max.(row_name)(id_idx)

save fd_max fd_max