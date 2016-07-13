clear all;
%% Grabbing files
data_dir_str= '/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior';
data_dump_str = '/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior';

if ~exist(data_dump_str,'file')
    mkdir(data_dump_str)
    fprintf('Creating id specific reg folder in: %s\n\n',data_dump_str);
end

cd(data_dir_str)
files = dir('trust*.mat');
num_of_subjects = length(files);

conditions = double(zeros(num_of_subjects,5));

%for index = 7
for    index=1:num_of_subjects
    filename = files(index).name;
    fprintf('File processing: %s\n', filename);
    string_ID = filename(isstrprop(filename,'digit'));
    if not(str2double(string_ID)==219956||str2double(string_ID)==220017||str2double(string_ID)==208572)
        load(filename);
        %% 1 = good; 2 = bad; 3 = neutral; 4 = computer
        trustee = ones(length(b.identity),1);
        trustee(strcmp(b.identity,'computer')) = 4; 
        trustee(strcmp(b.identity,'neutral')) = 3; 
        trustee(strcmp(b.identity,'bad')) = 2; 
        block1 = trustee(1);
        block2 = trustee(49);
        block3 = trustee(97);
        block4 = trustee(145);
        responses(index,:)=[str2double(string_ID) block1 block2 block3 block4];
    end
end

ID = responses(:,1);
Block1 = responses(:,2);
Block2 = responses(:,3);
Block3 = responses(:,4);
Block4 = responses(:,5);

OrderTrustee = table(ID, Block1, Block2, Block3, Block4);
save(data_dump_str, 'OrderTrustee');