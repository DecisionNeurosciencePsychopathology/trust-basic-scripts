%%%This script will check if any new data exists on skinner (Google
%%%Drive) and process the ePrime .txt to .mat files. The .mat files should
%%%then be available on box? Why are we still using box again???

%TODO 
%     -- encorporate this into the automatic pipeline
%     -- automate the fd_max.mat generation

%Start log
t=datetime('now','TimeZone','local','Format','yyyy-MM-dd');
diaryfile = ['log_' sprintf('%s',t)];
diary(diaryfile)

%Load in the current trust game scan ids
load('trust_ids.mat') 

%Grab current scan subjs from Skinner 
path_to_skinner = 'E:\Box Sync\'; %(adjust path as needed)
trust_data_path = 'skinner\projects_analyses\Project Trust\data\eprime\'; %This should always be the same!
full_data_path = [path_to_skinner trust_data_path];

%List the data dirs, reduce to only subject dirs
skinner_ids=dir(full_data_path);
skinner_ids={skinner_ids.name}';
skinner_ids=regexp(skinner_ids,'[0-9]+','match'); %Filter ids
skinner_ids=cell2mat(cellfun(@str2num,[skinner_ids{:}],'UniformOutput',0))'; %Convert to matrix

%Are there any new ids
unprocessed_ids = skinner_ids(~ismember(skinner_ids,trust_ids));

%Does the local unporcessed dir exist
local_dir = 'C:/kod/trust_basic_scripts/unprocessed_data/';
if ~exist(local_dir,'file')
    mkdir(local_dir)
    mkdir([local_dir filesep 'archive'])
end

%If any subjects are missing data
bad_id_idx=[];

if ~isempty(unprocessed_ids)
    %Move contents of unprocessed ids to local dir
    for i=1:length(unprocessed_ids)
        unprocessed_id = unprocessed_ids(i);
        file_checker = []; %Does the proper file exist, we might need to be more specific here!
        remote_data=[full_data_path num2str(unprocessed_id)];
        file_checker = glob([remote_data filesep '*rs*.txt']);
        if ~isempty(file_checker)
            copyfile(remote_data, [local_dir filesep num2str(unprocessed_id)])
        else
            fprintf('%d''s folder contains no .txt file!!! \n\n',unprocessed_id)
            bad_id_idx=[bad_id_idx; i];
        end
    end
    unprocessed_ids(bad_id_idx)=[];
else
    fprintf('No new trust ids could be found! \n\n')
end

%Where to save the .mat files
destin_path = 'E:/Box Sync/Project Trust Game/data/processed/scan_behavior/';

%Make the mat files
processed_ids=makeTrustMatFiles(unprocessed_ids,destin_path,local_dir);

%Update the trust_ids matrix
trust_ids = sort([trust_ids; processed_ids]);
save('trust_ids','trust_ids');

%Log off
diary off 