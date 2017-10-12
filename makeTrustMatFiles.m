function proc_ids=makeTrustMatFiles(trust_ids,destin_path,local_dir)
%This script is an intermediary that will attempt to make .mat files for
%any new subjects. We used to run makeTrustIDList, but this pulls from box,
%which is where all the subjects data is eventually stored, but not where
%new subjects are stored. Therefore I suggest makeTrustIDList only be used
%if trust_ids becomes corrupted or deleted.

%NOTE: for behaviour subjects only put the raw data files into folders and
%run this command with trust_ids= new subject id's. There shouldn't be a
%lot of them.

global skip_subjs

%load('trust_ids.mat')
%load('hallquist_trust_ids.mat')

current_dir = pwd; %Set a way back point since we cd into different dirs on trustbehavior --This is bad, we should not do this fix in future--

proc_ids = trust_ids; %assume everyone will be processed

for i = 1:length(trust_ids)
    try
        b = trustbehavior(trust_ids(i),destin_path,local_dir); %Make sure the paths are correct
    catch
        continue
        proc_ids(i)=[];
    end
end

cd(current_dir)

%If trustbehavior crashed for someone report them
if ~isempty(skip_subjs) 
    fprintf('Skipped subjects!:\n\n')
    dsip(skip_subjs)
end
    