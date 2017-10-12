function trust_arc_proc
%Very superficial way to implement LEARN trust into ARC

%Load the sceptic config files and initialize the tracking data.
%As a new user you will have to create the config files (see
%https://github.com/DecisionNeurosciencePsychopathology/temporal_instrumental_agent
%for more help) & set the paths to said file.
task_data=initialize_task_tracking_data('trust');

%Trust has it;s own version of updating subjects for now...
processNewTGScanSubjs

%File I/o - grab all scan subjects' .mat files
data_dir_str= 'E:/Box Sync/Project Trust Game/data/processed/scan_behavior/';
files = dir([data_dir_str 'trust*.mat']);

%So trust is in a weird state right now, just to get this rolling create a
%id_list and use another for loop to move all the subjects regs
id_list=[];
currfolder='E:/data/trust/regs/asterisked_regs';
%Dir pointer to Throndike (server) reg housing
dest_folder='/Volumes/bek/trust_analyses/regs';

for i = 1:length(files)
    
    %Pull id number from mat file
    id = regexp(files(i).name,'\d{4,6}', 'match');
    id=str2double(id{:});
    id_list = [id_list; id];
    
    %Update task_tracking data if they have a .mat they were processed
    %successfully
    task_data.behave_completed=1;
    
    %processNewTGScanSubjs processes behav data
            
    %Update task_tracking data
    task_data.behave_processed=1;
    
    %Make and move regressors last
    if i==length(files)
        %TODO automate the fx_max process
        
        %Construct regressors
        trustmakeregressor_group()
        asterisk()
        
        %move regs          
        moveregs(currfolder,id_list(1),dest_folder)
    end
    
    %write the task data to file
    record_subj_to_file(id,task_data)
end