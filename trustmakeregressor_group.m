function b = trustmakeregressor_group()
%Grabbing all subjects that have already been processed and turned 
%into .mat file. Then making regressors for each one.

%grabbing files
%data_dir_str= '/Volumes/bek/trust_analyses/regs/scan_behavior/';
%data_dump_str = '/Volumes/bek/trust_analyses/regs/';
%grabbing files
data_dir_str= '/Users/polinavanyukov/Box Sync/Project Trust Game/data/trust data recent/scan_behavior/';
data_dump_str = '/Users/polinavanyukov/Box Sync/Project Trust Game/regs/';

if ~exist(data_dump_str,'file')
    mkdir(data_dump_str)
    fprintf('Creating id specific reg folder in: %s\n\n',data_dump_str);
end

cd(data_dir_str)
files = dir('*.mat');
num_of_subjects = length(files);

%hard coded times
dur_choice_display = 300;
dur_feedback = 1200; %1100msec in the file?

%scanning parameters
scan_tr = 1.67;
block_length = 155;
block_end = block_length*scan_tr*1000; %in msec
hemoir = spm_hrf(scan_tr, [6,16,1,1,6,0,32]); % better than resampling and smoothing 
frequency_scale_hz = 10;
% this scale is in msec, but it is separated into bins of X
% Hz (defined by 'frequency_scale' above. the resulting
% output will be in the scale of X Hz.
bin_size = 1/frequency_scale_hz*1000; % convert Hz to ms

    

for index=1:num_of_subjects
    filename = files(index).name;
    fprintf('File processing: %s\n', filename);
    load(filename);

%     if str2num(filename(89:93)) == 46069
%     data_dump_str = strcat('/Volumes/bek/trust_analyses/regs/', filename(89:93));
%     else
%     data_dump_str = strcat('/Volumes/bek/trust_analyses/regs/', filename(88:93));
%     end
    %data_dump_str = strcat('/Volumes/bek/trust_analyses/regs/', num2str(id));
    data_dump_str = strcat('/Users/polinavanyukov/Box Sync/Project Trust Game/regs/', num2str(id));
    b.regs = [];    
%Subjects may be using different trial versions of the task, 
%which only affects how the timings for the scanner should be calculated. 
%For some participants first fixation onset is the ITI fixation(1). 
%For others, first fixation onset is labeled as such, 2nd fixation onset 
%= ITI fixation (1), etc.
    
    if b.ITIfixation_OnsetTime(192) == -999 || not(iscell(b.firstFixation_OnsetTime))
        firstfix_Onset = b.firstFixation_OnsetTime(1);
        trial2_ITI = 1;
        trial48_ITI = 47;
    else firstfix_Onset = b.ITIfixation_OnsetTime(1);
        trial2_ITI = 2;
        trial48_ITI = 48;
    end
    
    %for those participants for whom the partnerchoice offset time was not
    %recorded by e-prime
    if iscell(b.partnerchoice_OffsetTime)
        partnerchoice_OffsetTime=b.displaychoice_OnsetTime-1;
    else
        partnerchoice_OffsetTime=b.partnerchoice_OffsetTime;
    end
    
    fixations = [];
    fixations.event_beg=zeros(48,4);
    fixations.event_end=zeros(48,4); 
    taskness.event_beg =zeros(48,4);
    taskness.event_end=zeros(48,4);
    decision.event_beg=zeros(48,4);
    decision.event_end=zeros(48,4);
    feedback.event_beg=zeros(48,4);
    feedback.event_end=zeros(48,4);
    b.missed_trials = (b.partnerchoice_RESP==-999);
    
    trial1_index = 1;
    trial48_index = 48;
    
    for block= 1:4
        %for fixation screens
        if block == 1 && (b.ITIfixation_OnsetTime(192) == -999 || not(isempty([b.firstFixation_OnsetTime])))
            fixations.event_beg(:,block) = [0; b.ITIfixation_OnsetTime(trial2_ITI:trial48_ITI)-firstfix_Onset];
        else
            %fixations.event_beg(:,block) = b.ITIfixation_OnsetTime(trial2_ITI-1:trial48_ITI)-firstfix_Onset+block_end*(block-1);
            fixations.event_beg(:,block) = b.ITIfixation_OnsetTime(trial2_ITI-1:trial48_ITI)-firstfix_Onset;
        end        
        %fixations.event_end(:,block) = b.partnerchoice_OnsetTime(trial1_index:trial48_index) - firstfix_Onset+block_end*(block-1);
        fixations.event_end(:,block) = b.partnerchoice_OnsetTime(trial1_index:trial48_index) - firstfix_Onset;
        
        %for trial onset to offset; Taskness
%         taskness.event_beg(:,block) = b.partnerchoice_OnsetTime(trial1_index:trial48_index)-firstfix_Onset+block_end*(block-1);
%         taskness.event_end(:,block) = b.outcome_OffsetTime(trial1_index:trial48_index)-firstfix_Onset+block_end*(block-1);
        
        taskness.event_beg(:,block) = b.partnerchoice_OnsetTime(trial1_index:trial48_index)-firstfix_Onset;
        taskness.event_end(:,block) = b.outcome_OffsetTime(trial1_index:trial48_index)-firstfix_Onset;
        
        %for decision onset to response (motor response)
%         decision.event_beg(:,block) = b.partnerchoice_OnsetTime(trial1_index:trial48_index)-firstfix_Onset+block_end*(block-1);
%         decision.event_end(:,block) = partnerchoice_OffsetTime(trial1_index:trial48_index)-firstfix_Onset+block_end*(block-1); 
        decision.event_beg(:,block) = b.partnerchoice_OnsetTime(trial1_index:trial48_index)-firstfix_Onset;
        decision.event_end(:,block) = partnerchoice_OffsetTime(trial1_index:trial48_index)-firstfix_Onset; 
        
        %for feedback onset to offset
%         feedback.event_beg(:,block) = b.outcome_OnsetTime(trial1_index:trial48_index)-firstfix_Onset+block_end*(block-1);
%         feedback.event_end(:,block) = b.outcome_OffsetTime(trial1_index:trial48_index)-firstfix_Onset+block_end*(block-1);
        feedback.event_beg(:,block) = b.outcome_OnsetTime(trial1_index:trial48_index)-firstfix_Onset;
        feedback.event_end(:,block) = b.outcome_OffsetTime(trial1_index:trial48_index)-firstfix_Onset;
        
        %epoch window + missed trials + to censor regressor
        epoch_window = 0:bin_size:decision.event_end(48, block);   
       
        if any(b.missed_trials(trial1_index:trial48_index))==0
           tmp_reg.(['regressors' num2str(block)]).to_censor = ones(size(createSimpleRegressor(decision.event_beg, decision.event_end, epoch_window))); 
        else
           tmp_reg.(['regressors' num2str(block)]).to_censor = createSimpleRegressor(decision.event_beg,decision.event_end, epoch_window, b.missed_trials(trial1_index:trial48_index)); 
           tmp_reg.(['regressors' num2str(block)]).to_censor = ones(size(tmp_reg.(['regressors' num2str(block)]).to_censor)) - tmp_reg.(['regressors' num2str(block)]).to_censor;
        end
        
        if block < 4
            firstfix_Onset = b.ITIfixation_OnsetTime(trial2_ITI-1+48);
        end
        trial2_ITI=trial2_ITI+48;
        trial48_ITI=trial48_ITI+48;
        trial1_index = trial1_index+48;
        trial48_index = trial48_index+48;
        
        tmp = gsresample( ...
         [zeros(50,1)' tmp_reg.(['regressors' num2str(block)]).to_censor(1:end-51)], ...
         10,1./scan_tr);
        tmp = floor(tmp);
        tmp = [tmp ones(1, (block_length-1)-length(tmp))];
        tmp = [tmp zeros(1,155-length(tmp))];
        tmp_reg.(['regressors' num2str(block)]).to_censor = tmp;      
    end
    fixations.event_beg=reshape(fixations.event_beg,[192,1]);
    fixations.event_end=reshape(fixations.event_end,[192,1]);
    taskness.event_beg =reshape(taskness.event_beg,[192,1]);
    taskness.event_end=reshape(taskness.event_end,[192,1]);
    decision.event_beg=reshape(decision.event_beg,[192,1]);
    decision.event_end=reshape(decision.event_end,[192,1]);
    feedback.event_beg=reshape(feedback.event_beg,[192,1]);
    feedback.event_end=reshape(feedback.event_end,[192,1]);
    
    %concatenating
    b.to_censor = [tmp_reg.regressors1.to_censor tmp_reg.regressors2.to_censor tmp_reg.regressors3.to_censor tmp_reg.regressors4.to_censor]; 
    b.to_censor = transpose(b.to_censor);
    %plot(b.to_censor);
    
    %exclude missed trials
    b.notmissed_trials = (b.partnerchoice_RESP~=-999);
      
    [b.stim_times.fix_fsl,b.stim_times.fix_spmg]=write3Ddeconv_startTimes(data_dump_str,fixations.event_beg,fixations.event_end,'fixation_Times',b.notmissed_trials,0);
    [b.stim_times.task_fsl,b.stim_times.task_spmg]=write3Ddeconv_startTimes(data_dump_str,taskness.event_beg,taskness.event_end,'task_Times',b.notmissed_trials,0);
    [b.stim_times.resp_fsl,b.stim_times.resp_spmg]=write3Ddeconv_startTimes(data_dump_str,decision.event_beg,decision.event_end,'decision_Times',b.notmissed_trials,0);
   
    %censor missed trials
    b.censor_missed = createSimpleRegressor(taskness.event_beg,taskness.event_end,epoch_window,b.missed_trials);
 
    % right = 2; left = 7;
    b.leftVSright = zeros(size(b.partnerchoice_RESP));
    b.leftVSright(b.partnerchoice_RESP==7)=1;
    b.leftVSright(b.partnerchoice_RESP==2)=-1; 
       
    % share/keep;
    b.shareVSkeep = zeros(size(b.partnerchoice_RESP));
    b.shareVSkeep(b.decisions==1 & b.partnerchoice_RESP ~= -999) = 1;
    b.shareVSkeep(b.decisions==-1 & b.partnerchoice_RESP ~= -999) = -1;
    
    [b.stim_times.left_fsl,b.stim_times.left_spmg]=write3Ddeconv_startTimes(data_dump_str,decision.event_beg,decision.event_end,'leftVSright',b.leftVSright,0);
    [b.stim_times.share_fsl,b.stim_times.share_spmg]=write3Ddeconv_startTimes(data_dump_str,decision.event_beg,decision.event_end,'p_shareVSkeep',b.shareVSkeep,0);
       
    share =~cellfun(@isempty,strfind(b.TrusteeDecides(1:192),'share'));
    keep =~cellfun(@isempty,strfind(b.TrusteeDecides(1:192),'keep'));
    b.t_share = zeros(192,1);
    b.t_keep = zeros(192,1);
    b.t_share(share) = 1;
    b.t_share(keep) = -1;
   
    [b.stim_times.feedback_fsl,b.stim_times.feedback_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'feedback_Times',b.notmissed_trials,0);
    [b.stim_times.t_shared_fsl,b.stim_times.t_shared_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'tshared_Times',b.t_share,0);
    
    gdlmwrite(strcat(data_dump_str, 'to_censor'),[b.to_censor],'\t');
end



return


