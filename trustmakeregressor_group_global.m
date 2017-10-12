function b = trustmakeregressor_group_global()
%% Grabbing all subjects that have already been processed and turned 
%into .mat file. Then making regressors for each one.


%% Grabbing files
data_dir_str= '/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior';
data_dump_str = '/Users/polinavanyukov/Box Sync/Project Trust Game/regs/July 18/';

if ~exist(data_dump_str,'file')
    mkdir(data_dump_str)
    fprintf('Creating id specific reg folder in: %s\n\n',data_dump_str);
end

cd(data_dir_str)
files = dir('trust*.mat');
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

%% Subject loop
%for index = 8
for index=46
%for index=1:num_of_subjects
    filename = files(index).name;
    fprintf('File processing: %s\n', filename);
    id = filename(isstrprop(filename,'digit'));
    if not(str2double(id)==219956||str2double(id)==220017||str2double(id)==208572)
        load(filename);
        data_dump_str = strcat('/Users/polinavanyukov/Box Sync/Project Trust Game/regs/', num2str(id));
        b.regs = [];    
        
    %Subjects may be using different trial versions of the task, 
    %which only affects how the timings for the scanner should be calculated. 
    %For some participants first fixation onset is the ITI fixation(1). 
    %For others, first fixation onset is labeled as such, 2nd fixation onset 
    %= ITI fixation (1), etc.
        if not(exist('decisions'))
            share =~cellfun(@isempty,strfind(b.PartDecides,'share'));
            keep =~cellfun(@isempty,strfind(b.PartDecides,'keep'));
            noresponse = ~cellfun(@isempty,strfind(b.PartDecides,'noresponse'));
            b.decisions = zeros(192, 1);
            b.decisions(share) = 1;
            b.decisions(keep) = -1;
            b.decisions(noresponse) = 0;
        end
        
        if length(b.ITIfixation_OnsetTime) < 192
            b.ITIfixation_OnsetTime(length(b.ITIfixation_OnsetTime)+1:192)=-999;
            b.ITIfixation_OffsetTime(145:192)=num2cell(-999);
        end
        
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
        if iscell(b.partnerchoice_RESP)
            b.missed_trials = b.decisions == 0;
        else       
            b.missed_trials = (b.partnerchoice_RESP==-999);
        end
        b.notmissed_trials = ~(b.missed_trials); 

        %plot(b.notmissed_trials);
        trial1_index = 1;
        trial48_index = 48;
        
        for block= 1:4
            %for fixation screens
            if block == 1 && (b.ITIfixation_OnsetTime(192) == -999 || not(isempty([b.firstFixation_OnsetTime])))
                fixations.event_beg(:,block) = [0; b.ITIfixation_OnsetTime(trial2_ITI:trial48_ITI)-firstfix_Onset];
            else
                %fixations.event_beg(:,block) = b.ITIfixation_OnsetTime(trial2_ITI-1:trial48_ITI)-firstfix_Onset+block_end*(block-1);
                fixations.event_beg(:,block) = [155*1670*(block-1); b.ITIfixation_OnsetTime(trial2_ITI:trial48_ITI)-firstfix_Onset+155*1670*(block-1)];
            end        
            %fixations.event_end(:,block) = b.partnerchoice_OnsetTime(trial1_index:trial48_index) - firstfix_Onset+block_end*(block-1);
            fixations.event_end(:,block) = b.partnerchoice_OnsetTime(trial1_index:trial48_index) - firstfix_Onset+155*1670*(block-1);

            %for trial onset to offset; Taskness
            taskness.event_beg(:,block) = b.partnerchoice_OnsetTime(trial1_index:trial48_index)-firstfix_Onset+155*1670*(block-1);
            taskness.event_end(:,block) = b.outcome_OffsetTime(trial1_index:trial48_index)-firstfix_Onset+155*1670*(block-1);

            %for decision onset to response (motor response)
            decision.event_beg(:,block) = b.partnerchoice_OnsetTime(trial1_index:trial48_index)-firstfix_Onset+155*1670*(block-1);
            decision.event_end(:,block) = partnerchoice_OffsetTime(trial1_index:trial48_index)-firstfix_Onset+155*1670*(block-1); 

            %for feedback onset to offset
            feedback.event_beg(:,block) = b.outcome_OnsetTime(trial1_index:trial48_index)-firstfix_Onset+155*1670*(block-1);
            feedback.event_end(:,block) = b.outcome_OffsetTime(trial1_index:trial48_index)-firstfix_Onset+155*1670*(block-1); 
           
            %% New for aligning value: RT aligned to feedback
%             decision.event_beg(:,block) = b.partnerchoice_RTTime(trial1_index:trial48_index)-firstfix_Onset;
%             decision.event_end(:,block) = partnerchoice_OffsetTime(trial1_index:trial48_index)-firstfix_Onset; 

            %epoch window + missed trials + to censor regressor 
            epoch_window = 0:bin_size:taskness.event_end(48, block);
            %% There is a strange glitch, likely owing to the createSimpleRegressor code;
            %% modified by AD to include ITI
            if(sum(b.missed_trials(trial1_index:trial48_index))) > 0
                % write 1s for missed trials to censor
                tmp_reg.(['regressors' num2str(block)]).to_censor = createSimpleRegressor(taskness.event_beg,taskness.event_end, epoch_window, b.missed_trials(trial1_index:trial48_index));
            else
                % write a vector of 0s the size of regressors
                tmp_reg.(['regressors' num2str(block)]).to_censor = zeros(size(createSimpleRegressor(taskness.event_beg,taskness.event_end, epoch_window, b.notmissed_trials(trial1_index:trial48_index))));
            end
            % flip to_censor
            tmp_reg.(['regressors' num2str(block)]).to_censor = 1-tmp_reg.(['regressors' num2str(block)]).to_censor;

            %% GSR resample
            tmp = gsresample( ...
             [zeros(50,1)' tmp_reg.(['regressors' num2str(block)]).to_censor(1:end-51)], ...
             10,1./scan_tr);
            tmp = ceil(tmp);
            tmp = [tmp ones(1, (block_length-1)-length(tmp))];
            tmp = [tmp zeros(1,155-length(tmp))];
            tmp_reg.(['regressors' num2str(block)]).to_censor = tmp;
            %plot(tmp_reg.(['regressors' num2str(block)]).to_censor);
            
            %% Censoring blocks w/ movement
            tmp_reg=censorMovement(id, tmp_reg, block); 
            tmp_reg=censorComputer(id, tmp_reg, block);
            %plot(tmp_reg.(['regressors' num2str(block)]).to_censor);
            
            %for next loop iteration, reinitilize variables
            if block < 4
                firstfix_Onset = b.ITIfixation_OnsetTime(trial2_ITI-1+48);
            end
            trial2_ITI=trial2_ITI+48;
            trial48_ITI=trial48_ITI+48;
            trial1_index = trial1_index+48;
            trial48_index = trial48_index+48;       
        end
        fixations.event_beg=reshape(fixations.event_beg,[192,1]);
        fixations.event_end=reshape(fixations.event_end,[192,1]);
        taskness.event_beg =reshape(taskness.event_beg,[192,1]);
        taskness.event_end=reshape(taskness.event_end,[192,1]);
        decision.event_beg=reshape(decision.event_beg,[192,1]);
        decision.event_end=reshape(decision.event_end,[192,1]);
        feedback.event_beg=reshape(feedback.event_beg,[192,1]);
        feedback.event_end=reshape(feedback.event_end,[192,1]);
        %plotting to censor regressors
%         plot(tmp_reg.regressors3.to_censor);
%         figure(2);
%         plot(b.notmissed_trials(192-48:192));
        %concatenating
        b.to_censor = [tmp_reg.regressors1.to_censor tmp_reg.regressors2.to_censor tmp_reg.regressors3.to_censor tmp_reg.regressors4.to_censor]; 
 %       plot(b.to_censor);
        b.to_censor = transpose(b.to_censor);
        
      [b.stim_times.fix_fsl,b.stim_times.fix_spmg]=write3Ddeconv_startTimes(data_dump_str,fixations.event_beg,fixations.event_end,'gt_fixation_Times',b.notmissed_trials',0);
      [b.stim_times.task_fsl,b.stim_times.task_spmg]=write3Ddeconv_startTimes(data_dump_str,taskness.event_beg,taskness.event_end,'gt_task_Times',b.notmissed_trials',0);
      [b.stim_times.resp_fsl,b.stim_times.resp_spmg]=write3Ddeconv_startTimes(data_dump_str,decision.event_beg,decision.event_end,'gt_decision_Times',b.notmissed_trials',0);
   
        % right = 2; left = 7;
        b.leftVSright = zeros(size(b.partnerchoice_RESP));
        b.leftVSright(b.partnerchoice_RESP==7)=1;
        b.leftVSright(b.partnerchoice_RESP==2)=-1; 

        % subject's decisions, share/keep;
        b.shareVSkeep = zeros(size(b.partnerchoice_RESP));
        b.shareVSkeep(b.decisions==1 & b.partnerchoice_RESP ~= -999) = 1;
        b.shareVSkeep(b.decisions==-1 & b.partnerchoice_RESP ~= -999) = -1;

       [b.stim_times.left_fsl,b.stim_times.left_spmg]=write3Ddeconv_startTimes(data_dump_str,decision.event_beg,decision.event_end,'gt_leftVSright',b.leftVSright',0);
%       [b.stim_times.share_fsl,b.stim_times.share_spmg]=write3Ddeconv_startTimes(data_dump_str,decision.event_beg,decision.event_end,'p_shareVSkeep',b.shareVSkeep,0);
%       [b.stim_times.share_fsl,b.stim_times.share_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'gt_p_shareVSkeep',b.shareVSkeep,0);
      
        share =~cellfun(@isempty,strfind(b.TrusteeDecides(1:192),'share'));
        keep =~cellfun(@isempty,strfind(b.TrusteeDecides(1:192),'keep'));
        b.t_share = zeros(192,1);
        b.t_keep = zeros(192,1);
        b.t_share(share) = 1;
        b.t_share(keep) = -1;

        [b.stim_times.feedback_fsl,b.stim_times.feedback_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'gt_feedback_Times',b.notmissed_trials',0);
%        [b.stim_times.t_shared_fsl,b.stim_times.t_shared_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'tshared_Times',b.t_share,0);

       %% congruent vs. incongruent trials
       congruent10 = b.shareVSkeep == b.t_share;
       %[b.stim_times.congruent_fsl,b.stim_times.congruent_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'congruent10',congruent10,0);
       
        %% trustee contrasts 
        b.trustee_GB = zeros(length(b.identity),1);
        b.trustee_G0 = zeros(length(b.identity),1);
        b.trustee_HC = ones(length(b.identity),1);
%        b.trustee_BA = ones(length(b.identity),1)*-1;
        b.trustee_BA = zeros(length(b.identity),1);
        b.trustee_survey_GB = zeros(length(b.identity),1);
        b.trustee_action = ones(length(b.identity),1);
        
        % good VS bad
        b.trustee_GB(strcmp(b.identity,'good')) = 1;         
        b.trustee_G0(strcmp(b.identity,'good')) = 1;
        b.trustee_GB(strcmp(b.identity,'bad')) = -1;
        
%         [b.stim_times.trusteeGB_fsl,b.stim_times.trusteeGB_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'trusteeGB',b.trustee_GB,0);        
%         [b.stim_times.trusteeG0_fsl,b.stim_times.trusteeG0_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'trusteeG0',b.trustee_G0,0);        
%         [b.stim_times.trusteeGB_fsl,b.stim_times.trusteeGB_spmg]=write3Ddeconv_startTimes(data_dump_str,decision.event_beg,decision.event_end,'trusteeGB',b.trustee_GB,0);        
%         [b.stim_times.trusteeG0_fsl,b.stim_times.trusteeG0_spmg]=write3Ddeconv_startTimes(data_dump_str,decision.event_beg,decision.event_end,'trusteeG0',b.trustee_G0,0);        

        % human vs non-human
        b.trustee_HC(strcmp(b.identity,'computer')) = 0;
%       [b.stim_times.trusteeHC_fsl,b.stim_times.trusteeHC_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'trusteeHC',b.trustee_HC,0);
%       [b.stim_times.trusteeHC_fsl,b.stim_times.trusteeHC_spmg]=write3Ddeconv_startTimes(data_dump_str,decision.event_beg,decision.event_end,'trusteeHC',b.trustee_HC,0);        
       
       % bad vs rest
       b.trustee_BA(strcmp(b.identity,'bad')) = 1;
%       [b.stim_times.trusteeBA_fsl,b.stim_times.trusteeBA_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'trusteeBA',b.trustee_BA,0); 
%       [b.stim_times.trusteeBA_fsl,b.stim_times.trusteeBA_spmg]=write3Ddeconv_startTimes(data_dump_str,decision.event_beg,decision.event_end,'trusteeBA',b.trustee_BA,0);
%       [b.stim_times.trusteeBA_fsl,b.stim_times.trusteeBA_spmg]=write3Ddeconv_startTimes(data_dump_str,decision.event_beg,feedback.event_end,'gt_trusteeBA',b.trustee_BA',0);
            
       % trustee action regressor
       b.trustee_action(strcmp(b.TrusteeDecides,'keep'))=-1;
%       [b.stim_times.trustee_action_fsl,b.stim_times.trustee_action_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'trusteeAction',b.trustee_action,0);
       
       
       % GB contrast needs a censor regressor that censors other trustee
       % blocks
       b.trusteeBYtactionGB = b.trustee_GB.*b.trustee_action;
%       [b.stim_times.trusteeBYactionGB_fsl,b.stim_times.trusteeBYactionGB_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'GBxtaction',b.trusteeBYtactionGB,0);
       b.trusteeBAxp_shareVSkeep = b.trustee_BA.*b.shareVSkeep;  
       b.trusteeGBxp_shareVSkeep = b.trustee_GB.*b.shareVSkeep;      
       b.trusteeG0xp_shareVSkeep = b.trustee_G0.*b.shareVSkeep;     
       b.trusteeHCxp_shareVSkeep = b.trustee_HC.*b.shareVSkeep;
       
%        [b.stim_times.trusteeBYactionG0_fsl,b.stim_times.trusteeBYactionG0_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'G0xp_decision',b.trusteeG0xp_shareVSkeep,0);
%        [b.stim_times.trusteeBYactionGB_fsl,b.stim_times.trusteeBYactionGB_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'GBxp_decision',b.trusteeGBxp_shareVSkeep,0);
%        [b.stim_times.trusteeBYactionBA_fsl,b.stim_times.trusteeBYactionBA_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'BAxp_decision',b.trusteeBAxp_shareVSkeep,0);
%        [b.stim_times.trusteeBYactionHC_fsl,b.stim_times.trusteeBYactionHC_spmg]=write3Ddeconv_startTimes(data_dump_str,feedback.event_beg,feedback.event_end,'HCxp_decision',b.trusteeHCxp_shareVSkeep,0);
        
%        [b.stim_times.trusteeBYactionG0_fsl,b.stim_times.trusteeBYactionG0_spmg]=write3Ddeconv_startTimes(data_dump_str,decision.event_beg,decision.event_end,'G0xp_decision',b.trusteeG0xp_shareVSkeep,0);
%        [b.stim_times.trusteeBYactionGB_fsl,b.stim_times.trusteeBYactionGB_spmg]=write3Ddeconv_startTimes(data_dump_str,decision.event_beg,decision.event_end,'GBxp_decision',b.trusteeGBxp_shareVSkeep,0);
%        [b.stim_times.trusteeBYactionBA_fsl,b.stim_times.trusteeBYactionBA_spmg]=write3Ddeconv_startTimes(data_dump_str,decision.event_beg,decision.event_end,'BAxp_decision',b.trusteeBAxp_shareVSkeep,0);
%        [b.stim_times.trusteeBYactionHC_fsl,b.stim_times.trusteeBYactionHC_spmg]=write3Ddeconv_startTimes(data_dump_str,decision.event_beg,decision.event_end,'HCxp_decision',b.trusteeHCxp_shareVSkeep,0);
 
%       b.surveyGBbytactionGB = b.trustee_survey_GB * b.trustee_action;
%         gdlmwrite(strcat(data_dump_str, 'to_censor'),[b.to_censor],'\t');
%         gdlmwrite(strcat(data_dump_str, 'to_censor_motion_corr'),[b.to_censor],'\t');          
%         gdlmwrite(strcat(data_dump_str, 'to_censor_PEs'),[b.to_censor_PEmaps],'\t');
         gdlmwrite(strcat(data_dump_str, 'to_censor_computer'),[b.to_censor],'\t');
%% block regressors
        block = [];
        block.event_beg = [0 fixations.event_beg(49) fixations.event_beg(97) fixations.event_beg(145)]';
        block.event_end = [feedback.event_end(48) feedback.event_end(96) feedback.event_end(144) feedback.event_end(192)]';
        block.BA = [b.trustee_BA(1) b.trustee_BA(49) b.trustee_BA(97) b.trustee_BA(145)]';
        block.G0 = [b.trustee_G0(1) b.trustee_G0(49) b.trustee_G0(97) b.trustee_G0(145)]';
        block.HC = [b.trustee_HC(1) b.trustee_HC(49) b.trustee_HC(97) b.trustee_HC(145)]';
         [b.stim_times.blockBA_fsl, b.stim_times.blockBA_spmg] = write3Ddeconv_startTimes(data_dump_str,block.event_beg,block.event_end,'gt_blockBA',block.BA',0);
%         [b.stim_times.blockG0_fsl, b.stim_times.blockG0_spmg] = write3Ddeconv_startTimes(data_dump_str,block.event_beg,block.event_end,'blockG0',block.G0,0);
%         [b.stim_times.blockHC_fsl, b.stim_times.blockHC_spmg] = write3Ddeconv_startTimes(data_dump_str,block.event_beg,block.event_end,'blockHC',block.HC,0);
        %% condition order regressor
        block.CO = [1 2 3 4]';
        block.BAxCO = block.BA.*block.CO;
        block.CO = block.CO - mean(block.CO);
         [b.stim_times.blockOrder_fsl, b.stim_times.blockOrder_spmg] = write3Ddeconv_startTimes(data_dump_str,block.event_beg,block.event_end,'gt_blockOrder',block.CO',0);
        block.BAxCO = block.BAxCO - mean(block.BA);
%        [b.stim_times.BAxOrder_fsl, b.stim_times.BAxOrder_spmg] = write3Ddeconv_startTimes(data_dump_str,block.event_beg,block.event_end,'BAxorder',block.BAxCO,0);
%% trial regressors
        b.trial = [1:192] - mean(1:192);
%       [b.stim_times.trial_fsl, b.stim_times.trial_spmg] = write3Ddeconv_startTimes(data_dump_str, fixations.event_beg, taskness.event_end, 'gt_trial1to192', b.trial,0);

%% rep consistent: BAD KEEP = 1; BAD SHARE = -1; GOOD SHARE = 1; GOOD KEEP = -1; NEUTRAL = 0
    bad_keeps1shareNEG1 = b.trustee_BA.* b.t_share*-1;
    good_shares1keepsNEG1 = b.trustee_G0.* b.t_share;
    consistent1elseNEG1 = bad_keeps1shareNEG1 + good_shares1keepsNEG1;
    consistent1elseNEG1 = consistent1elseNEG1 - mean(consistent1elseNEG1); 
    [b.stim_times.trial_fsl, b.stim_times.trial_spmg] = write3Ddeconv_startTimes(data_dump_str, feedback.event_beg, feedback.event_end, 'gt_consistWrep1elseNEG1',consistent1elseNEG1',0);

%% consistent w/ rep + sub_decision: BAD KEEP/SUBJECT KEPT = 1; BAD SHARE/SUBJECT SHARED = -1; MISMATCH IN BAD/GOOD CENSORED; NEUTRAL = 0;
    SgSs1KgKsNEG1 = zeros(length(good_shares1keepsNEG1),1);
    SgSs1KgKsNEG1(good_shares1keepsNEG1~=0) = good_shares1keepsNEG1(good_shares1keepsNEG1~=0)==b.shareVSkeep(good_shares1keepsNEG1~=0);
    SgSs1KgKsNEG1(good_shares1keepsNEG1~=0 & SgSs1KgKsNEG1 ==0) = -1;
    KbKs1SbSsNEG1 = zeros(length(bad_keeps1shareNEG1),1);
    b.shareNEG1VSkeep1 = b.shareVSkeep*-1;
    KbKs1SbSsNEG1(bad_keeps1shareNEG1~=0) = bad_keeps1shareNEG1(bad_keeps1shareNEG1~=0)==b.shareNEG1VSkeep1(bad_keeps1shareNEG1~=0);
    KbKs1SbSsNEG1(bad_keeps1shareNEG1~=0 & KbKs1SbSsNEG1 ==0) = -1;
    congruence1elseNEG1 = SgSs1KgKsNEG1+KbKs1SbSsNEG1;
    congruence1elseNEG1(b.trustee_BA==1|b.trustee_G0==1) = congruence1elseNEG1(b.trustee_BA==1|b.trustee_G0==1)  - mean(congruence1elseNEG1(b.trustee_BA==1|b.trustee_G0==1) );
%    [b.stim_times.trial_fsl, b.stim_times.trial_spmg] = write3Ddeconv_startTimes(data_dump_str, feedback.event_beg, feedback.event_end, 'gt_consistWrepNsd1elseNEG1',congruence1elseNEG1',0);
%% BAD KEEP/SUBJECT KEPT = 2; BAD SHARE/SUBJECT SHARED = -2; NEUTRAL SHARE/SUBJECT SHARED = 1; NEUTRAL = 0; MISMATCH IN BAD/GOOD/NEUTRAL CENSORED; 

%% PROSOCIAL CODING: TRUSTEE SHARES/SUBJECT KEPT = -1; TRUSTEE SHARES/SUBJECT SHARES = 0; TRUSTEE KEEPS/SUBJECT KEEPS = 0; TRUSTEE KEEPS/SUBJECT SHARES = 1
prosocial = zeros(length(b.trustee_action),1);
prosocial(b.trustee_action == 1 & b.shareVSkeep == -1) = -1;
prosocial(b.trustee_action == -1 & b.shareVSkeep == 1) = 1;
% plot(b.shareVSkeep(1:48));
% hold on;
% plot(b.trustee_action(1:48));
% figure(2)
% plot(prosocial(1:48));
prosocial(b.shareVSkeep~=0) = prosocial(b.shareVSkeep~=0) - mean(prosocial(b.shareVSkeep~=0));
[b.stim_times.trial_fsl, b.stim_times.trial_spmg] = write3Ddeconv_startTimes(data_dump_str, feedback.event_beg, feedback.event_end, 'gt_prosocial',prosocial',0);

    end
end



return


