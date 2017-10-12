function makeBaselineForTrust(b,contrast,suffix)
%Code to make the baseline model for afni, normalize the trustee_BA
%data for each subject and use -stim_base to load it in the model
%(See michaels email for more info)


trial1_index = 1;
trial48_index = 48;
blocks = 4;

volumewise_baseline=[];
tmp=[];
for i = 1:blocks
    tmp = repmat(contrast(trial1_index),b.block_length,1)';

    volumewise_baseline = [volumewise_baseline tmp];
    trial1_index = trial1_index+48;
end

%Convolove it with the HRF
convolved_baseline = conv(1*volumewise_baseline,b.hemoir);

%Downsample?
%tmp_baseline = gsresample(convolved_baseline,10,1./b.scan_tr);

%Cut off the tail
convolved_baseline = convolved_baseline(1:length(volumewise_baseline));

%Normalize to a max of 1
convolved_baseline = (convolved_baseline - min(convolved_baseline)) ./ ( max(convolved_baseline) - min(convolved_baseline) );

% Debug
% if max(convolved_baseline)~=1
%     stop=0;
% end

%fname = [num2str(b.id) '_baselineBA.dat'];
fname = [num2str(b.id) suffix];
dump_str = 'E:\data\trust\regs\baselines\';
gdlmwrite([dump_str fname],convolved_baseline','\t')











% for block_n= 1:blocks
%     %epoch_window = 0:b.bin_size:taskness.event_end(48, block_n);
%     epoch_window = taskness.event_beg(1,block_n):b.bin_size:taskness.event_end(48, block_n);
%    epoch_window_2 = b.stim_OnsetTime(x(block_n,1)):b.bin_size:b.stim_OnsetTime(x(block_n,1))+b.scan_tr*b.block_length*1000;
%     
%     
%     %Get event times
%     event_beg=taskness.event_beg(trial1_index:trial48_index);
%     event_end=taskness.event_end(trial1_index:trial48_index);
%     tmp_reg.(['regressors' num2str(block_n)]).BA_trustee =...
%         createSimpleRegressor(event_beg,event_end, epoch_window, b.trustee_BA(trial1_index:trial48_index));
%     
%     
%     % HRF-convolve all the event regressors
%     hrfregs = fieldnames(tmp_reg.regressors1);
%     for n = 1:numel(hrfregs)
%         % b.hrfreg1.RT
%         tmp_reg.(['hrfreg' num2str(block_n)]).(hrfregs{n}) = ...
%             conv(1*tmp_reg.(['regressors' num2str(block_n)]).(hrfregs{n}),b.hemoir);
%         
%         % cut off the tail after convolution and downsample
%         tmp = gsresample(tmp_reg.(['hrfreg' num2str(block_n)]).(hrfregs{n}),10,1./b.scan_tr);
%         tmp_reg.(['hrfreg' num2str(block_n)]).(hrfregs{n}) = tmp(1:b.block_length);
%     end
%     
%     trial1_index = trial1_index+48;
%     trial48_index = trial48_index+48;
%     
% end
% % concatenate everything
% fnm = fieldnames(tmp_reg.regressors1)';
% 
% %Added switch case for subjects with irregular trials
% ct=1:length(fnm);
% switch num_blocks
%     case 1
%         for ct=1:length(fnm)
%             b.hrf_regs.(fnm{ct}) = [tmp_reg.hrfreg1.(fnm{ct})];
%         end
%     case 2
%         for ct=1:length(fnm)
%             b.hrf_regs.(fnm{ct}) = [tmp_reg.hrfreg1.(fnm{ct}) tmp_reg.hrfreg2.(fnm{ct})];
%         end
%     case 3
%         for ct=1:length(fnm)
%             b.hrf_regs.(fnm{ct}) = [tmp_reg.hrfreg1.(fnm{ct}) tmp_reg.hrfreg2.(fnm{ct}) tmp_reg.hrfreg3.(fnm{ct})];
%         end        
%     case 4
%         for ct=1:length(fnm)
%             b.hrf_regs.(fnm{ct}) = [tmp_reg.hrfreg1.(fnm{ct}) tmp_reg.hrfreg2.(fnm{ct}) tmp_reg.hrfreg3.(fnm{ct}) tmp_reg.hrfreg4.(fnm{ct})];
%         end
%     otherwise
%         disp('Error occured somewhere')
% end
% 
% 
% baselineBA_data = zscore(b.trustee_BA); %Make sure this is volume wise not trial wise!!!
% 
% fname = [num2str(id) '_baselineBA.dat'];
% gdlmwrite([dump_str fname],baselineBA_data,'\t')
% 
% 
% 
% % function foo = createSimpleRegressor(event_begin,event_end,epoch_window,conditional_trial)
% % 
% % if(~exist('conditional_trial','var') || isempty(conditional_trial))
% %     conditional_trial = ones(length(event_begin),1);
% % end
% % 
% % % create epoch windows for each trial
% % epoch = arrayfun(@(a,b) a:b,event_begin,event_end,'UniformOutput',false);
% % 
% % foo = zeros(size(epoch_window));
% % 
% % for n = 1:numel(epoch)
% %     if(conditional_trial(n))
% %         foo = logical(foo + histc(epoch{n},epoch_window));
% %     end
% % end
% % 
% % return



% % % 
% % % baseline_data = zscore(baseline_condition); %Make sure this is volume wise not trial wise!!!
% % % 
% % % 
% % % 
% % % %fname = [num2str(id) '_baselineBA.dat'];
% % % gdlmwrite([dump_str fname],baseline_data,'\t')