clear all;

data_dir_str= '/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior';
cd(data_dir_str)
matfiles = dir('trust*.mat');
num_of_subjects = length(matfiles);
ct = 1;
errorarray = double.empty;
for index=1:num_of_subjects
    filename = matfiles(index).name;
    load(filename);
    start = 0;
    trials = length(b.TrialNumber);
    id = filename(isstrprop(filename,'digit'));
    
    id = str2double(id);
    if (trials<192 || trials>192)
        errorarray(ct) = double(id);
        fprintf('File processing: %s\n', filename);
        ct=ct+1;
    end
        
end

load('trust212856.mat');
subject_trials = length(b.ConditionOrder);
trials = 192;
to_fill = [subject_trials+1:trials];
b.ConditionOrder(to_fill) = 4;
b.identity(to_fill)={'-999'};
b.partnerchoice_RESP(to_fill) = -999;
b.TrusteeDecides(to_fill) = {'-999'};
b.TrialNumber(to_fill) = to_fill;
b.exchangeNum(to_fill)=[1:48];
b.partnerchoice_OnsetTime(to_fill) = -999;
b.partnerchoice_RTTime(to_fill) = -999;
b.partnerchoice_RT(to_fill) = -999;
b.partnerchoice_Duration(to_fill) = -999;
b.PartDecides(to_fill)={'-999'};
b.displaychoice_Duration(to_fill) = -999;
b.outcome_OnsetTime(to_fill)= -999;
b.outcome_OffsetTime(to_fill)=-999;
b.Order_RS(to_fill)=-999;
b.ITIfixation_OnsetTime(to_fill)=-999;
b.ITIfixation_Duration(to_fill)=-999;
b.ITIfixation_OffsetTime(to_fill)=-999;
b.ITIjitter(to_fill)=-999;
b.partnerchoice_OffsetTime(to_fill)=-999;
b.displaychoice_OnsetTime(to_fill)=-999;
b.displaychoice_OffsetTime(to_fill)=-999;
b.firstFixation_OnsetTime(to_fill)=-999;
b.decisions(to_fill)=-999;
save('trust212856.mat','b');

 if not(exist('decisions'))
            share =~cellfun(@isempty,strfind(b.PartDecides,'share'));
            keep =~cellfun(@isempty,strfind(b.PartDecides,'keep'));
            noresponse = ~cellfun(@isempty,strfind(b.PartDecides,'noresponse'));
            b.decisions = zeros(192, 1);
            b.decisions(share) = 1;
            b.decisions(keep) = -1;
            b.decisions(noresponse) = 0;
            save('trust212857.mat','b');
 end
 
 