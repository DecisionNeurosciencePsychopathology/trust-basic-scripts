function b = trustbehavior(id)
% Polina Vanyukov, Jon Wilson & Alex Dombrovski
% 2014-05: revision of the program for the trust game paradigm
% Program that reads in a single data file



data_dir_str= '/Users/polinavanyukov/Box Sync/Suicide studies/data/';
filename = sprintf('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/trust%d.mat',id);

%% Find the eprime file - MODIFY PATHS IF NEEDED
%cd('/Users/polinavanyukov/Box Sync/Project Trust Game/eprime')
%cd(['C:\Users\',username ,'\Box Sync\Suicide studies\data'])
%cd(['/Volumes/bek/trust_analyses/temp_trust_data'])
cd(data_dir_str)

if id>209999 && id<300000
    id5=sprintf('%d',(id-200000));
elseif id<209999 && id>200999
    id5=sprintf('0%d',(id-200000));
elseif id<200999 && id>200099
    id5=sprintf('00%d',(id-200000));
elseif id<200099 && id>200009
    id5=sprintf('000%d',(id-200000));
elseif id<200009 && id>200000
    id5=sprintf('0000%d',(id-200000));
elseif id<100000
    id5=sprintf('%d',id);
elseif id>109999 && id<200000
    id5=sprintf('%d',(id-100000));
elseif id>880000
    stringid=num2str(id);
    id5=stringid(3:end);
else
    stringid=num2str(id);
    id5=stringid(2:end);
end
subdir=dir(sprintf('%d*',id));

cd(subdir.name)

file = dir(sprintf('trust*%d*.txt', id));
if isempty(file)
  file = dir(sprintf('trust*%s*.txt', id5));  
end

fname = file.name;
%fname = sprintf('trust_05222015_scan_rs3-%s-1.txt', id5);
%fname = 'trust_05192015_scan_rs4-0064-1.txt';
%Added more fields for timing checking
debug_fields = {'displaychoice.OnsetTime', 'displaychoice.OffsetTime', ...
    'firstFixation.OnsetTime'};

%JW: Might be better to put these in order from txt file
fields = {'ConditionOrder','identity','Running','partnerchoice.RESP',...
     'TrusteeDecides','TrialNumber','exchangeNum','Reversal', ...
     'partnerchoice.OnsetTime','partnerchoice.RTTime','partnerchoice.RT'...
     ,'partnerchoice.Duration','partnerchoice.RESP','PartDecides','displaychoice.Duration',...
     'outcome.OnsetTime','outcome.OffsetTime','Order_RS', ...
     'ITIfixation.OnsetTime','ITIfixation.Duration',...
     'ITIfixation.OffsetTime', 'ITIjitter','partnerchoice.OffsetTime', debug_fields{:}};
startoffset = 14; % trust: from TrialProc occurrence to the beginning of the block containing relevant info
endoffset = 32;   % trust: from TrialProc occurrence to the end of the block containing relevant info
b = read_in_trust(fname,'TrialProc', fields, startoffset, endoffset);

trials = 192;
start = 0;

%identifying the beginning of nonpractice trials
for i=1:trials
    try
        if (b.TrialNumber(i) < 0)
            start = i;
        end
    catch
        throw_trust_error(id)
        %return
    end
end

if length(b.TrialNumber) > 192
    b = structfun(@(x) (x(13:end)), b, 'UniformOutput', false);
end

%Decisions share = 1; keep = -1; no response = 0;
share =~cellfun(@isempty,strfind(b.PartDecides(start+1:trials),'share'));
keep =~cellfun(@isempty,strfind(b.PartDecides(start+1:trials),'keep'));
noresponse = ~cellfun(@isempty,strfind(b.PartDecides(start+1:trials),'noresponse'));
b.decisions = zeros(trials-start, 1);
b.decisions(share) = 1;
b.decisions(keep) = -1;
b.decisions(noresponse) = 0;

trustee1 = 'none';
trustee2 = 'none';
trustee3 = 'none';
trustee4 = 'none';


%get identity for each phase of the task
for ind=start+1:trials
    if strcmp(trustee1, 'none')||strcmp(trustee1, b.identity(ind))
        trustee1 = b.identity(ind);
    elseif strcmp(trustee2, 'none')||strcmp(trustee2, b.identity(ind))
        trustee2 = b.identity(ind);
    elseif strcmp(trustee3, 'none')||strcmp(trustee3, b.identity(ind))
        trustee3 = b.identity(ind);
    else
        trustee4 = b.identity(ind);
    end
end

share1_50 = 0;
share2_50 = 0;
share3_50 = 0;
share4_50 = 0;

share1_88 = 0;
share2_88 = 0;
share3_88 = 0;
share4_88 = 0;

share1_25 = 0;
share2_25 = 0;
share3_25 = 0;
share4_25 = 0;

share1Total = 0;
share2Total = 0;
share3Total = 0;
share4Total = 0;

b.trustee1 = {trustee1, share1_50, share1_88, share1_25, share1Total};
b.trustee2 = {trustee2, share2_50, share2_88, share2_25, share2Total};
b.trustee3 = {trustee3, share3_50, share3_88, share3_25, share3Total};
b.trustee4 = {trustee4, share4_50, share4_88, share4_25, share4Total};

% %share responses for different phases of the task
for ind=start+1:trials
    if strcmp(trustee1, b.identity(ind)) && strcmp(b.PartDecides(ind), 'share')
        share1Total = share1Total + 1;
        if b.Reversal(ind) == 50
            share1_50 = share1_50 + 1;
        end
        if b.Reversal(ind) == 25
            share1_25 = share1_25+1;
        end
        if b.Reversal(ind) == 88
            share1_88 = share1_88+1;
        end
        b.trustee1 = {trustee1, share1_50, share1_88, share1_25, share1Total};
    elseif strcmp(trustee2, b.identity(ind)) && strcmp(b.PartDecides(ind), 'share')
        share2Total = share2Total + 1;
        if b.Reversal(ind) == 50
            share2_50 = share2_50 + 1;
        end
        if b.Reversal(ind) == 25
            share2_25 = share2_25+1;
        end
        if b.Reversal(ind) == 88
            share2_88 = share2_88+1;
        end
        b.trustee2 = {trustee2, share2_50, share2_88, share2_25, share2Total};
    elseif strcmp(trustee3, b.identity(ind)) && strcmp(b.PartDecides(ind), 'share')
        share3Total = share3Total + 1;
        if b.Reversal(ind) == 50
            share3_50 = share3_50 + 1;
        end
        if b.Reversal(ind) == 25
            share3_25 = share3_25+1;
        end
        if b.Reversal(ind) == 88
            share3_88 = share3_88+1;
        end
        b.trustee3 = {trustee3, share3_50, share3_88, share3_25, share3Total};
    elseif strcmp(trustee4, b.identity(ind)) && strcmp(b.PartDecides(ind), 'share')
        share4Total = share4Total + 1;
        if b.Reversal(ind) == 50
            share4_50 = share4_50 + 1;
        end
        if b.Reversal(ind) == 25
            share4_25 = share4_25+1;
        end
        if b.Reversal(ind) == 88
            share4_88 = share4_88+1;
        end
        b.trustee4 = {trustee4, share4_50, share4_88, share4_25, share4Total};
    end
end

save(filename, 'b');

%%

return

%JW: Functions
%If erroneous subjs creep up
function throw_trust_error(id)
global skip_subjs
    skip_subjs = [skip_subjs  id];
    fprintf(['\n\nTo many files for this subj review for later!',...
        'OR to something is wrong with the trials!\n'])
    save skip_subjs skip_subjs



