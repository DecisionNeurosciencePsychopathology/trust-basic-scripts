function b = trustsurvey(id)

%reading in one subject's survey pre and post data;

%Quick username check, and path setting, this may have to change depending
%on the machine you are currently working on!
os = computer;
if strcmp(os(1:end-2),'PCWIN')
    datalocation = glob('?');
else
    [~, me] = system('whoami');
    me = strtrim(me);
    if strcmp(me,'polinavanyukov')==1
        datalocation = glob('/Users/polinavanyukov/Box Sync/Suicide studies/data/');
    else
        datalocation = glob('?');
    end
end


filename = sprintf('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/surveys/trust_survey%d.mat',id);

%% Find the eprime file - MODIFY PATHS IF NEEDED

cd(datalocation{1})

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
subdir=sprintf('%d',id);

cd(subdir)
files = dir('trust*.mat');
num_of_subjects = length(files);

% fname_1 = sprintf('trust_10212014_prescan-%s-1.txt', id5(1:end));
% fname_2 = sprintf('trust_10212014_postscan-%s-1.txt', id5(1:end));

fname_1 = dir('*prescan*.txt');
fname_2 = dir('*postscan*.txt');

fields = {'procedure','Name','Question','Rating'};
startoffset = 2; % trust survey: from procedure occurrence to the beginning of the block containing relevant info
endoffset = 18;   % trust survey: from procedure occurrence to the end of the block containing relevant info
b = read_in_trust(fname_1.name,'procedure', fields, startoffset, endoffset);
startoffset = 2; % trust survey: from procedure occurrence to the beginning of the block containing relevant info
endoffset = 18;   % trust survey: from procedure occurrence to the end of the block containing relevant info
fields = {'procedure','Name', 'Question','Rating'};
c = read_in_trust(fname_2.name,'procedure', fields, startoffset, endoffset);
b.trustee = strrep(b.procedure,'ratings','');
b.trustee = cellstr(b.trustee([1:8]));
b.Rating = b.Rating([1:8],:);
b.Question = b.Question([1:8]);
c.trustee = strrep(c.procedure,'ratings','');
c.trustee = c.trustee([1:24],:);
c.Rating = c.Rating([1:24],:); 
c.Question = c.Question([1:24],:);
T = table(repmat(id,size(b.trustee)), b.trustee, b.Question, b.Rating, 'VariableNames',{'ID','Trustee', 'Question', 'Rating_Pre'});
for index=1:2:size(T.Trustee)
    rows = strcmp(c.trustee,T.Trustee(index));
    d(index:index+1)=c.Rating(rows);
end
T.Rating_Post = d';
save(filename,'T'); 
return
