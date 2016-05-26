%Pulling all the surveys into one group survey and performing some basic
%operations
%Quick username check, and path setting, this may have to change depending
%on the machine you are currently working on!
os = computer;
if strcmp(os(1:end-2),'PCWIN')
    datalocation = glob('?');
else
    [~, me] = system('whoami');
    me = strtrim(me);
    if strcmp(me,'polinavanyukov')==1
        datalocation = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/surveys/');
    else
        datalocation = glob('?');
    end
end

cd(datalocation{1})

files = dir('trust_survey*.mat');
num_of_subjects = length(files);
group_survey_table = table([],[],[],[],[],'VariableNames',{'ID','Trustee', 'Question', 'Rating_Pre','Rating_Post'});

for index=1:num_of_subjects
    load(files(index).name);
    group_survey_table = [group_survey_table;T];
    
end
save('trust_group_survey','group_survey_table'); 

b = group_survey_table;
func = @mean;
varMeans = varfun(func, b, 'InputVariables', {'Rating_Pre','Rating_Post'}, 'GroupingVariables',{'Trustee','Question'});
