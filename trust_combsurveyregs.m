clear variables;
close all;

%Quick username check, and path setting, this may have to change depending
%on the machine you are currently working on!
os = computer;
if strcmp(os(1:end-2),'PCWIN')
    datalocation = glob('?');
else
    [~, me] = system('whoami');
    me = strtrim(me);
    if strcmp(me,'polinavanyukov')==1
        datalocation = glob('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/beha_behavior/');
    else
        datalocation = glob('?');
    end
end

cd(datalocation{1});
files_regs = dir('trust*.mat');
num_of_subjects = length(files_regs);
masterlist = dlmread('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/hcscansubjs_all_ages.txt');

load(strcat(datalocation{1},'surveys/trust_group_survey.mat'));
files_survey = dir(strcat(datalocation{1},'surveys/','trust_survey*.mat'));
survey_ids = arrayfun(@(x) (x.name(isstrprop(x.name,'digit'))), files_survey, 'UniformOutput', false);
survey_ids = arrayfun(@(x) str2double(x), survey_ids, 'UniformOutput', false);
survey_ids = cell2mat(survey_ids);

surveys_exist = intersect(masterlist, survey_ids);
surveys_NTM = setxor(masterlist, surveys_exist);

for ct = 1:length(surveys_NTM)
    trustsurvey(surveys_NTM(ct));
end
trust_group_survey();
load(strcat(datalocation{1},'surveys/trust_group_survey.mat'));

for ct = 1:length(surveys_exist)
    survey_ID = surveys_exist(ct);
    if not(survey_ID==881209||survey_ID ==217909)
        load(strcat(datalocation{1},'trust',num2str(survey_ID),'.mat'));
        b.rating_pre.good = group_survey_table{(group_survey_table.ID == survey_ID...
            & strcmp(group_survey_table.Trustee, 'good') & strcmp(group_survey_table.Question, 'trustworthy')),'Rating_Pre'};
        b.rating_post.good = group_survey_table{(group_survey_table.ID == survey_ID...
            & strcmp(group_survey_table.Trustee, 'good') & strcmp(group_survey_table.Question, 'trustworthy')),'Rating_Post'};
        b.rating_pre.bad = group_survey_table{(group_survey_table.ID == survey_ID...
            & strcmp(group_survey_table.Trustee, 'bad') & strcmp(group_survey_table.Question, 'trustworthy')),'Rating_Pre'};
        b.rating_post.bad = group_survey_table{(group_survey_table.ID == survey_ID...
            & strcmp(group_survey_table.Trustee, 'bad') & strcmp(group_survey_table.Question, 'trustworthy')),'Rating_Post'};
        b.rating_pre.neutral = group_survey_table{(group_survey_table.ID == survey_ID...
            & strcmp(group_survey_table.Trustee, 'neutral') & strcmp(group_survey_table.Question, 'trustworthy')),'Rating_Pre'};
        b.rating_post.neutral = group_survey_table{(group_survey_table.ID == survey_ID...
            & strcmp(group_survey_table.Trustee, 'neutral') & strcmp(group_survey_table.Question, 'trustworthy')),'Rating_Post'};
        b.rating_pre.computer = group_survey_table{(group_survey_table.ID == survey_ID...
            & strcmp(group_survey_table.Trustee, 'computer') & strcmp(group_survey_table.Question, 'trustworthy')),'Rating_Pre'};
        b.rating_post.computer = group_survey_table{(group_survey_table.ID == survey_ID...
            & strcmp(group_survey_table.Trustee, 'computer') & strcmp(group_survey_table.Question, 'trustworthy')),'Rating_Post'};
    end
    upd_filename = strcat('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/beha_behavior/trust',num2str(survey_ID),'.mat');
    save(upd_filename, 'b');
end


