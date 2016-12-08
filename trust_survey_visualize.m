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

masterlist = dlmread('/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/hcscansubjs_all_ages.txt');

load(strcat(datalocation{1},'surveys/trust_group_survey.mat'));
files_survey = dir(strcat(datalocation{1},'surveys/','trust_survey*.mat'));
survey_ids = arrayfun(@(x) (x.name(isstrprop(x.name,'digit'))), files_survey, 'UniformOutput', false);
survey_ids = arrayfun(@(x) str2double(x), survey_ids, 'UniformOutput', false);
survey_ids = cell2mat(survey_ids);

rating_good = [];
rating_bad=[]
rating_neutral=[];
rating_computer = [];
for ct = 1:length(surveys_exist)
        rating_good = [rating_good group_survey_table{(group_survey_table.ID == surveys_exist(ct)...
            & strcmp(group_survey_table.Trustee, 'good') & strcmp(group_survey_table.Question, 'trustworthy')),'Rating_Pre'}];
        rating_bad = [rating_bad group_survey_table{(group_survey_table.ID == surveys_exist(ct)...
            & strcmp(group_survey_table.Trustee, 'bad') & strcmp(group_survey_table.Question, 'trustworthy')),'Rating_Pre'}];
        rating_neutral = [rating_neutral group_survey_table{(group_survey_table.ID == surveys_exist(ct)...
            & strcmp(group_survey_table.Trustee, 'neutral') & strcmp(group_survey_table.Question, 'trustworthy')),'Rating_Pre'}];
        rating_computer = [rating_computer group_survey_table{(group_survey_table.ID == surveys_exist(ct)...
            & strcmp(group_survey_table.Trustee, 'computer') & strcmp(group_survey_table.Question, 'trustworthy')),'Rating_Pre'}];
end

ave_good = mean(rating_good);
ave_bad = mean(rating_bad);
ave_neutral = mean(rating_neutral);
ave_computer = mean(rating_computer);
figure(1)
ratings = [rating_good; rating_bad; rating_neutral; rating_computer];
boxplot(ratings','notch','on','labels',{'good','bad','neutral','computer'});
figure(2)
subplot(2,2,1);
histogram(rating_good);
title('good');
subplot(2,2,2);
histogram(rating_bad);
title('bad');
subplot(2,2,3);
histogram(rating_neutral);
title('neutral');
subplot(2,2,4);
histogram(rating_computer);
title('computer');
figure(3)
scatter(rating_good+0.5*rand(size(rating_good)), rating_bad);


