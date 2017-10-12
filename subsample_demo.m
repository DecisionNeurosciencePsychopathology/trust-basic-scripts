
%read in database files
[num,txt,raw]=xlsread('demos07192016.xlsx');

beha=load('beha_behavior_data.mat');
scan=load('scan_behavior_data.mat');

beha.ready=cat(2,beha.x,ones(size(beha.x,1),1)); %column 13: 1 = behavior
scan.ready=cat(2,scan.x,zeros(size(scan.x,1),1));%column 13: 0 = scan
x = cat(1,beha.ready,scan.ready);
save('beha_n_scan.mat','x');

%create empty array for subsample demographics
demoarray=cell(0,size(raw,2));

all_ids = cell2mat(raw(2:size(raw,1),1));
sample_ids = x(:,1);
id_double=0;
for i=2:size(raw,1)
    id = raw{i,1};
    if not(id==id_double)
        %find if this id is in the datatable
        if not(isempty(find(sample_ids==id,1)))
            %copy the relevant row to the demoarray
            demoarray=cat(1,demoarray,raw(i,:));
        end
        id_double=id;
    end
end

save('demoarray.mat','demoarray');

%converting the data (array of doubles) from participants to a cell array
test = num2cell(x,size(x));

%match the ids from test and demoarray
%using the info from demo array, fill in the columns in the test array for
%that id
[d_rows, d_columns]=size(demoarray);
[test_rows, test_columns] = size(test);
for i=1:d_rows
    id = demoarray{i};
    id_rows = find(cell2mat(test(:,1))==id);
    test(id_rows,test_columns+1:test_columns+d_columns)=repmat(demoarray(i,:),size(id_rows,1),1);
end

save('all_behavior_demo.mat','test');

%may be different depending on what type of data is being combined with
%demographics
Names = {'subject','trialnum','trustee','exchange','reward_schedule',...
    'exch2','s_decision','t_decision','decision_Onset','decision_RT',...
    'feedback_Onset','feedback_Offset','beha1scan0','id2','consent_age',...
    'today_age','group','group1245','group12467','group12','sext','racet','maxleth',...
    'maxofEDUC','experiment_t','qualityt','ideation_screen','ideation_total','WTAR_scal',...
    'MDRS','EXIT','PHYS_ILL_CIS','HRSD_no_suic','NEUROTICISM', 'EXTRAVERSION','OPENNESS',...
    'AGREEABLENESS','CONSCIENTIOUSNESS','BECK_HOPELSS','INTENT_PLAN', 'INTENT_TOTAL',...
    'SPSI_ICS','BIS_NONPLAN','BIS_MOTOR','BIS_COGN', 'BIS_TOTMEAN','ARSTOTAL','IIP15INTSENS',...
    'IIP15INTAMBV','IIP15AGRESS','UPPSP_NEG_URG','UPPSP_POS_URG','UPPSP_LACKOFPREMED',...
    'UPPSP_LACKOFPERSEV','WCSTERR','WERR','ATHF_STRENGTH','ATHF_CUMSTRENGTH',...
    'SUBSTANCE_LIFE','SUBSTANCE_CURR','ANXIETY_LIFE', 'ANXIETY_CURR'};

behavior_n_demographics=cell2table(test,'VariableNames',Names);
filename = '/Users/polinavanyukov/Box Sync/Project Trust Game/analyses/poster/group_behavior_demos'; 
%no_ideators = behavior_n_demographics(not(cell2mat(behavior_n_demographics.group1245)=='4'),:);
%writetable(no_ideators,filename);
writetable(behavior_n_demographics, '/Users/polinavanyukov/Box Sync/Project Trust Game/analyses/manuscript/all_behavior_demos');