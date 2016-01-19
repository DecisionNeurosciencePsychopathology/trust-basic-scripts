
%read in database files
[num,txt,raw]=xlsread('splashDemo2.xlsx');

beha=load('allsubjects.mat');
scan=load('allsubjects2.mat');

beha.ready=cat(2,beha.x,ones(size(beha.x,1),1));
scan.ready=cat(2,scan.x,zeros(size(scan.x,1),1));
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
    'feedback_Onset','feedback_Offset','beha1scan0','id2','study',...
    'initials','consent','consent_age','today_age','dateTermin','pattype',...
    'comment','group1245','group12467','sext','ethnicity','racet','maxleth',...
    'maxofEDUC'};

behav_demo=cell2table(test,'VariableNames',Names);
