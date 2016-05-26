%% this script is for removing unwanted subjects (e.g., IDEATORS) from the dataset
% and for calculating IQR for RTs (to filter the outliers in SPSS)
all_ids=unique(behavior_n_demographics.subject);

S = table2struct(behavior_n_demographics, 'ToScalar', true);
non_ideators=cellfun(@(x) ~strcmp(x, 'IDEATOR'),S.pattype,'UniformOutput',false);
non_ideators=cell2mat(non_ideators);

ids = unique(S.subject);

S.iqr = double(zeros(length(S.subject),1));
S.q1 = double(zeros(length(S.subject),1));
S.q3 = double(zeros(length(S.subject),1));
S.upperL = double(zeros(length(S.subject),1));
S.lowerL = double(zeros(length(S.subject),1));
S.pt_decision = double(zeros(length(S.subject),1));
pt_decisions = S.t_decision;
pt_decisions2 = [0; S.t_decision(1:end-1)];
S.pt_decision = pt_decisions2;
for ct=1:length(ids)
    current=ids(ct);
    indices = S.subject==current;
    IQR = iqr(S.decision_RT(indices));
    S.iqr(indices) = IQR;
    S.q1(indices) = prctile(S.decision_RT(indices),25);
    S.q3(indices) = prctile(S.decision_RT(indices),75);
    S.upperL(indices) = S.q3(indices)+1.5*IQR;
    S.lowerL(indices) = S.q1(indices)-1.5*IQR;
end
indices_bt = S.trialnum==1;
S.pt_decision(indices_bt) = 0;
indices_trustee = (S.exchange==1);
S.pt_decision(indices_trustee) = 0;


T = struct2table(S);
writetable(T, '/Users/polinavanyukov/Box Sync/Project Trust Game/analyses/poster/all_behavior_demos');