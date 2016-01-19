%%plot share probabilities by trustee, reward schedule sequence and trial

behav_demo.s_decision10 = behav_demo.s_decision==1;
bad = behav_demo.s_decision10(behav_demo.trustee==2);
good = behav_demo.s_decision10(behav_demo.trustee==1);
neutral = behav_demo.s_decision10(behav_demo.trustee==3);
comp = behav_demo.s_decision10(behav_demo.trustee==4);


t = table2struct(behav_demo, 'ToScalar','1');
idlist = unique(behav_demo.subject);
mat_data = zeros(length(unique(behav_demo.trialnum)),7,length(idlist));
for ct = 1:length(idlist)
    sub = idlist(ct);
    
    %% code block order 1 = 25 first, 2 = 88 first
    if t.reward_schedule(t.subject==sub & t.trialnum==17)==25 
    t.block_order(t.subject==sub & t.trialnum<49)=1;
    elseif t.reward_schedule(t.subject==sub & t.trialnum==17)==88
    t.block_order(t.subject==sub & t.trialnum<49)=2;
    end
    if t.reward_schedule(t.subject==sub & t.trialnum==65)==25
    t.block_order(t.subject==sub & t.trialnum>48 & t.trialnum<97)=1;
    elseif t.reward_schedule(t.subject==sub & t.trialnum==65)==88
    t.block_order(t.subject==sub & t.trialnum>48 & t.trialnum<97)=2;
    end   
    if t.reward_schedule(t.subject==sub & t.trialnum==113)==25
    t.block_order(t.subject==sub & t.trialnum>96 & t.trialnum<145)=1;
    elseif t.reward_schedule(t.subject==sub & t.trialnum==113)==88
    t.block_order(t.subject==sub & t.trialnum>96 & t.trialnum<145)=2;
    end   
    if t.reward_schedule(t.subject==sub & t.trialnum==161)==25
    t.block_order(t.subject==sub & t.trialnum>144 )=1;
    elseif t.reward_schedule(t.subject==sub & t.trialnum==161)==88
    t.block_order(t.subject==sub & t.trialnum>144)=2;
    end   
    
    sub_data = [t.subject(t.subject==sub) t.trialnum(t.subject==sub) t.trustee(t.subject==sub) t.exchange(t.subject==sub) ...
       t.reward_schedule(t.subject==sub) t.s_decision10(t.subject==sub) t.block_order(t.subject==sub)'];
     
    mat_data(:,:,ct) = sub_data;
end
t.block_order = t.block_order';
good25 = [];
good88 = [];
bad25 = [];
bad88 = [];
neutral25 = [];
neutral88 = [];
computer25 = [];
computer88 = [];
for ct = 1:length(idlist)
    sub = idlist(ct);
    foo = t.s_decision10(t.subject==sub & t.trustee==1 & t.block_order==1);
    if length(foo)>1
    good25 = [good25 foo(1:47)];
    end
    foo = t.s_decision10(t.subject==sub & t.trustee==1 & t.block_order==2);
    if length(foo)>1
    good88 = [good88 foo(1:47)];
    end
    foo = t.s_decision10(t.subject==sub & t.trustee==2 & t.block_order==1);
    if length(foo)>1
    bad25 = [bad25 foo(1:47)];
    end
    foo = t.s_decision10(t.subject==sub & t.trustee==2 & t.block_order==2);
    if length(foo)>1
    bad88  = [bad88 foo(1:47)];
    end  
    foo = t.s_decision10(t.subject==sub & t.trustee==3 & t.block_order==1);
    if length(foo)>1
    neutral25  = [neutral25 foo(1:47)];
    end
    foo = t.s_decision10(t.subject==sub & t.trustee==3 & t.block_order==2);
    if length(foo)>1
    neutral88 = [neutral88 foo(1:47)];
    end
    foo = t.s_decision10(t.subject==sub & t.trustee==4 & t.block_order==1);
    if length(foo)>1
    computer25 = [computer25 foo(1:47)];
    end
    foo = t.s_decision10(t.subject==sub & t.trustee==4 & t.block_order==2);
    if length(foo)>1
    computer88  = [computer88 foo(1:47)];
    end
end

%% plot all trustees
figure(3); clf;
subplot(2,1,1)
hold on; plot(smooth(mean(good25,2))); plot(smooth(mean(bad25,2)));plot(smooth(mean(neutral25,2)));plot(smooth(mean(computer25,2)));
legend('Good', 'Bad', 'Neutral', 'Computer');
title('50% - 25% - 88%');
subplot(2,1,2)
hold on; plot(smooth(mean(good88,2))); plot(smooth(mean(bad88,2)));plot(smooth(mean(neutral88,2)));plot(smooth(mean(computer88,2)));
legend('Good', 'Bad', 'Neutral', 'Computer');
title('50% - 88% - 25%');
%% plot good and bad by contingency order only -- changed to smoothed version
figure(4); clf; 
subplot(2,1,1)
hold on; plot(smooth(mean(good25,2))); plot(smooth(mean(bad25,2)));
legend('Good', 'Bad');
title('Good vs. bad trustee, 50% - 25% - 88%')
subplot(2,1,2)
hold on; plot(smooth(mean(good88,2))); plot(smooth(mean(bad88,2)));
title('Good vs. bad trustee 50% - 88% - 25%')
legend('Good', 'Bad');


%% separate and plot choices according to the order of contingencies
choice = squeeze(mat_data(:,6,:));
block_order = squeeze(mat_data(:,7,:));
first25choice = choice(49:144,:);
first25choice1 = first25choice(1:48,:); 
first25choice2 = first25choice(49:end,:);
first25choice3d = cat(3,first25choice1,first25choice2);
first88choice = [choice(1:48,:);  choice(145:end,:)];
first88choice1 = first88choice(1:48,:); 
first88choice2 = first88choice(49:end,:);
first88choice3d = cat(3,first88choice1,first88choice2);

figure(1); clf; hold on; plot(mean(mean(first25choice3d,3),2)); plot(mean(mean(first88choice3d,3),2)); hold off; 

% figure(1); clf; 
% subplot(2,1,1)
% plot(smooth(mean(first25choice(1:48,:),2)));hold on;  plot(smooth(mean(first88choice(1:48,:),2)));
% subplot(2,1,2)
% plot(smooth(mean(first25choice(49:end,:),2))); hold on; plot(smooth(mean(first88choice(49:end,:),2)));


%% separate choices according to trustee and order of contingencies
trustee = squeeze(mat_data(:,3,:));
trustee1 = trustee==1;
trustee2 = trustee==2;
trustee3 = trustee==3;
trustee4 = trustee==4;
choiceGood=zeros(48,67);

% 
% for ind=1:67
%     choiceGood(:,ind)=choice(trustee1(:,ind),ind); %mismatched dimensions?
% end



mat_fields = {'ID' 'trial' 'trustee' 'exchange' 'contingency' 's_decision' 'block_order'};
    