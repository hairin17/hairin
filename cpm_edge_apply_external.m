clear;

%% read resting fmri_bul
no_nodes = 268
no_subs = 121

for i = 1:121;
    load(['resultsROI_Subject' num2str(i,'%03d') '_Condition001.mat']);
    all_mats_rs(:,:,i) = Z(:,1:no_nodes);
end 

% for i = 10:no_subs;
%     load(['resultsROI_Subject0' num2str(i) '_Condition001.mat']);
%     all_mats_rs(:,:,i) = Z(:,1:no_nodes);
% end

all_mats_rs(isnan(all_mats_rs)) = 0;
all_mats_rs = all_mats_rs(:,:,[1:8, 10:13, 15:28, 30:44, 46:60,62:68, 70, 72:73, 75:85, 87:95, 97:98, 99:121]);  %n 121 => 111 / 10 excluded motion (average rp subt >0.15 either in session 1 or 2)

%% remove cereb_ power
% powercereb=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0];
% networkmask=find(powercereb==0);
% all_mats_c=all_mats(networkmask,networkmask,:);
% all_mats_c(isnan(all_mats_c))=0;

%shen atlas remove cerebellum
cereb=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0];
networkmask=find(cereb==0); % Exclude cerebelum 
all_mats=all_mats_rs(networkmask,networkmask,:);
all_mats(isnan(all_mats))=0;


%% load behavioral data 
all_demo2=readtable('D:\matlab_code\cpm\cov.xlsx');
all_behav = table2array(all_demo2(:,52)); 

%% set structure 
for leftout = 1:no_subs;
    no_sub = size(all_mats,3);
    no_node = size(all_mats,1);
    behav_pred_pos = zeros(no_sub,1);
    behav_pred_neg = zeros(no_sub,1);

    train_mats = all_mats; % check
    train_mats(:,:,leftout) = [];

    train_behav = all_behav; %check
    train_behav(leftout) = [];

%% load selected overlapping edges from tsFC 
load selected_edges_ts_100.mat
pos_mask = zeros(no_node, no_node);
neg_mask = zeros(no_node, no_node);

pos_edges = overlap_mask_pos
neg_edges = overlap_mask_neg

pos_mask(pos_edges) = 1;
neg_mask(neg_edges) = 1; %% 여기 pos_edges => overlap으로 수정

% pos_mask = sum_mask_pos

%% number of selected feartures
    FeatSel_no_pos(leftout)=length(pos_edges);    
    FeatSel_no_neg(leftout)=length(neg_edges); 
 
%% sum all edges using pos/neg mask from ts FC
train_mats = all_mats;

train_sumpos = zeros(111,1); %train mat size external validation sample size로 수정
train_sumneg = zeros(111,1); 

for ss = 1:size (train_sumpos);
        train_sumpos(ss) = sum(sum(train_mats(:,:,ss).*overlap_loo_pos))/2;
        train_sumneg(ss) = sum(sum(train_mats(:,:,ss).*overlap_loo_neg))/2;
    end

%     for ss = 1:size (train_sumpos);
%         train_sumpos(ss) = sum(sum(train_mats(:,:,ss).*pos_mask))/2;
%         train_sumneg(ss) = sum(sum(train_mats(:,:,ss).*neg_mask))/2;
%     end

%% glm
test_mat = all_mats;

b = regress(train_behav, [train_sumpos, train_sumneg, ones(no_subs-1,1)]);
    test_mat = all_mats(:,:,leftout);
    test_sumpos = sum(sum(test_mat.*pos_mask))/2;
    test_sumneg = sum(sum(test_mat.*neg_mask))/2;

    behav_pred_glm(leftout) = b(1)*test_sumpos + b(2)*test_sumneg + b(3);
    
end

%% plot 
FeatSel_no = FeatSel_no_pos+FeatSel_no_neg;
FeatSel_perc = mean(FeatSel_no) / ((no_node)*(no_node-1)/2) * 100;
fprintf('\n Average number of feature selected = %6.3f among  %6.3f',mean(FeatSel_no),no_node*(no_node-1)/2);
fprintf('\n Percentage of feature selected = %6.3f \n',FeatSel_perc);
[R_glm, P_glm] = corr(behav_pred_glm',all_behav)

strglm = ['r = ', num2str(R_glm),', P = ', num2str(P_glm), ' , Feat = ', num2str(mean(FeatSel_no_pos+FeatSel_no_neg))];

figure(1); plot(behav_pred_glm,all_behav,'r.');lsline; xlabel('Predicted behavior');ylabel('Observed behavior')
text(mean(behav_pred_glm)-std(behav_pred_glm)*.9,mean(all_behav)-std(behav_pred_glm)*.9,strglm,'HorizontalAlignment','left');
axis([min(all_behav),max(all_behav), min(all_behav),max(all_behav)]);






