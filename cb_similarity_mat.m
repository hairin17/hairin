%% condition
% cyberball_condition
% Condition001 = rest
% Condition002 = inclusion
% Condition003 = exclusion
% Condition004 = re-inclusion
clear all;

%% settings---------------------------
for condition = 2:4    
cwd = 'E:\04_kshap2018\!analysis\stats_cb_95\conn_cb95\results\firstlevel\pow' %change shen, pow directory 
no_nodes = 264 %pow = 264; shen 268
no_subs = 95
folder = 'E:\04_kshap2018\!analysis\stats_cb_95\conn_cb95\pow_mat'%name
saveDir = 'E:\04_kshap2018\!analysis\stats_cb_95\conn_cb95\pow_mat'%name
filename = ['con_' num2str(condition)]
%------------------------------------
%% FC_matrix_extraction
cd (cwd)
[all_mats]=all_mats_extract(no_nodes,no_subs,condition);
all_mats(isnan(all_mats))=0;
% % shen cereb_rmove; network sorting
% shen_crb = load('shen_crb.mat');
% shen_crb = struct2array(shen_crb);
% networkmask = find(shen_crb == 0); % Exclude cerebellum
% all_mats = all_mats(networkmask,networkmask,:);
% all_mats(isnan(all_mats))=0;

% pow cereb_rmove; network sorting
pow_crb = load('E:\04_kshap2018\!analysis\stats_cb_95\conn_cb95\pow_mat\power_crb.mat');
pow_crb = struct2array(pow_crb);
networkmask = find(pow_crb == 0); % Exclude cerebellum
all_mats = all_mats(networkmask,networkmask,:);
all_mats(isnan(all_mats))=0;

%save_outputs
save (fullfile(saveDir,filename), 'all_mats')

cd(saveDir)
end

%% sub_network_shen
% load whole brain mat
cd 'E:\04_kshap2018\!analysis\stats_cb_95\conn_cb95\shen_mat'
for conNum = 1:4
load(['con_' num2str(conNum)])

folder = 'E:\04_kshap2018\!analysis\stats_cb_95\conn_cb95\shen_mat'
saveDir = 'E:\04_kshap2018\!analysis\stats_cb_95\conn_cb95\shen_mat'
filename = ['ment_con_' num2str(conNum)] %change file name as network mask

% network masks
net = table2array(readtable('shen_socog_rewment.xlsx'))
indices = find(net(:,2)==1);
net(indices,:) = [];
ment = net(:,3)';
rew = net(:,4)';
networkmask = find(ment == 1); % change ment or rew
all_mats = all_mats(networkmask,networkmask,:);
all_mats(isnan(all_mats))=0;
save (fullfile(saveDir,filename), 'all_mats')

cd(saveDir)

end

%% sub_network_power
% load whole brain mat
cd 'E:\04_kshap2018\!analysis\stats_cb_95\conn_cb95\pow_mat'
for conNum = 2:4
load(['con_' num2str(conNum)])

folder = 'E:\04_kshap2018\!analysis\stats_cb_95\conn_cb95\pow_mat'
saveDir = 'E:\04_kshap2018\!analysis\stats_cb_95\conn_cb95\pow_mat'
filename = ['ment_con_' num2str(conNum)] %change file name as network mask

% network masks
ment = struct2array(load('pow_ment'))
networkmask = find(ment == 1); % change ment or rew
all_mats = all_mats(networkmask,networkmask,:);
all_mats(isnan(all_mats))=0;
save (fullfile(saveDir,filename), 'all_mats')

cd(saveDir)

end

%% vectorize
for i = 2:4 
saveDir = 'E:\04_kshap2018\!analysis\stats_cb_95\conn_cb95\pow_mat'%change name as atlas
filename = ['ment_vec_' num2str(i)] %change name as network mask
load (['ment_con_' num2str(i)]) %change name as network mask
[all_vec] = lowervec(all_mats);

save(fullfile(saveDir,filename), 'all_vec')
end

%% matrix correlation_fisher'z
%con_nums I(2) E(3) RE(4)
j = 2
k = 4 %3,4
saveDir = 'E:\04_kshap2018\!analysis\stats_cb_95\conn_cb95\pow_mat'%change name
filename = ['ment_cor_' num2str(j) '_' num2str(k) ] %change name as network mask

vec_j  = struct2array(load(['ment_vec_' num2str(j)])); %change name as network mask
vec_k  = struct2array(load(['ment_vec_' num2str(k)])); %change name as network mask
[corr_val] = vec_corr(vec_j,vec_k);
corr_val = corr_val';
save(fullfile(saveDir,filename), 'corr_val');


%% correlation w/ behavior
 [r_mat,p_mat] = corr(train_vcts',train_behav); 


