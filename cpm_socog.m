
%% connectome-based predictive modeling
%Shen et al., 2017
%all_mats task, resting check
%path check
clear;

%% read behavioral data;
all_demo2=readtable('E:\04_kshap2018\analysis\cb\cpm\cbcpm_cov.xlsx');

%% fMRI data Inputs
%cbcond
% Condition001 = inc
% Condition002 = exc
% Condition003 = reinc
% Condition004 = fix1
% Condition005 = fix2

%check matlab current directory;
% of ROIs Shen(268), HO(132), HCP(32)  Power(264)
no_nodes = 268 
no_subs = 58

for i=1:9;
    load (['resultsROI_Subject00' num2str(i) '_Condition002.mat']); %conn_condition number, 002 = I, 003 = C, 004 = IC
    all_mats(:,:,i)=Z(:,1:no_nodes);
end

for i=10:no_subs;
    load (['resultsROI_Subject0' num2str(i) '_Condition002.mat']); %conn_condition number, 002 = I, 003 = C, 004 = IC
    all_mats(:,:,i)=Z(:,1:no_nodes);

end
all_mats(isnan(all_mats))=0;
% 
% powercereb=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0];
% networkmask=find(powercereb==0);
% all_mats=all_mats(networkmask,networkmask,:);
% all_mats(isnan(all_mats))=0;
%  
% % fmri data inputs_resting; DIRECTORY!!!!!!!!
% check matlab current directory;
% # of ROIs Shen(268), HO(132), HCP(32)  Power(264)
% no_nodes = 268 
% no_subs = 58
% 
% for i=1:9;
%     load (['resultsROI_Subject00' num2str(i) '_Condition001.mat']); %rest
%     all_mats2(:,:,i)=Z(:,1:no_nodes);
% end
% 
% for i=10:no_subs;
%     load (['resultsROI_Subject0' num2str(i) '_Condition001.mat']); %rest
%     all_mats2(:,:,i)=Z(:,1:no_nodes);
% 
% end
% all_mats2(isnan(all_mats2))=0;
% % power 264 atlas remove cerebellum
% powercereb=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0];
% networkmask=find(powercereb==0);
% all_mats2=all_mats2(networkmask,networkmask,:);
% all_mats2(isnan(all_mats2))=0;
%  

%% remove cerebellum
%shen atlas remove cerebellum
cereb=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0];
networkmask=find(cereb==0); % Exclude cerebelum 
all_mats=all_mats(networkmask,networkmask,:);
all_mats(isnan(all_mats))=0;



% % power 264 atlas remove cerebellum
% powercereb=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0];
% networkmask=find(powercereb==0);
% all_mats2=all_mats2(networkmask,networkmask,:);
% all_mats2(isnan(all_mats2))=0;
%  
%% behavioral data selection
all_behav = table2array(all_demo2(:,9));  % 51(IICV) 52(CICV) 50(CSD)  49(ISD) 33(CRT) 34(IRT) 141(MEAN MOTION), select a column 

%% --------------options---------------- %%
%% including covariates
cov = table2array(all_demo2(:,[2:4]));
%% missing exclude
% missing=find(isnan(all_behav));
% all_behav(missing,:)=[];
% all_mats(:,:,missing)=[];
%% motion exclude 
% maxmotion= table2array(all_demo2(:,156)); 
% motionexc=find(maxmotion>3);
% all_behav(motionexc,:)=[];
% all_mats(:,:,motionexc)=[];
% vct=reshape(all_mats_shen_net,[],size(all_mats_shen_net,3));
% net_mean=mean(vct);
%% split data
% sex=table2array(all_demo2(:,10)); % female=0
% sexwhere=find(sex==1);  %to exclude
% all_behav(sexwhere,:)=[];
% all_mats(:,:,sexwhere)=[];
%% Random sampling data split
% id = 1:size(all_mats,3);
% testid=id;
% trainid=randsample(id,round(length(all_behav)*.667));
% testid(trainid)=[];
% test_mats=all_mats(:,:,testid);
% test_behav=all_behav(testid);
% all_mats=all_mats(:,:,trainid);
% all_behav=all_behav(trainid);

%% ------------- options end ------------------%%
%% modeling part 
%% feature selection
thresh = 0.01;
%
no_sub=size(all_mats,3);
%no_node = size(all_mats,1);

behav_pred_pos=zeros(no_sub,1);
behav_pred_neg=zeros(no_sub,1);

%
overlap_pos = zeros(227,25); % added to find overlap edges
overlap_neg = zeros(227,25);

for leftout = 1:no_sub;
    fprintf('\n Leaving out subj # %6.3f',leftout);

    % leave out subject from matrices and behavior
    
    train_mats = all_mats;% all_mats; %%mats2 = rest training
    train_mats = all_mats(:,[1 2 3	5 10 16	19	48	49	55	64	65	90	117	118	120	125	128	163	165	167	170	202	205	211],:);
    
    train_mats(:,:,leftout) = [];
    train_vcts=reshape(train_mats,[],size(train_mats,3)); %reshape(X,[M,N]) returns the M-by-N matrix 
    
    train_behav = all_behav;
    train_behav(leftout) = [];

%coorrelation with edges; correlation type
    % correlate all edges with behavior
%     [r_mat,p_mat] = corr(train_vcts',train_behav); 
%     r_mat=reshape(r_mat,no_node,no_node);
%     p_mat=reshape(p_mat,no_node,no_node);
    
    % correlate all edges with behavior using rank correlateion
%     [r_mat,p_mat]=corr(train_vcts',train_behav,'type','Spearman');
%     r_mat=reshape(r_mat,no_node,no_node);
%     p_mat=reshape(p_mat,no_node,no_node);    

     %correlate all edges using partial correlation
    train_cov=cov;
    train_cov(leftout,:)=[];
    [r_mat, p_mat] = partialcorr(train_vcts',train_behav,train_cov); %trainvcts 5675 = 227*25
    r_mat=reshape(r_mat,227,25);
    p_mat=reshape(p_mat,227,25);  
    

 %set threshold and define masks
    pos_mask = zeros(227,25);
    neg_mask = zeros(227,25);
    
    pos_edges = find(r_mat>0 & p_mat<thresh);
    neg_edges = find(r_mat<0 & p_mat<thresh);
    
    pos_mask(pos_edges)=1;
    neg_mask(neg_edges)=1;
    
%     pos_mask = pos_mask(:,[1 2	3 5 10 16	19	48	49	55	64	65	90	117	118	120	125	128	163	165	167	170	202	205	211]);
%     neg_mask = neg_mask(:,[1 2	3 5 10 16	19	48	49	55	64	65	90	117	118	120	125	128	163	165	167	170	202	205	211]);


%set overlapping masks
    overlap_pos = overlap_pos + pos_mask;
    overlap_neg = overlap_neg + neg_mask;
    
 % Number of features selected
    FeatSel_no_pos(leftout)=length(pos_edges);    
    FeatSel_no_neg(leftout)=length(neg_edges);   
    
 % get sum of all edges in TRAIN subs (divide by 2 to control for the fact that matrices are symmetric)
    train_sumpos=zeros(no_sub-1,1);
    train_sumneg=zeros(no_sub-1,1);
    for ss= 1:size(train_sumpos);
        train_sumpos(ss) = sum(sum(train_mats(:,:,ss).*pos_mask))/2;  
        train_sumneg(ss) = sum(sum(train_mats(:,:,ss).*neg_mask))/2;
    end
 % akin to computing the dot product btw an individual connectivity matrix
 % and the binary feature masks generated above. Alternatively, a sigmoid
 % function centered at the P value threshold can be used as a weighting
 % function, and a seighted sum can be used as the single-subject summary
 % value; this oviates a binary threshold 
 
 %-----------Sigmoidal weighting----------------%
%  pos_edges = find(r_mat > 0);
%  neg_edges = find(r_mat < 0);
%  
%  convert p threshold to r threshold 
%  T = tinv(thresh/2, no_sub-1-2);
%  R = sqrt(T^2/(no_sub-1-2+T^2));
%  
%  create a weighted mask using sigmoidal function
%  weight = 0.5, when corrlated = R/3;
%  weight = 0.88, when correlation = R;
%  pos_mask(pos_edges) = sigmf( r_mat(pos_edges), [3/R, R/3]);
%  neg_maask(neg_edges) = sigmf( r_mat(neg_edges), [-3/R, R/3]);
 %----------Sigmoida weighting-----------------%

%% build model on TRAIN subs & test model
%% ---------pos/neg edges separated----------%
    fit_pos = polyfit(train_sumpos, train_behav,1);
    fit_neg = polyfit(train_sumneg, train_behav,1);
%    fit_pos = polyfit(train_sumpos, train_behav,2);
%     fit_neg = polyfit(train_sumneg, train_behav,2);    

    fitted_coef_pos(leftout,1:2)=fit_pos;
    fitted_coef_neg(leftout,1:2)=fit_neg;
% run model on TEST sub
%     test_mat = all_mats2(:,:,leftout);%resting    
%     test_mat = all_mats(:,:,leftout); %all_mats(:,:,leftout); 
    test_mat = all_mats(:,[1 2 3	5 10 16	19	48	49	55	64	65	90	117	118	120	125	128	163	165	167	170	202	205	211],leftout);
    test_sumpos = sum(sum(test_mat.*pos_mask))/2;
    test_sumneg = sum(sum(test_mat.*neg_mask))/2;
    
    behav_pred_pos(leftout) = fit_pos(1)*test_sumpos + fit_pos(2);
    behav_pred_neg(leftout) = fit_neg(1)*test_sumneg + fit_neg(2);
%     behav_pred_pos(leftout) = fit_pos(1)*test_sumpos^2 + fit_pos(2)*test_sumpos + fit_pos(3);
%     behav_pred_neg(leftout) = fit_neg(1)*test_sumneg^2 + fit_neg(2)*test_sumneg + fit_neg(3);    
%---------pos/neg edges separated----------%
%% combining both pos and neg features (glm)    
    b = regress(train_behav, [train_sumpos, train_sumneg, ones(no_sub-1,1)]);
% run model on TEST sub

% % test_mat = all_mats(:,:,leftout);  
% test_mat = all_mats2(:,:,leftout); %resting validation 
%     test_sumpos = sum(sum(test_mat.*pos_mask))/2;
%     test_sumneg = sum(sum(test_mat.*neg_mask))/2;
%     
%     behav_pred_glm(leftout) = b(1)*test_sumpos + b(2)*test_sumneg + b(3);

    
end

%% plot_glm
FeatSel_no = FeatSel_no_pos+FeatSel_no_neg;
FeatSel_perc = mean(FeatSel_no) / ((no_node)*(no_node-1)/2) * 100;
fprintf('\n Average number of feature selected = %6.3f among  %6.3f',mean(FeatSel_no),no_node*(no_node-1)/2);
fprintf('\n Percentage of feature selected = %6.3f \n',FeatSel_perc);
[R_glm, P_glm] = corr(behav_pred_glm',all_behav)

strglm = ['r = ', num2str(R_glm),', P = ', num2str(P_glm), ' , Feat = ', num2str(mean(FeatSel_no_pos+FeatSel_no_neg))];

figure(1); plot(behav_pred_glm,all_behav,'r.');lsline; xlabel('Predicted behavior');ylabel('Observed behavior')
text(mean(behav_pred_glm)-std(behav_pred_glm)*.9,mean(all_behav)-std(behav_pred_glm)*.9,strglm,'HorizontalAlignment','left');
axis([min(all_behav),max(all_behav), min(all_behav),max(all_behav)]);


%% plot_pos/neg
FeatSel_no = FeatSel_no_pos+FeatSel_no_neg;
% FeatSel_perc = mean(FeatSel_no) / ((no_node)*(no_node-1)/2) * 100;
% fprintf('\n Average number of feature selected = %6.3f among  %6.3f',mean(FeatSel_no),no_node*(no_node-1)/2);
% fprintf('\n Percentage of feature selected = %6.3f',FeatSel_perc);
[R_pos, P_pos] = corr(behav_pred_pos,all_behav)
[R_neg, P_neg] = corr(behav_pred_neg, all_behav)
strpos = ['r = ', num2str(R_pos),', P = ', num2str(P_pos), ' , Feat = ', num2str(mean(FeatSel_no_pos))];
strneg = ['r = ', num2str(R_neg),', P = ', num2str(P_neg), ' , Feat = ', num2str(mean(FeatSel_no_neg))];
figure(1); plot(behav_pred_pos,all_behav,'r.');lsline; xlabel('Predicted behavior');ylabel('Observed behavior')
text(mean(behav_pred_pos)-std(behav_pred_pos)*.9,mean(all_behav)-std(behav_pred_pos)*.9,strpos,'HorizontalAlignment','left');
axis([min(all_behav),max(all_behav), min(all_behav),max(all_behav)]);
figure(2); plot(behav_pred_neg,all_behav,'b.');lsline; xlabel('Predicted behavior');ylabel('Observed behavior')
text(mean(behav_pred_neg)-std(behav_pred_neg)*.9,mean(all_behav)-std(behav_pred_neg)*.9,strneg,'HorizontalAlignment','left');    
axis([min(all_behav),max(all_behav), min(all_behav),max(all_behav)]);

%% find overlapping edges
overlap_mask_pos = find(overlap_pos > 57);
overlap_mask_neg = find(overlap_neg > 57);
overlap_loo_pos = zeros(227,25);
overlap_loo_neg = zeros(227,25);
overlap_loo_pos(overlap_mask_pos) = 1;
overlap_loo_neg(overlap_mask_neg) = 1;

save('overlap_pos_100per.txt','overlap_loo_pos','-ascii');
save('overlap_neg_100per.txt','overlap_loo_neg','-ascii');
%% compare predicted and observed scores using mean squared error
MSE_pos = sum((behav_pred_pos-all_behav).^2)/(no_sub-length(fit_pos)-1)
MSE_neg = sum((behav_pred_neg-all_behav).^2)/(no_sub-length(fit_neg)-1)
%% header

data = rand(4,3);
header = {'Col 1','Col 2','Col 3'};
output = [header; num2cell(data)]

sample = rand(3,3);
rowNames = {'a','b','c'};
colNames = {'x','y','z'};
sTable = array2table(sample,'RowNames',rowNames,'VariableNames',colNames)

% [R_pos, P_pos] = corr(behav_pred_pos,all_behav)
% [R_neg, P_neg] = corr(behav_pred_neg, all_behav)
% % strpos = ['r = ', num2str(R_pos),', P = ', num2str(P_pos), ' , Feat = ', num2str(mean(FeatSel_no_pos))];
% % strneg = ['r = ', num2str(R_neg),', P = ', num2str(P_neg), ' , Feat = ', num2str(mean(FeatSel_no_neg))];
% figure(1); plot(behav_pred_pos,all_behav,'r.');lsline; xlabel('Predicted behavior');ylabel('Observed behavior')
% axis([10,45, 10,45]);
% figure(2); plot(behav_pred_neg,all_behav,'b.','linewidth',25);lsline; xlabel('Predicted PU score');ylabel('Self-rated PU score')
% axis([10,45, 10,45]);


%return to 227*227 
% invsocog = [4 6	7	8	9	11	12	13	14	15	17	18	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	40	41	42	43	44	45	46	47			50	51	52	53	54		56	57	58	59	60	61	62	63			66	67	68	69	70	71	72	73	74	75	76	77	78	79	80	81	82	83	84	85	86	87	88	89		91	92	93	94	95	96	97	98	99	100	101	102	103	104	105	106	107	108	109	110	111	112	113	114	115	116			119		121	122	123	124		126	127		129	130	131	132	133	134	135	136	137	138	139	140	141	142	143	144	145	146	147	148	149	150	151	152	153	154	155	156	157	158	159	160	161	162		164		166		168	169		171	172	173	174	175	176	177	178	179	180	181	182	183	184	185	186	187	188	189	190	191	192	193	194	195	196	197	198	199	200	201		203	204		206	207	208	209	210		212	213	214	215	216	217	218	219	220	221	222	223	224	225	226	227
% ]

a = zeros(227,227)
a(:, [1 2 3 5 10 16 19	48	49	55	64	65	90	117	118	120	125	128	163	165	167	170	202	205	211]) = overlap_loo_neg(:,[1:25])
