% % Inputs
% all_demo1=xlsread('C:\Users\Seyul\Dropbox\Researches\Researches\A_journal_Neuromarker\rsFC63_allsub.xlsx');
% all_demo2=readtable('C:\Users\Seyul\Dropbox\Researches\Researches\A_journal_Neuromarker\rsFC63_allsub.xlsx');
% all_mats_shen = all_mats_extract(268,63);  % all_mats_extract(no_nodes, no_subs) in "resultsROI_Condition001.mat" - Z
% 
% all_mats = all_mats_shen;
% all_mats(isnan(all_mats))=0;
% all_behav = table2array(all_demo2(:,23));  % 23 ~ 27: UPPSP
%     
% no_sub=size(all_mats,3);

% calculate true prediction correlation
[true_prediction_r_glm]=predict_behavior_glm(all_mats,all_behav,3);  % (mat,behav,condition,thresh)

% number of iteration for permutation testing
no_iterations=5000;
prediction_r=zeros(no_iterations,2);
prediction_r(1,1)=true_prediction_r_glm;

% create estimate distribution of the test statistic
% via random shuffle of data labelsspm fmri
for it=2:no_iterations
    fprintf('\n Performing iteration %d out of %d', it, no_iterations);
    new_behav=all_behav(randperm(no_sub));
%     [prediction_r(it,1),prediction_r(it,2)] = predict_behavior(all_mats,new_behav); 
    [prediction_r(it,1)] = predict_behavior_glm(all_mats,new_behav);
end

sorted_prediction_r_glm = sort(prediction_r(:,1),'descend');
position_glm = find(sorted_prediction_r_glm==true_prediction_r_glm);
pval_glm = position_glm(1)/no_iterations;



figure(1);hist(sorted_prediction_r_glm);title('Positive CPM Permuation Result');
hold on
plot([true_prediction_r_glm,true_prediction_r_glm],ylim,'r--','LineWidth',2)
hold off