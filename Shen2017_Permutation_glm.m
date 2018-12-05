
all_mats = all_mats;
all_behav = table2array(all_demo2(:,52));
cov = table2array(all_demo2(:,141));

no_sub=size(all_mats,3);

% calculate true prediction correlation
[true_prediction_r_glmcov] = pred_glm_cov(all_mats,all_behav);  % (mat,behav,condition,thresh) %%use pred_glmcov %%true_prediction_r_glm => true_prediction_r_glmcov

% number of iteration for permutation testing
no_iterations=10;
prediction_r=zeros(no_iterations,1);
prediction_r(1,1)= true_prediction_r_glmcov


% create estimate distribution of the test statistic
% via random shuffle of data labelsspm fmri
for it=2:no_iterations
    fprintf('\n Performing iteration %d out of %d', it, no_iterations);
    new_behav=all_behav(randperm(no_sub));
    [prediction_r(it,1)] = predict_behavior(all_mats,new_behav);
end

sorted_prediction_r_glm = sort(prediction_r(:,1),'descend');
position_glm = find(sorted_prediction_r_glm==true_prediction_r_glmcov);
pval_glm = position_glm(1)/no_iterations;

% sorted_prediction_r_neg = sort(prediction_r(:,2),'descend');
% position_neg = find(sorted_prediction_r_neg==true_prediction_r_neg);
% pval_neg = position_neg(1)/no_iterations;


figure(1);hist(sorted_prediction_r_glm);title('GLM CPM Permuation Result');
hold on
plot([true_prediction_r_glmcov,true_prediction_r_glmcov],ylim,'r--','LineWidth',2)
hold off
% figure(2); hist(sorted_prediction_r_neg);title('Negative CPM Permuation Result');
% hold on
% plot([true_prediction_r_neg,true_prediction_r_neg],ylim,'r--','LineWidth',2)
% hold off