
all_mats = all_mats;
all_behav = table2array(all_demo2(:,52));
    
no_sub=size(all_mats,3);

% calculate true prediction correlation
[true_prediction_r_pos, true_prediction_r_neg] = pred_glm(all_mats,all_behav,1);  % (mat,behav,condition,thresh)

% number of iteration for permutation testing
no_iterations=10;
prediction_r=zeros(no_iterations,2);
prediction_r(1,1)=true_prediction_r_pos;
prediction_r(1,2)=true_prediction_r_neg;

% create estimate distribution of the test statistic
% via random shuffle of data labelsspm fmri
for it=2:no_iterations
    fprintf('\n Performing iteration %d out of %d', it, no_iterations);
    new_behav=all_behav(randperm(no_sub));
    [prediction_r(it,1),prediction_r(it,2)] = predict_behavior(all_mats,new_behav);
end

sorted_prediction_r_pos = sort(prediction_r(:,1),'descend');
position_pos = find(sorted_prediction_r_pos==true_prediction_r_pos);
pval_pos = position_pos(1)/no_iterations;
sorted_prediction_r_neg = sort(prediction_r(:,2),'descend');
position_neg = find(sorted_prediction_r_neg==true_prediction_r_neg);
pval_neg = position_neg(1)/no_iterations;


figure(1);hist(sorted_prediction_r_pos);title('Positive CPM Permuation Result');
hold on
plot([true_prediction_r_pos,true_prediction_r_pos],ylim,'r--','LineWidth',2)
hold off
figure(2); hist(sorted_prediction_r_neg);title('Negative CPM Permuation Result');
hold on
plot([true_prediction_r_neg,true_prediction_r_neg],ylim,'r--','LineWidth',2)
hold off