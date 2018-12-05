-
% left out  for loop 바로 앞에:
overlap_pos = zeros(no_node,no_node);
overlap_neg = zeros(no_node,no_node);

-
 pos_mask(pos_edges)=1;
 neg_mask(neg_edges)=1;
% 여기 다음에 :
        
    % set overlapping masks
    overlap_pos = overlap_pos + pos_mask;
    overlap_neg = overlap_neg + neg_mask;

-
% for loop 끝나고
overlap_mask_pos = find(overlap_pos > 62);
overlap_mask_neg = find(overlap_neg > 62);
overlap_loo_pos = zeros(no_node,no_node);
overlap_loo_neg = zeros(no_node,no_node);
overlap_loo_pos(overlap_mask_pos) = 1;
overlap_loo_neg(overlap_mask_neg) = 1;