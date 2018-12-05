% Input: Conn FC matrix
% Output: extract lower triangle vector
% cd 'E:\04_kshap2018\!analysis\stats_cb_95\conn_cb95\shen_mat';
% 
% load('all_mats_e.mat')
% all_mats = all_mats_e;

function [all_vec] = lowervec(all_mats)
no_sub = size(all_mats,3);
no_node = size(all_mats,1);
submatmask = tril(true(size(all_mats(:,:,1))),-1);  % first subject to lower tri

lowervec = zeros(no_sub,no_node*(no_node-1)/2);
  for s=1:no_sub
        submat = all_mats(:,:,s); 
        subvec = submat(submatmask);
        lowervec(s,:) = subvec;
  end
  all_vec = lowervec;
end