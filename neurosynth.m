%% Neurosynth map
P= spm_select(Inf,'image');
V = spm_vol(P);
mentalizing = spm_read_vols(V);

P= spm_select(Inf,'image');
V = spm_vol(P);
ment = spm_read_vols(V);

% atlas % neurosynth mask와 차원이 같은지 확인*
P= spm_select(Inf,'image');
V = spm_vol(P);
shen = spm_read_vols(V);

% 
% socog_shenmask = shen(k10_bin);
% shen_mask_rois = unique(socog_shenmask);


% k10 mask
mentalizing_k10 = find(mentalizing > 0);
ment_k10 = find(ment > 0);

% shen overlap
shen_mentalizing_mask = shen(mentalizing_k10);
shen_mentmask = shen(ment_k10);

% cluster size
shen_mentalizing_tab = tabulate(shen_mentalizing_mask);
shen_ment_tab = tabulate(shen_mentmask);


% aal label
aal = spm_read_vols(V);