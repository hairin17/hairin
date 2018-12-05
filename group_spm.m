%% Information about contrast
% con1 = 'Exclusion > Inclusion';
% con2 = 'Re-inclusion > inclusion';
% con3 = 'Inclusion > Exclusion';
% con4 = 'Inclusion';
% con5 = 'Exclusion';
% con6 = 'Reinclusion';

%% basic path setting
clear all;
files = dir(['E:\04_kshap2018\!analysis\stats_cb_75\'])';  %t1, field_map_0002
files = files(3:length(files));
data_path = 'E:\04_kshap2018\!analysis\stats_cb_75'

for i = 1:75
conList{i,1} = fullfile(data_path, files(i).name, 'con_0002.nii,1'); %
end

%%
c = load('cov.mat')
cname = 'sad'
dirSave = 'sum_stat'
pathSave = fullfile(data_path, dirSave, 'reincl-incl');

%%
matlabbatch{1}.spm.stats.factorial_design.dir = {'E:\04_kshap2018\!analysis\stats_cb_75\sum_stat\reincl-incl'}; %
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = conList;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov = struct('c', {}, 'cname', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;


spm_jobman('interactive', matlabbatch);
spm_jobman('run',matlabbatch);
 
clear matlabbatch;
% 
% 
% matlabbatch{1}.spm.stats.fmri_est.spmmat = {'E:\04_kshap2018\!analysis\stats_cb_95\sum_stats\glm_sad\SPM.mat'};
% spm_jobman('run',matlabbatch);
% % 

