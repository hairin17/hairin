%% mixed subject ANOVA for GLM 1
% using flexible factorial design 
% 
% 2018 PNU social neuroscience lab fMRI workshop.
% Minhee Yoo, 20180520
%%%%%%%%%%%%%%%%%

%% Information about contrast
% con1~3 : choice (self, fair, alt)
% con4~6 : reward - no reward (self, fair, alt)

%% 
clearvars;
clc;

%% basic path setting

tmp = fileparts(pwd);
pathBase    = fullfile(tmp,'kshap2018_dicom');
dirStat     = 'stat';
numModel    = 'GLM1';
dirSave     = 'sum_stat';
dirJobSave  = 'jobs';

% cov 
load(fullfile(pathBase, 'behavData', 'fMRIAnalysis_subjList.mat'));
subjList    = subjCondInfo(:,1); %col 1 = sub number
subjCond    = subjCondInfo(:,2); %col 2 = sub cond 


numCon = [1,3];
pathSave = fullfile(pathBase, dirSave, numModel, 'mANOVAFlex_choice');

%%

matlabbatch{1}.spm.stats.factorial_design.dir = {pathSave};
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'subject';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'group';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).name = 'target';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).ancova = 0;


for sIndex = 1:length(subjList)
    subjID = subjList(sIndex);
    condMatrix = [subjCond(sIndex) 1; subjCond(sIndex) 2]; %1 = proself, 2 = prosocial
    
    for i = 1:length(numCon)
        conList{i,1} = fullfile(pathBase, sprintf('subject%02d', subjID), dirStat, numModel, sprintf('con_%04d.nii,1', numCon(i)));
    end
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(sIndex).scans = conList;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(sIndex).conds = condMatrix;
end

matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.inter.fnums = [2;3];

matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) =...
    cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

spm_jobman('interactive', matlabbatch);
