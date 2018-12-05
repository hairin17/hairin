%% first level modeling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLM SPECIFICATION, ESTIMATION, INFERENCE, RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = [7	9	11	12	13	14	17	18	19	20	21	22	23	24	25	27	28	29	30	34	35	36	37	38	39	41	42	45	46	47	48	49	50	51	52	53	54	55	56	57	59	60	61	62	63	64	65	68	69	70	71	73	76	77	78	79	80	82	83	84	87	88	89	90	91	92	94	95	96	97	98	99	101	102	104	105	106	108	109	110	111	112	113	114	115	116	117	119	120	121	122	123	124	125	126]
%k = cb_sub95
    todir = ['E:\04_kshap2018\!analysis\cb_spm_prep_vdm\' num2str(k,'%03d') '\cb'];
        cd(todir);

% data_path = ['E:\04_kshap2018\!analysis\cb_conn18a_volbased_vdm\' num2str(k,'%03d') '\cb']

%directory containing data
data_path = fileparts(mfilename(['E:\04_kshap2018\!analysis\cb_spm_prep_vdm\' num2str(k,'%03d') '\'])); %data_path에 current directory 들어감
if isempty(data_path), data_path = pwd; end %initialize해야지 Data_path들어감

% Initialise SPM 
%--------------------------------------------------------------------------
spm('Defaults','fMRI');
spm_jobman('initcfg');

% spm_get_defaults('cmdline',true);


%modeling
f = spm_select('FPList', fullfile(data_path), '^s.*\.nii$'); %1:9
rp = spm_select('FPList', fullfile(data_path), '^rp.*\.txt$');
fileID = fopen(['bt_' num2str(k) '.txt'],'r');
formatSpec = '%f';
bt = fscanf(fileID, formatSpec);

clear matlabbatch

% Output Directory
%--------------------------------------------------------------------------
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = cellstr(data_path);
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'GLM_cb_bt';

% Model Specification
%--------------------------------------------------------------------------
matlabbatch{2}.spm.stats.fmri_spec.dir = cellstr(fullfile(data_path,'GLM_cb_bt'));
matlabbatch{2}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{2}.spm.stats.fmri_spec.timing.RT = 2;
matlabbatch{2}.spm.stats.fmri_spec.sess.scans = cellstr(f);
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).name = 'incl';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).onset = 11;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).duration = 115;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).name = 'excl';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).onset = 126;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).duration = 115;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).orth = 1;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).name = 'reincl';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).onset = 241;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).duration = 115;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).tmod = 0;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).orth = 1;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(4).name = 'fix';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(4).onset = 0;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(4).duration = 11;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(4).tmod = 0;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(4).orth = 1;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(5).name = 'fix2';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(5).onset = 356;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(5).duration = 4;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(5).tmod = 0;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(5).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(5).orth = 1;
%button response event
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(6).name = 'bt';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(6).onset = bt;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(6).duration = 0;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(6).tmod = 0;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(6).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(6).orth = 1;

matlabbatch{2}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{2}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{2}.spm.stats.fmri_spec.sess.multi_reg = cellstr(rp); %subject number
matlabbatch{2}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{2}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{2}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{2}.spm.stats.fmri_spec.volt = 1;
matlabbatch{2}.spm.stats.fmri_spec.global = 'None';
matlabbatch{2}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{2}.spm.stats.fmri_spec.mask = {''};
matlabbatch{2}.spm.stats.fmri_spec.cvi = 'AR(1)';

% Model Estimation
%--------------------------------------------------------------------------
matlabbatch{3}.spm.stats.fmri_est.spmmat = cellstr(fullfile(data_path,'GLM_cb_bt','SPM.mat'));

% Contrasts %E-I, RI-I, E-RI 
%--------------------------------------------------------------------------
matlabbatch{4}.spm.stats.con.spmmat = cellstr(fullfile(data_path,'GLM_cb_bt','SPM.mat'));
matlabbatch{4}.spm.stats.con.consess{1}.tcon.name = 'Exclusion > Inclusion';
matlabbatch{4}.spm.stats.con.consess{1}.tcon.weights = [-1 1 0];

matlabbatch{4}.spm.stats.con.consess{2}.tcon.name = 'Re-inclusion > inclusion';
matlabbatch{4}.spm.stats.con.consess{2}.tcon.weights = [-1 0 1];

matlabbatch{4}.spm.stats.con.consess{3}.tcon.name = 'Inclusion > Exclusion';
matlabbatch{4}.spm.stats.con.consess{3}.tcon.weights = [1 -1 0];

matlabbatch{4}.spm.stats.con.consess{4}.tcon.name = 'Inclusion';
matlabbatch{4}.spm.stats.con.consess{4}.tcon.weights = [1 0 0];

matlabbatch{4}.spm.stats.con.consess{5}.tcon.name = 'Exclusion';
matlabbatch{4}.spm.stats.con.consess{5}.tcon.weights = [0 1 0];

matlabbatch{4}.spm.stats.con.consess{6}.tcon.name = 'Reinclusion';
matlabbatch{4}.spm.stats.con.consess{6}.tcon.weights = [0 0 1];
% Inference Results
%--------------------------------------------------------------------------
matlabbatch{5}.spm.stats.results.spmmat = cellstr(fullfile(data_path,'GLM_cb_bt','SPM.mat'));
matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{5}.spm.stats.results.conspec.extent = 0;
matlabbatch{5}.spm.stats.results.print = false;

% Rendering
%--------------------------------------------------------------------------
matlabbatch{6}.spm.util.render.display.rendfile = {fullfile(spm('Dir'),'canonical','cortex_20484.surf.gii')};
matlabbatch{6}.spm.util.render.display.conspec.spmmat = cellstr(fullfile(data_path,'GLM_cb_bt','SPM.mat'));
matlabbatch{6}.spm.util.render.display.conspec.contrasts = 1;
matlabbatch{6}.spm.util.render.display.conspec.threshdesc = 'FWE';
matlabbatch{6}.spm.util.render.display.conspec.thresh = 0.05;
matlabbatch{6}.spm.util.render.display.conspec.extent = 0;

spm_jobman('run',matlabbatch);
end


%% modeling w/ button 
%% first level modeling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLM SPECIFICATION, ESTIMATION, INFERENCE, RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = [1:31, 33:57, 59:60]
%     todir = ['E:\kshap2018_dicom\' num2str(i) '\'];
%         cd(todir);

i = 1
todir = ['E:\kshap2018_conn\cb\s\'];
cd(todir);

%directory containing data
data_path = fileparts(mfilename('E:\kshap2018_dicom\')); %data_path에 current directory 들어감
if isempty(data_path), data_path = pwd; end %initialize해야지 Data_path들어감

% Initialise SPM 
%--------------------------------------------------------------------------
spm('Defaults','fMRI');
spm_jobman('initcfg');
% spm_get_defaults('cmdline',true);

%modeling

f = spm_select('FPList', fullfile(data_path,'cb'), '^s.*\.nii$'); %1:9
rp = spm_select('FPList', fullfile(data_path, 'cb'), '^rp.*\.txt$');
fileID = fopen('bt.txt','r');
formatSpec = '%f';
bt = fscanf(fileID, formatSpec);

clear matlabbatch

% Output Directory
%--------------------------------------------------------------------------
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = cellstr(data_path);
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'GLM_bt3';

% Model Specification
%--------------------------------------------------------------------------
matlabbatch{2}.spm.stats.fmri_spec.dir = cellstr(fullfile(data_path,'GLM_bt3'));
matlabbatch{2}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{2}.spm.stats.fmri_spec.timing.RT = 2;
matlabbatch{2}.spm.stats.fmri_spec.sess.scans = cellstr(f);
%inclusion
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).name = 'incl';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).onset = 11;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).duration = 115;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
%exclusion
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).name = 'excl';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).onset = 126;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).duration = 115;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(2).orth = 1;
%re-inclusion
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).name = 'reincl';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).onset = 241;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).duration = 115;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).tmod = 0;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(3).orth = 1;
%fixation-1
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(4).name = 'fix';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(4).onset = 0;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(4).duration = 11;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(4).tmod = 0;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(4).orth = 1;
%fixation-2
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(5).name = 'fix2';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(5).onset = 356;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(5).duration = 4;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(5).tmod = 0;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(5).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(5).orth = 1;
%button response event
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(6).name = 'bt';
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(6).onset = bt;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(6).duration = 0;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(6).tmod = 0;
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(6).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{2}.spm.stats.fmri_spec.sess.cond(6).orth = 1;
%multiple regressors 
matlabbatch{2}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{2}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{2}.spm.stats.fmri_spec.sess.multi_reg = cellstr(rp); %subject number
matlabbatch{2}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{2}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{2}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{2}.spm.stats.fmri_spec.volt = 1;
matlabbatch{2}.spm.stats.fmri_spec.global = 'None';
matlabbatch{2}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{2}.spm.stats.fmri_spec.mask = {''};
matlabbatch{2}.spm.stats.fmri_spec.cvi = 'AR(1)';

% Model Estimation
%--------------------------------------------------------------------------
matlabbatch{3}.spm.stats.fmri_est.spmmat = cellstr(fullfile(data_path,'GLM_bt3','SPM.mat'));

% Contrasts %E-I, RI-I, E-RI 
%--------------------------------------------------------------------------
matlabbatch{4}.spm.stats.con.spmmat = cellstr(fullfile(data_path,'GLM_bt3','SPM.mat'));
matlabbatch{4}.spm.stats.con.consess{1}.tcon.name = 'Exclusion > Inclusion';
matlabbatch{4}.spm.stats.con.consess{1}.tcon.weights = [-1 1 0];

matlabbatch{4}.spm.stats.con.consess{2}.tcon.name = 'Re-inclusion > inclusion';
matlabbatch{4}.spm.stats.con.consess{2}.tcon.weights = [-1 0 1];

matlabbatch{4}.spm.stats.con.consess{2}.tcon.name = 'Exclusion > Re-inclusion';
matlabbatch{4}.spm.stats.con.consess{2}.tcon.weights = [0 1 -1];

% Inference Results
%--------------------------------------------------------------------------
matlabbatch{5}.spm.stats.results.spmmat = cellstr(fullfile(data_path,'GLM_bt3','SPM.mat'));
matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{5}.spm.stats.results.conspec.extent = 0;
matlabbatch{5}.spm.stats.results.print = false;

% Rendering
%--------------------------------------------------------------------------
matlabbatch{6}.spm.util.render.display.rendfile = {fullfile(spm('Dir'),'canonical','cortex_20484.surf.gii')};
matlabbatch{6}.spm.util.render.display.conspec.spmmat = cellstr(fullfile(data_path,'GLM_bt3','SPM.mat'));
matlabbatch{6}.spm.util.render.display.conspec.contrasts = 1;
matlabbatch{6}.spm.util.render.display.conspec.threshdesc = 'FWE';
matlabbatch{6}.spm.util.render.display.conspec.thresh = 0.05;
matlabbatch{6}.spm.util.render.display.conspec.extent = 0;

spm_jobman('run',matlabbatch);
end