%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Edited by Hairin 2018
%directory set up new 
data_path = 'E:\04_kshap2018\!analysis\stats_cb_95\'
dirfiles1 = dir('E:\04_kshap2018\!analysis\stats_cb_95\'); 
dirfiles = dirfiles1(3:length(dirfiles1)); %excl . ..

for id  = 1:length(dirfiles)

fid = id; 
idfile = spm_select('FPList', fullfile(data_path, dirfiles(fid).name, 'cb'), '^uaf.*\.nii$'); % 불러들일 id

BATCH.Setup.functionals{id} = idfile; % 넣고 싶은 자리 id

end


%% T1, field map rename 
for k = [104:126]
todir= ['E:\04_kshap2018\fulldata\' num2str(k,'%03d') '\field_map_0003'];   %t1, field_map_0002/field_map_0003
cd(todir);

files = dir(['E:\04_kshap2018\fulldata\' num2str(k,'%03d') '\field_map_0003\*.nii'])';  %t1, field_map_0002
for id=1:length(files)
    
    movefile(files(id).name, sprintf('s00-%03d.nii', id));
end
end

%% preprocess
% Directory containing the data
%--------------------------------------------------------------------------
for i = [1:9,11:31,33:57,59:66,68:80,82:85,87:99,101:102,104:106,108:117,119:126] %116 [1:9,11:31,33:57,59:66,68:80,82:85,87:99,101:102,104:106,108:117,119:126]
    todir = ['E:\04_kshap2018\!analysis\cb_spm_prep_vdm\' num2str(i, '%03d')]
        cd(todir);

%directory containing data
% data_path = fileparts(mfilename('E:\04_kshap2018\fulldata\')); %data_path에 current directory = todir 들어감
% if isempty(data_path), data_path = pwd; end %initialize해야지 Data_path들어감
%add data path moddata
data_path = ['E:\04_kshap2018\!analysis\cb_spm_prep_vdm\' num2str(i, '%03d')]
data_path2v = ['E:\04_kshap2018\moddata\fmap_cb\' num2str(i, '%03d') '\']
data_path2a = ['E:\04_kshap2018\moddata\t1\' num2str(i, '%03d') '\']
data_path2f = ['E:\04_kshap2018\moddata\cb\' num2str(i, '%03d') '\']

%%
% for sess = 1:2;
% for id = 1:length(dirfiles)
%     
%     % Specify 1)protocol name    2) nii prefix
%     idfile = spm_select('FPList', fullfile(data_path, dirfiles(id).name, ['rs' num2str(sess)]), '^f.*\.nii$');
%     
%     BATCH.Setup.functionals{id}{sess} = idfile;
% 
% end
% end




% Initialise SPM 
%--------------------------------------------------------------------------
spm('Defaults','fMRI');
spm_jobman('initcfg');
% spm_get_defaults('cmdline',true);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPATIAL PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%datapath = fulldata
%datapath2 = moddata
f = spm_select('FPList', fullfile(data_path, 'cb'), '^u.*\.nii$');
a = spm_select('FPList', fullfile(data_path,'t1'), '^s.*\.nii$');
%m = spm_select('FPList', fullfile(data_path, 'field_map_0002'), '^s00-001.*\.nii$');  %shorter TE -> renmame needed
%s = spm_select('FPList', fullfile(data_path, 'field_map_0003'), '^s.*\.nii$');     %subtract file
v = spm_select('FPList', fullfile(data_path2v, 'field_map_0003'), '^vdm5_.*\.nii$');

clear matlabbatch


% Field_map (done)
%--------------------------------------------------------------------------
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = cellstr(spm_file(s,'ext','nii'));
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = cellstr(spm_file(m,'ext','nii'));
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsfile = {'D:\matlab_code\pm_defaults_CB.m'};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session.epi = cellstr(spm_file(f,'ext','nii'));
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'cb';
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;

spm_jobman('run',matlabbatch);
clear matlabbatch

% Slice timing Correction
%--------------------------------------------------------------------------
matlabbatch{1}.spm.temporal.st.scans = {cellstr(spm_file(f,'ext','nii'))};
matlabbatch{1}.spm.temporal.st.nslices = 30;
matlabbatch{1}.spm.temporal.st.tr = 2;
matlabbatch{1}.spm.temporal.st.ta = 1.93333333333333;
matlabbatch{1}.spm.temporal.st.so = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30];
matlabbatch{1}.spm.temporal.st.refslice = 1;
matlabbatch{1}.spm.temporal.st.prefix = 'a';

spm_jobman('run',matlabbatch);
clear matlabbatch

% Realign & unwarp
%--------------------------------------------------------------------------
% matlabbatch{1}.spm.spatial.realign.estwrite.data = {cellstr(f)};
% matlabbatch{1}.spm.spatial.realignunwarp.data.pmscan = cellstr(spm_file(v,'ext','nii'));
% matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1];
af = spm_select('FPList', fullfile(data_path2f, 'cb'), '^af.*\.nii$');

matlabbatch{1}.spm.spatial.realignunwarp.data.scans = cellstr(af);
matlabbatch{1}.spm.spatial.realignunwarp.data.pmscan = cellstr(spm_file(v,'ext','nii'));
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 0;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 5; %2 to 5
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';

spm_jobman('run',matlabbatch);
clear matlabbatch
% Segment
%--------------------------------------------------------------------------
matlabbatch{1}.spm.spatial.preproc.channel.vols = cellstr(a);
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'C:\Program Files\MATLAB\R2017b\toolbox\spm12\tpm\TPM.nii,1'};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [1 0];%wc1
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'C:\Program Files\MATLAB\R2017b\toolbox\spm12\tpm\TPM.nii,2'};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [1 0];%wc2
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'C:\Program Files\MATLAB\R2017b\toolbox\spm12\tpm\TPM.nii,3'};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [1 0];%wc3
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'C:\Program Files\MATLAB\R2017b\toolbox\spm12\tpm\TPM.nii,4'};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'C:\Program Files\MATLAB\R2017b\toolbox\spm12\tpm\TPM.nii,5'};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'C:\Program Files\MATLAB\R2017b\toolbox\spm12\tpm\TPM.nii,6'};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 1];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'eastern';

spm_jobman('run',matlabbatch);
clear matlabbatch

%image calc
%--------------------------------------------------------------------------
aa = spm_select('FPList', fullfile(data_path,'t1'));
 
% todir = ['E:\kshap2018_dicom\' num2str(i) '\t1\'];
%         cd(todir);
matlabbatch{1}.spm.util.imcalc.input = cellstr(aa([1:3 6],:))
matlabbatch{1}.spm.util.imcalc.output = 'brain';
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2+i3).*i4';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

spm_jobman('run',matlabbatch);
clear matlabbatch

% Coregister
%--------------------------------------------------------------------------
b = spm_select('FPList', fullfile(data_path), '^b.*\.nii$');
uaf = spm_select('FPList',fullfile(data_path,'cb'), '^uaf.*\.nii$'); 
m = spm_select('FPList',fullfile(data_path,'cb'), '^mean.*\.nii$'); 

matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr(b);
matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr(m);
matlabbatch{1}.spm.spatial.coreg.estimate.other = cellstr(uaf);
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
spm_jobman('run',matlabbatch);
clear matlabbatch

% Normalise: Write f
%--------------------------------------------------------------------------
%
uaf = spm_select('FPList',fullfile(data_path,'cb'), '^uaf.*\.nii$'); 
df =  spm_select('FPList', fullfile(data_path,'t1'), '^y.*\.nii$');

matlabbatch{1}.spm.spatial.normalise.write.subj.def      = cellstr(df);
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(uaf);
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                             78 76 85]; %bounding box f 
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

spm_jobman('run',matlabbatch);
clear matlabbatch

% Normalise: Write structure %%bounding box option change
st =  spm_select('FPList', fullfile(data_path), '^b.*\.nii$');

matlabbatch{1}.spm.spatial.normalise.write.subj.def      = cellstr(df);
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(st);
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
                                                             78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

spm_jobman('run',matlabbatch);
clear matlabbatch

% Smooth
s = spm_select('FPList',fullfile(data_path,'cb'),'^w.*\.nii$');
matlabbatch{1}.spm.spatial.smooth.data = cellstr(s); %cellstr(spm_file(f,'prefix','w'));
matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
spm_jobman('run',matlabbatch);
clear matlabbatch


end

%% first level modeling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLM SPECIFICATION, ESTIMATION, INFERENCE, RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = [7	9	11	12	13	14	17	18	19	20	21	22	23	24	25	27	28	29	30	34	35	36	37	38	39	41	42	45	46	47	48	49	50	51	52	53	54	55	56	57	59	60	61	62	63	64	65	68	69	70	71	73	76	77	78	79	80	82	83	84	87	88	89	90	91	92	94	95	96	97	98	99	101	102	104	105	106	108	109	110	111	112	113	114	115	116	117	119	120	121	122	123	124	125	126]
%k = cb_sub95
    todir = ['E:\04_kshap2018\!analysis\cb_spm_prep_vdm\' num2str(k,'%03d') '\'];
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
f = spm_select('FPList', fullfile(data_path,'cb'), '^sw.*\.nii$');
rp = spm_select('FPList', fullfile(data_path,'cb'), '^rp.*\.txt$'); %fullfile(data_path,'cb'), '^sw.
clear matlabbatch

% Output Directory
%--------------------------------------------------------------------------
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = cellstr(data_path);
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'GLM_cb';

% Model Specification
%--------------------------------------------------------------------------
matlabbatch{2}.spm.stats.fmri_spec.dir = cellstr(fullfile(data_path,'GLM_cb'));
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
matlabbatch{3}.spm.stats.fmri_est.spmmat = cellstr(fullfile(data_path,'GLM_cb','SPM.mat'));

% Contrasts %E-I, RI-I, E-RI 
%--------------------------------------------------------------------------
matlabbatch{4}.spm.stats.con.spmmat = cellstr(fullfile(data_path,'GLM_cb','SPM.mat'));
matlabbatch{4}.spm.stats.con.consess{1}.tcon.name = 'Exclusion > Inclusion';
matlabbatch{4}.spm.stats.con.consess{1}.tcon.weights = [-1 1 0];

matlabbatch{4}.spm.stats.con.consess{2}.tcon.name = 'Re-inclusion > inclusion';
matlabbatch{4}.spm.stats.con.consess{2}.tcon.weights = [-1 0 1];

matlabbatch{4}.spm.stats.con.consess{3}.tcon.name = 'Inclusion > Exclusion';
matlabbatch{4}.spm.stats.con.consess{3}.tcon.weights = [1 -1 0];

% Inference Results
%--------------------------------------------------------------------------
matlabbatch{5}.spm.stats.results.spmmat = cellstr(fullfile(data_path,'GLM_cb','SPM.mat'));
matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{5}.spm.stats.results.conspec.extent = 0;
matlabbatch{5}.spm.stats.results.print = false;

% Rendering
%--------------------------------------------------------------------------
matlabbatch{6}.spm.util.render.display.rendfile = {fullfile(spm('Dir'),'canonical','cortex_20484.surf.gii')};
matlabbatch{6}.spm.util.render.display.conspec.spmmat = cellstr(fullfile(data_path,'GLM_cb','SPM.mat'));
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