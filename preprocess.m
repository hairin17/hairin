%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging
% Edited by Hairin 2019

%2019 Yangsa longitudinal w2 data 
%imaging protocol (default)
%field map/MSIT*2/RS*2/T1/T2FLAIR/TOM*2

%%preprocess pipeline
%fieldmap/slice timing
%correction/realign/segmentation/imcalc/normalize/smoothing8mm/art using
%conn2018

%select data
clear all;
dataPath = '/Volumes/CNS_3/kshap2019_nifti'%temp
list = dir(dataPath)
list=list(~ismember({list.name},{'.','..'}));

subList = size(list)
subList = subList(:,1)

%rename
cd(dataPath);
for id=1:2%length(subList)
    filename = list(id).name;
    underlineLocation = strfind(filename,'_');
    filename = filename(1:underlineLocation - 1);
    movefile(list(id).name, sprintf(filename, id));

end

%

for i = [12:27]

    mod = 'TOM_1';
    
    subId = num2str(i, '%03d');
    toDir = fullfile(dataPath,subId,mod);
    cd(toDir);

% Initialise SPM 
%--------------------------------------------------------------------------
spm('Defaults','fMRI');
spm_jobman('initcfg');
spm_get_defaults('cmdline',true);


%select data with surfix
f = spm_select('FPList', fullfile(dataPath,subId,mod), '^f.*\.nii$'); %epi
a = spm_select('FPList', fullfile(dataPath,subId,'T1'), '^s.*\.nii$'); %structure

%m = spm_select('FPList', fullfile(dataPath, 'field_map_0002'), '^s00-001.*\.nii$');  %shorter TE -> renmame needed
%s = spm_select('FPList', fullfile(dataPath, 'field_map_0003'), '^s.*\.nii$');     %subtract file
%v = spm_select('FPList', fullfile(dataPath, 'field_map_0003'), '^vdm5_.*\.nii$');


% Slice timing Correction
%--------------------------------------------------------------------------
%KSHAP2019_yangsa_w2
%interleaved bottom up
%MSIT(TR/TE/TA);REST(TR/TE/TA);TOM(TR/TE/TA);
matlabbatch{1}.spm.temporal.st.scans = {cellstr(spm_file(f,'ext','nii'))};
matlabbatch{1}.spm.temporal.st.nslices = 38;
matlabbatch{1}.spm.temporal.st.tr = 2;
matlabbatch{1}.spm.temporal.st.ta = 2-(2/38); %TR-(TR/nSlices)
matlabbatch{1}.spm.temporal.st.so = [1:2:38 2:2:38];
matlabbatch{1}.spm.temporal.st.refslice = 1;
matlabbatch{1}.spm.temporal.st.prefix = 'a';

spm_jobman('run',matlabbatch);
clear matlabbatch

% Realign & unwarp /no vdm ver
%--------------------------------------------------------------------------

af = spm_select('FPList', fullfile(dataPath,subId,mod), '^af.*\.nii$'); %epi

matlabbatch{1}.spm.spatial.realignunwarp.data.scans = cellstr(af);
matlabbatch{1}.spm.spatial.realignunwarp.data.pmscan = ''%cellstr(spm_file(v,'ext','nii'));
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
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/Applications/matlab_tool/spm12/tpm/TPM.nii,1'};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/Applications/matlab_tool/spm12/tpm/TPM.nii,2'};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/Applications/matlab_tool/spm12/tpm/TPM.nii,3'};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'/Applications/matlab_tool/spm12/tpm/TPM.nii,4'};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'/Applications/matlab_tool/spm12/tpm/TPM.nii,5'};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'/Applications/matlab_tool/spm12/tpm/TPM.nii,6'};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'eastern';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 1];


spm_jobman('run',matlabbatch);
clear matlabbatch

%image calc
%--------------------------------------------------------------------------
cs = spm_select('FPList', fullfile(dataPath,subId,'T1'));
cd(fullfile(dataPath,subId,'T1'));
matlabbatch{1}.spm.util.imcalc.input = cellstr(cs([1:3 6],:))
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
b = spm_select('FPList', fullfile(dataPath,subId,'T1'), '^b.*\.nii$'); %brain
uaf = spm_select('FPList',fullfile(dataPath,subId,mod), '^uaf.*\.nii$'); 
m = spm_select('FPList',fullfile(dataPath,subId,mod), '^mean.*\.nii$'); 

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
uaf = spm_select('FPList',fullfile(dataPath,subId,mod), '^uaf.*\.nii$'); 
df =  spm_select('FPList', fullfile(dataPath,subId,'T1'), '^y.*\.nii$');

matlabbatch{1}.spm.spatial.normalise.write.subj.def      = cellstr(df);
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(uaf);
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                             78 76 85]; %bounding box f 
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

spm_jobman('run',matlabbatch);
clear matlabbatch


% Smooth
s = spm_select('FPList',fullfile(dataPath,subId,mod),'^w.*\.nii$');
matlabbatch{1}.spm.spatial.smooth.data = cellstr(s); %cellstr(spm_file(f,'prefix','w'));
matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
spm_jobman('run',matlabbatch);
clear matlabbatch


end

