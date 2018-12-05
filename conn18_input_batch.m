%run conn 
%fill conn set up options (n, TR, sessions)
%current folder = data_path
%% cb functionals
data_path = pwd
dirfiles1 = dir('E:\04_kshap2018\!analysis\stats_cb_95')
dirfiles = dirfiles1(3:97); %exclusion

for id  = 1:length(dirfiles)

% idfile = spm_select('FPList', fullfile(data_path, dirfiles(id).name, 'cb'), '^f.*\.nii$'); 
idfile = spm_select('FPList', fullfile(data_path, dirfiles(id).name), '^sw.*\.nii$'); 
BATCH.Setup.functionals{id} = idfile; 

end

conn_batch(BATCH);

%% cb structurals
data_path = pwd
dirfiles1 = dir('E:\04_kshap2018\!analysis\stats_cb_95')
dirfiles = dirfiles1(3:97); %exclusion

for id  = 1:length(dirfiles)

idfile = spm_select('FPList', fullfile(data_path, dirfiles(id).name), '^wbrain.*\.nii$'); 

BATCH.Setup.structurals{id} = idfile; 
end

conn_batch(BATCH);

%% cb mask
data_path = pwd
dirfiles1 = dir('E:\04_kshap2018\!analysis\stats_cb_95')
dirfiles = dirfiles1(3:97); %exclusion

for id  = 1:length(dirfiles)

idfile = spm_select('FPList', fullfile(data_path, dirfiles(id).name), '^wc3.*\.nii$'); 

BATCH.Setup.masks.CSF.files{id} = idfile; 
end

conn_batch(BATCH);

%% cb realignment
data_path = pwd
dirfiles1 = dir('E:\04_kshap2018\!analysis\stats_cb_95')
dirfiles = dirfiles1(3:97); %exclusion

for id  = 1:length(dirfiles)

idfile = spm_select('FPList', fullfile(data_path, dirfiles(id).name), '^rp.*\.txt$'); 
BATCH.Setup.covariates.names={'realignment'};
BATCH.Setup.covariates.files{1}{id}{1} = idfile; 
end
conn_batch(BATCH);


%% rs functional 2sessions
data_path = 'E:\04_kshap2018\!analysis\rs_conn18a_defaul_NoVdm\rs'
dirfiles1 = dir('E:\04_kshap2018\data\'); 
dirfiles = dirfiles1(3:length(dirfiles1)); %excl . ..
dirfiles = dirfiles([1:45,47:57,59:99,101:102,104:126])
%rs session1
for id  = 1:length(dirfiles)


idfile = spm_select('FPList', fullfile(data_path, dirfiles(id).name, 'rs1'), '^f.*\.nii$'); 
idfile2 = spm_select('FPList', fullfile(data_path, dirfiles(id).name, 'rs2'), '^f.*\.nii$'); 

BATCH.Setup.functionals{id}{1} = idfile; % 넣고 싶은 자리 id BATCH.Setup.functionals{sub#}{session}
BATCH.Setup.functionals{id}{2} = idfile2; % 넣고 싶은 자리 id BATCH.Setup.functionals{sub#}{session}

end

conn_batch(BATCH);
