
%% head_motion_check

%% rp_abs_count
data_path = 'E:\04_kshap2018\!analysis\cb_conn18a_volbased_vdm\cb_volbased_vdm'
dirfiles1 = dir('E:\04_kshap2018\!analysis\cb_conn18a_volbased_vdm\cb_volbased_vdm'); 
dirfiles = dirfiles1(3:length(dirfiles1)); %excl . ..

for id  = 1:117

idfile = spm_select('FPList', fullfile(data_path, dirfiles(id).name, 'cb'), '^rp.*\.txt$'); % 불러들일 id
rp = load(idfile);
ncount(id,:) = sum(abs(rp) >= 3);
mov_id(id) = any(ncount(id,:)) > 0;

end

%% rp_dspl_mean
for id = [1:80, 82:117];

idfile = spm_select('FPList', fullfile(data_path, dirfiles(id).name, 'cb'), '^rp.*\.txt$'); % 불러들일 id
rp = load(idfile);

    for k = 1:179
            rp_subt(k,:) = rp(k+1,:) - rp(k,:);
    end

       rp_subt = mean(abs(rp_subt));
       rp_sub_id(id,:) = rp_subt
       %rp_displacement(id,:) = rp_subt;
end

%% rp_dspl_2mm prior
for id = 1:length(dirfiles)

idfile = spm_select('FPList', fullfile(data_path, dirfiles(id).name, 'cb'), '^rp.*\.txt$'); % 불러들일 id
rp = load(idfile);

    for k = 1:179
            rp_subt(k,:) = rp(k+1,:) - rp(k,:);
    end

      rp_subt = abs(rp_subt);

    for j = 1:179
        count(j) = any(rp_subt(j,:) > 2);
    end

rp_displace_2mm(id) = sum(count);
end

