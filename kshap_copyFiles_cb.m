%copy file to another directory
%make sure current directory is the directory you want to copy to
%this code skips folders(sub_datapath) that do not exist in datapath
E:\04_kshap2018\moddata\fmap_cb\001\field_map_0003
% 001KCS	0
% 032SMS	0
% 058NSJ	0
% 067KJO	0
% 081JAJ	0
% 086KMH	0
% 100NSY	0
% 103KBS	0
% 107JNS	0
% 118JGS	0
% cb & t1 available 116 id 
for i= [1:9,11:31,33:57,59:66,68:80,82:85,87:99,101:102,104:106,108:117,119:126]
%    [status, message, messageid]=mkdir(num2str(i,'%03d'),'cb');
    % [status, message, messageid]=mkdir(num2str(i,'%03d'),'t1');
    datapath=['E:\04_kshap2018\!analysis\cb_conn18a_volbased_vdm\cb_volbased_vdm\' num2str(i,'%03d') '\cb']; %E:\04_kshap2018\!analysis\cb_conn18a_volbased_vdm\t1_cb_vol\001\t1
    copypath=['E:\04_kshap2018\!analysis\cb_spm_prep_vdm\' num2str(i,'%03d') '\cb']; %t1
    
    %copy by folder
    type='rp*.*'; %s
    sub_datapath=[datapath '\' type];
%     if exist(sub_datapath, 'file')~=0
%         sub_copypath=[copypath '\' type];
        copyfile(sub_datapath, copypath);
%     end
    %add more sets
%     type='cb';
%     sub_datapath=[datapath '\' type];
%     if exist(sub_datapath, 'file')~=0
%         sub_copypath=[copypath '\' type];
%         copyfile(sub_datapath, sub_copypath);
%     end
end

%% copy sfmri 95sub
for i = [ 7 9   11  12  13  14  17  18  19  20  21  22  23  24  25  27  28  29  30  34  35  36  37  38  39  41  42  45  46  47  48  49  50  51  52  53  54  55  56  57  59  60  61  62  63  64  65  68  69  70  71  73  76  77  78  79  80  82  83  84  87  88  89  90  91  92  94  95  96  97  98  99  101 102 104 105 106 108 109 110 111 112 113 114 115 116 117 119 120 121 122 123 124 125 126]

datapath=['E:\04_kshap2018\!analysis\cb_spm_prep_vdm\' num2str(i,'%03d') '\cb']; %E:\04_kshap2018\!analysis\cb_conn18a_volbased_vdm\t1_cb_vol\001\t1
copypath=['E:\04_kshap2018\!analysis\stats_cb_95\' num2str(i,'%03d') '\']; %t1
    
    %copy by folder
    type = 'rp*.*'; %s
    sub_datapath=[datapath '\' type];
%     if exist(sub_datapath, 'file')~=0
%         sub_copypath=[copypath '\' type];
%     end
        copyfile(sub_datapath, copypath);
end

        
        