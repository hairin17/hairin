%copy file to another directory
%make sure current directory is the directory you want to copy to
%this code skips folders(sub_datapath) that do not exist in datapath

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


        
        
