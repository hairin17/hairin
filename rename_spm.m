% 
% for k = 1
% files = dir(['F:\kshap2018\' num2str(k) '\cb\*.nii'])';
% for id=1:length(files)
% %     f = fileparts(files(id).name);
%     
%     movefile(files(id).name, sprintf('%03d.nii', id));
% end
% end

for k = 59:60
todir= ['E:\kshap2018_dicom\' num2str(k) '\field_map_0002'];   %cb, msit, t1, rs1, rs2, t2fl
cd(todir);

files = dir(['E:\kshap2018_dicom\' num2str(k) '\field_map_0002\*.nii'])';  %cb, msit, t1, rs1, rs2, t2fl
for id=1:length(files)
    
    movefile(files(id).name, sprintf('00-%03d.nii', id));
end
end