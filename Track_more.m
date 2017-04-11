clc
clear

[FileName, PathName] = uigetfile(pwd,filesep,'*.mp4','MultiSelect','on');

t_start = tic;
for j = 1:length(FileName)
   v =[FileName{j}(1:end-4) '_tracked'];
   mouse = Track_Video(FileName{j},PathName);
   save(v,'mouse');
end
 toc(t_start)