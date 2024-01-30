% rename the LFP files

folder = '/Users/zli328/Documents/Research/Research/Project_Hippocampus/raw_data/403/NLX files/403-045_TMS/2017-10-23_14-54-26/LFP';
filelist = dir(folder);

for i = 3:numel(filelist)
    
    oldName = filelist(i).name;
    oldpath = fullfile(folder,oldName);
    
    %     newName = filelist(i).name(1:end-9);
    newName = [filelist(i).name,'.ncs'];
    newpath = fullfile(folder,newName);
    
    movefile(oldpath,newpath)
end


