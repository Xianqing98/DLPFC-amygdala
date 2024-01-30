close all
clear

%% Set up directories
dataDir =  'C:\Users\xjl19\Desktop\TMS-iEEG\xl_Data_Processed';
queueFName = 'C:\Users\xjl19\Desktop\TMS-iEEG\Processed_xlsx\xl_FixTrigs_Queue.xlsx';


%% ESTABLISH DIRECTORY AND ADD PATH
% cd(fileparts(which(mfilename))); %Navigate to the script's location
addpath 'C:\Users\xjl19\Desktop\Matlab\fieldtrip-20230206'
ft_defaults;

addpath(genpath('.\Jeff_Scripts'), ...
    genpath('.\dannyIEPScripts'), ...
    genpath('.\dannyGenScripts'), ...
    genpath('.\MatlabImportExport_v6.0.0'), ...
    genpath('.\util'));

sham = true;
pre_epoch = 1.0; % Length to integrate beforehand
length_epoch = 2.0; % Length of epoch in Seconds

%% Import file names to be processed
queue = readtable(queueFName);

for i = 1:height(queue)
    disp("-------------------------------------" + newline + ...
        "Working on Patient " + num2str(queue{i,'Patient'}) + ...
        ", File " + queue{i,'File_Name'}{1})
    
    if queue{i,'Fixed'} == 1
        disp('Patient already processed. Moving on.')
        continue;
    end
    tic;
    %Save Directory for everything
    dataDirSess = fullfile(dataDir, num2str(queue{i,'Patient'}), queue{i,'Folder_Name'}{1});
    
    %% IMPORT DATA
    disp('Loading File')
    fName = fullfile(dataDirSess,'3-Subsampled',strcat(queue{i,'File_Name'}{1},'_3A-CleanSubsampled.mat'));
    load(fName);
    
    %% REMOVE BAD TRIGGERS
    disp('Removing Triggers...')
    remove_starts = queue{i,'Remove_Starts'};
    remove_ends = queue{i,'Remove_Ends'};
    
    mask = zeros(size(trigTimes));
    
    % If only one window to remove, then will be numerical
    if isnumeric(remove_starts)
        if ~isnan(remove_starts)
            mask = (trigTimes >= remove_starts) & (trigTimes <= remove_ends);
            trigTimes(mask) = [];
            trigLengths(mask) = [];
        end
    else
        % Otherwise need to split by semicolons
        remove_starts = remove_starts{1};
        remove_ends = remove_ends{1};
        
        if ~isempty(remove_starts)
            remove_starts = str2double(split(remove_starts,';'));
            remove_ends = str2double(split(remove_ends,';'));
            for window = 1:length(remove_starts)
                start_window = remove_starts(window);
                end_window = remove_ends(window);
                mask = (trigTimes >= start_window) & (trigTimes <= end_window);
                trigTimes(mask) = [];
                trigLengths(mask) = [];
            end
        end
    end
    
    disp(strcat(num2str(sum(mask)), " Triggers Removed"))
    
    trigsChanged = false;
    
    if sum(mask)
        trigsChanged = true;
    end
    
    %% Save New Trigger Timings
    if trigsChanged
        disp('Saving Trigger Timings')
        save(fName, 'trigTimes', 'trigShift', '-append')
        
        fName = fullfile(dataDirSess,'0-Raw',strcat(queue{i,'File_Name'}{1},'_0-Raw.mat'));
        save(fName, 'trigTimes', 'trigShift', '-append')
        
        fName = fullfile(dataDirSess,'1-Cleaned',strcat(queue{i,'File_Name'}{1},'_1-Cleaned.mat'));
        save(fName, 'trigTimes', 'trigShift', '-append')
        
        fName = fullfile(dataDirSess,'2-Filtered',strcat(queue{i,'File_Name'}{1},'_2A-CleanFiltered.mat'));
        save(fName, 'trigTimes', 'trigShift', '-append')
        
        fName = fullfile(dataDirSess,'2-Filtered',strcat(queue{i,'File_Name'}{1},'_2B-Unfiltered.mat'));
        save(fName, 'trigTimes', 'trigShift', '-append')
        
        fName = fullfile(dataDirSess,'3-Subsampled',strcat(queue{i,'File_Name'}{1},'_3B-UnfilterSubsampled.mat'));
        save(fName, 'trigTimes', 'trigShift', '-append')
    else
        disp('Skipping Save since trigTimes has no change.')
    end
    

    
    %% Re-Epoch, splitting files as necessary
    disp('Re-Epoching...')
    fName = fullfile(dataDirSess,'3-Subsampled',strcat(queue{i,'File_Name'}{1},'_3A-CleanSubsampled.mat'));
    load(fName);
    ftData_orig = ftData;
    trigTimes_orig = trigTimes;
    split_points = [];
    
    % Check if "Split_Points" is text or numeric
    if iscellstr(queue{i,'Split_Points'})
        % "Split_Points" is text, split it
        if ~isempty(queue{i,'Split_Points'}{1})
            split_points = str2double(split(queue{i,'Split_Points'}{1},';'));
        end
    else
        % "Split_Points" is numeric, assign it directly
        split_points = [queue{i,'Split_Points'}];
    end

    if ~isempty(queue{i,'Split_Points'}) 
        split_points = [0; split_points; ftData.time{1}(end-1)]; % Need to add the start and end
    else
        split_points = [0; ftData.time{1}(end-1)];
    end
    conditions = split(queue{i,'Condition_Names'}{1},';');
    
    %% Extract only triggers within each window
    for window = 1:length(split_points) - 1
        
        %% Extract trigger times that are only within the window
        window_start = split_points(window);
        window_end = split_points(window + 1);
        mask = (trigTimes_orig >= window_start) & (trigTimes_orig <= window_end);
        
        trigTimes = trigTimes_orig(mask);
        if isempty(trigTimes)
            trigTimes = window_start + (1:50) * 2;
        end
        
        %% Now perform the Epoching
        
        trigIdx = round(trigTimes * ftData.fsample);
        trialDef = zeros(length(trigIdx),3);
        trialDef(:,1) = round(trigIdx - ftData.fsample * pre_epoch); % Starting Index of Epoch
        trialDef(:,2) = round(trigIdx + ftData.fsample * length_epoch); % Ending Index
        trialDef(:,3) = -round(ftData.fsample * pre_epoch); % How far past the trigger the start of the epoch is (hence negative for starting before trigger)
        
        cfg = [];
        cfg.trl = trialDef;
        ftData = ft_redefinetrial(cfg,ftData_orig);
        
        %% SAVE EPOCHED RESULTS
        saveFile = fullfile(dataDirSess,'4-Epoched');
        mkdir(saveFile)
        
        if isempty(conditions) || length(conditions) == 1
            stem = strcat(queue{i,'File_Name'}{1},...
                '_4A-CleanEpoched.mat');
        else
            stem = strcat(num2str(queue{i,'Patient'}),'_',...
                conditions{window}, '_',...
                queue{i,'Folder_Name'}{1},...
                '_4A-CleanEpoched.mat');
        end
        
        saveFile = fullfile(saveFile,stem);
        saveFailed = true;
        attempts = 0;
        
        while saveFailed
            try
                save(saveFile, 'sessionDir', 'ftData', 'trigTimes', 'trigLengths', 'channelInfo', '-v7.3')
                saveFailed = false;
            catch ME
                saveFailed = true;
                attempts = attempts + 1;
                if attempts < 5
                    disp('Save Failed, Attempting Again...')
                else
                    disp('Save Failed, Giving up')
                    breakN
                end
            end
        end
    end
    
    %% REPEAT FOR UNFILTERED
    fName = fullfile(dataDirSess,'3-Subsampled',strcat(queue{i,'File_Name'}{1},'_3B-UnfilterSubsampled.mat'));
    load(fName,'ftData')
    ftData_orig = ftData;
    
    %% Extract only triggers within each window
    for window = 1:length(split_points) - 1
        window_start = split_points(window);
        window_end = split_points(window + 1);
        mask = (trigTimes_orig >= window_start) & (trigTimes_orig <= window_end);
        
        trigTimes = trigTimes_orig(mask);
        
        %% Now perform the Epoching
        pre_epoch = 1.0; % Length to integrate beforehand
        length_epoch = 2.0; % Length of epoch in Seconds
        
        trigIdx = round(trigTimes * ftData.fsample);
        trialDef = zeros(length(trigIdx),3);
        trialDef(:,1) = round(trigIdx - ftData.fsample * pre_epoch); % Starting Index of Epoch
        trialDef(:,2) = round(trigIdx + ftData.fsample * length_epoch); % Ending Index
        trialDef(:,3) = -round(ftData.fsample * pre_epoch); % How far past the trigger the start of the epoch is (hence negative for starting before trigger)
        
        cfg = [];
        cfg.trl = trialDef;
        ftData = ft_redefinetrial(cfg,ftData_orig);
        
        %% SAVE EPOCHED RESULTS
        saveFile = fullfile(dataDirSess,'4-Epoched');
        mkdir(saveFile)
        
        if isempty(conditions{window})
            stem = strcat(queue{i,'File_Name'}{1},...
            '_4B-UnfilterEpoched.mat');
        else
            stem = strcat(num2str(queue{i,'Patient'}),'_',...
                conditions{window}, '_',...
                queue{i,'Folder_Name'}{1},...
                '_4B-UnfilterEpoched.mat');
        end
        
        saveFile = fullfile(saveFile,stem);
        saveFailed = true;
        attempts = 0;
        
        while saveFailed
            try
                save(saveFile, 'sessionDir', 'ftData', 'trigTimes', 'trigLengths', 'channelInfo', '-v7.3')
                saveFailed = false;
            catch ME
                saveFailed = true;
                attempts = attempts + 1;
                if attempts < 5
                    disp('Save Failed, Attempting Again...')
                else
                    disp('Save Failed, Giving up')
                    break
                end
            end
        end
    end
end