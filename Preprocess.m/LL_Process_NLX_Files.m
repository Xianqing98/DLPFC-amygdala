%% START UP FILE FOR TMS-iEEG ANALYSIS
% DANNY HUANG, JEFF WANG, COREY KELLER
% 12/20/2022 XIANQING LIU
% 10/10/2023 ZHUORAN LI

% REQUIRED TOOLBOXES
% FIELDTRIP
% NEURALYNX V6.0
% DANNY's IEP SCRIPTS

clear
close all

%% Set up directories
dataDir = '/Users/zli328/Documents/Research/Research/Project_Hippocampus/raw_data';
queueFName = '/Users/zli328/Documents/Research/Research/Project_Hippocampus/raw_data/queque_file/LL_Process_Queue_test.xlsx';
saveDir =  '/Users/zli328/Documents/Research/Research/Project_Hippocampus/preprocess';

%% Parameters to use
preBuffer = 0.01; % How much pre and post to exclude from trigger for artifact
postBuffer = 0.025;
basefreq = 60; % Base Frequency for Line Noise
harmonics = 7; % How many harmonics of line noise to include in bandstop
bpfreqs = [1,35]; %Bandpass filter frequency cutoffs % mark: [2 for detrend]

%% ESTABLISH DIRECTORY AND ADD PATH
addpath /Users/zli328/Library/CloudStorage/OneDrive-UniversityofIowa/1_TMS_iEEG/0_Scripts/fieldtrip-20230206
ft_defaults         

% cd(fileparts(which(mfilename))); %Navigate to the script's location
addpath(genpath('/Users/zli328/Library/CloudStorage/OneDrive-UniversityofIowa/1_TMS_iEEG/0_Scripts/Jeff_Scripts'), ...
    genpath('./dannyIEPScripts'), ...
    genpath('./dannyGenScripts'), ...
    genpath('./MatlabImportExport_v6.0.0'), ...
    genpath('./util'));

%% Import file names to be processed
queue = readtable(queueFName);
trigShift = 0;

for i = 1:height(queue)
    disp("------Working on Patient " + num2str(queue{i,'Patient'}) + ...
        ", Session " + queue{i,'Session_Name'}{1} + " ------")
    
    %Save Directory for everything
    saveDirSess = fullfile(saveDir, num2str(queue{i,'Patient'}), queue{i,'Folder_Name'}{1});
    
    %% Check the status of processing
    % 0: Nothing Done
    % 1: Raw Imported
    % 2: Cleaned
    % 3: Filtered
    % 4: Subsampled
    % 5: Epoched/Finished
    
    if queue{i,'Completed'} == 5
        disp('Patient already processed. Moving on.')
        continue;
    end
    total_time = tic;
    
    %% Get the Directory of NLX Files
    sessionDir = fullfile(dataDir, num2str(queue{i,'Patient'}), 'NLX files', queue{i,'Session_Name'});
    files = dir(sessionDir{1});
    dirFlags = [files.isdir];
    subFolders = files(dirFlags);
    sessionDir = fullfile(sessionDir, subFolders(3).name);
    lfpDir = fullfile(sessionDir,'LFP');
    
    %% Random Initializations
    
    if queue{i,'Completed'} == 0
        %% LOAD NLX FILES
        disp('Loading NLX Files')
        subparttime = tic;
        ftData = jw_import_neuralynx(lfpDir);
        sessionDirOld = sessionDir;
        
        
        %% EXTRACT TRIGGER TIMINGS
        if isempty(queue{i,'TrigFile'}{1})
            trigFile = fullfile(sessionDir{1},'Mat','Filtered_subset.mat');
            load(trigFile,'trigTimes','trigLengths')
        elseif strcmp(queue{i,'TrigFile'}{1},'self')
            trigFile = fullfile(saveDirSess,'3-Subsampled');
            trigFile = fullfile(trigFile,strcat(queue{i,'File_Name'}{1},'_3A-CleanSubsampled.mat'));
            load(trigFile,'trigTimes','trigLengths')
        else
            trigFile = fullfile(sessionDir{1},queue{i,'TrigFile'}{1});
            [trigTimes, trigLengths] = jw_extract_trig_times(trigFile,queue{i,'Threshold'},false);
            trigTimes = trigTimes - 0.001;
            trigDiffs = [2,diff(trigTimes)];
            trigTimes = trigTimes(trigDiffs > 1); 
            trigTimes = trigTimes(trigTimes>0);
        end
        
        
        %% LOAD CHANNEL INFO
%         channelInfoFName = fullfile(dataDir, num2str(queue{i,'Patient'}), ...
%             strcat(num2str(queue{i,'Patient'}),'_Contacts_surfaceRAS.csv')); % mark
        channelInfoFName = fullfile(dataDir, num2str(queue{i,'Patient'}), ...
            strcat(num2str(queue{i,'Patient'}),'_Contacts.csv'));
        
        % Select only for channels that are connected
        channelInfo = readtable(channelInfoFName);
        validChannels = strtrim(cellstr(num2str(channelInfo.Channel)));
        validChannels = strcat('LFPx',validChannels);
        
        % Now select subset
        cfg = [];
        cfg.channel = validChannels;
        ftData = ft_selectdata(cfg, ftData);
        
        subparttime = toc(subparttime) / 60;
        disp(strcat(num2str(subparttime), ' Minutes have elapsed'))
        
        
        %% Plot an example channel and save
        saveFile = fullfile(saveDirSess,'5-Figures');
        mkdir(saveFile)
       
        if ~isempty(queue{i,'Example_Channel'}{1})
            split_channel = split(queue{i,'Example_Channel'}{1},';');
            channel_num = size(split_channel,1);
            
            for example_i = 1:channel_num
                for chn_i = 1:length(ftData.label)
                    if strcmp(ftData.label{chn_i,1},split_channel{example_i,1})
                        savename = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_0-Raw_',split_channel{example_i,1},'.png'));
                        plot_triggers(ftData,trigTimes,savename,chn_i);
                    end
                end
            end
        else
            savename = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_0-Raw_',ftData.label{1,1},'.png'));
            plot_triggers(ftData,trigTimes,savename,1); % if no example channel specified, plot the first channel in ftData.label
        end

        
        %% SAVE RAW FILES
        disp('Saving Raw Mat')
        subparttime = tic;
        
        saveFile = fullfile(saveDirSess,'0-Raw');
        mkdir(saveFile)
        saveFile = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_0-Raw.mat'));
        
        saveFailed = true;
        attempts = 0;
        channelInfoOrig = channelInfo;
        while saveFailed
            try
                save(saveFile, 'sessionDir', 'ftData', 'trigTimes', 'trigLengths', 'channelInfo', 'trigShift', 'channelInfoOrig', '-v7.3')
                saveFailed = false;
            catch ME
                saveFailed = true;
                attempts = attempts + 1;
                if attempts < 5
                    disp('Save Failed, Attempting Again...')
                    pause(5);
                else
                    disp('Save Failed, Giving up')
                    break
                end
            end
        end
        
        subparttime = toc(subparttime) / 60;
        disp(strcat(num2str(subparttime), ' minutes have elapsed'))
        
    elseif queue{i,'Completed'} == 1
        disp('Loading Raw MAT Directly')
        %% Prep to load raw files if need to clean
        loadFile = fullfile(saveDirSess,'0-Raw');
        loadFile = fullfile(loadFile,strcat(queue{i,'File_Name'}{1},'_0-Raw.mat'));
        load(loadFile)
        ftData = reselect_channels(ftData,channelInfo);
        
        %% Plot an example channel and save
        saveFile = fullfile(saveDirSess,'5-Figures');
        mkdir(saveFile)
       
        if ~isempty(queue{i,'Example_Channel'}{1})
            split_channel = split(queue{i,'Example_Channel'}{1},';');
            channel_num = size(split_channel,1);
            
            for example_i = 1:channel_num
                for chn_i = 1:length(ftData.label)
                    if strcmp(ftData.label{chn_i,1},split_channel{example_i,1})
                        savename = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_0-Raw_',split_channel{example_i,1},'.png'));
                        plot_triggers(ftData,trigTimes,savename,chn_i);
                    end
                end
            end
        else
            savename = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_0-Raw_',ftData.label{1,1},'.png'));
            plot_triggers(ftData,trigTimes,savename,1); % if no example channel specified, plot the first channel in ftData.label
        end
        
        
    end
    
    
    %% ARTIFACT REJECTION
    if queue{i,'Completed'} <= 1

        %% Remove Channels
        disp('Removing Channels...')
        subparttime = tic;
        remove_channels = queue{i,'Remove_Channels'};
        mask = zeros(size(trigTimes));
        idx_remove = [];
        
        % If only one window to remove, then will be numerical
        if isnumeric(remove_channels)
            if ~isnan(remove_channels) % mark
                idx_remove = find(channelInfo{:,'Channel'} == remove_channels);
                channelInfo(idx_remove,:) = [];
            end
        else
            % Otherwise need to split by semicolons
            remove_channels= remove_channels{1};
            
            if ~isempty(remove_channels)
                remove_channels = str2double(split(remove_channels,';'));
                idx_remove = find(ismember(channelInfo{:,'Channel'},remove_channels));
                channelInfo(idx_remove,:) = [];
            end
        end
        
        disp(strcat("Removing Channels: ",num2str(idx_remove)));
        ftData = reselect_channels(ftData,channelInfo);
        
        subparttime = toc(subparttime) / 60;
        disp(strcat(num2str(subparttime), ' minutes have elapsed'))
        channelInfoOrig = channelInfo;
        
        %% Re-reference
        disp('Rereferencing')
        [cortical_labels, bipolar_labelold, bipolar_labelnew, bipolar_tra, channelInfo] = generate_reref_matrix(channelInfoOrig);
        
        % Common Rereference for Corticals
        % Now select subset
        cfg = [];
        cfg.channel = cortical_labels;
        reref_cortical = ft_selectdata(cfg, ftData);
        
        if ~isempty(cortical_labels)
            cfg = [];
            cfg.reref       = 'yes';
            cfg.refchannel  = 'all';
            cfg.refmethod   = 'avg';
            
            reref_cortical = ft_preprocessing(cfg, reref_cortical);
        end
        
        %% Bipolar Montage for Depths
        cfg = [];
        cfg.channel = bipolar_labelold;
        reref_depth = ft_selectdata(cfg, ftData);
        
        if ~isempty(bipolar_labelold)
            cfg = [];
            cfg.channel = bipolar_labelold;
            cfg.montage = [];
            cfg.montage.tra = bipolar_tra;
            cfg.montage.labelold = bipolar_labelold;
            cfg.montage.labelnew = bipolar_labelnew;
            reref_depth = ft_preprocessing(cfg, reref_depth);
        end
        
        % Now append the Cortical and Depth Electrodes
        cfg = [];
        ftData = ft_appenddata(cfg, reref_cortical, reref_depth);
        
        
        %% REJECT ARTIFACTS
        disp('Rejecting Artifacts')
        subparttime = tic;
        
        % Clean with SIGNI interpolation
        ftData = dh_cleanartifact_interp(ftData,trigTimes, preBuffer, postBuffer,1); %Reject 25 ms after the trigger
        
        subparttime = toc(subparttime) / 60;
        disp(strcat(num2str(subparttime), ' Minutes have elapsed During Artifact Rejection'))
        
        saveFile = fullfile(saveDirSess,'5-Figures');
        mkdir(saveFile)
       
        if ~isempty(queue{i,'Example_Channel'}{1})
            split_channel = split(queue{i,'Example_Channel'}{1},';');
            channel_num = size(split_channel,1);
            
            for example_i = 1:channel_num
                for chn_i = 1:length(ftData.label)
                    if strcmp(ftData.label{chn_i,1},split_channel{example_i,1})
                        savename = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_1-Clean_',split_channel{example_i,1},'.png'));
                        plot_triggers(ftData,trigTimes,savename,chn_i);
                    end
                end
            end
        else
            savename = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_1-Clean_',ftData.label{1,1},'.png'));
            plot_triggers(ftData,trigTimes,savename,1); % if no example channel specified, plot the first channel in ftData.label
        end
       
        %% SAVE CLEANED RESULTS
        disp('Saving Cleaned MAT')
        subparttime = tic;
        
        saveFile = fullfile(saveDirSess,'1-Cleaned');
        mkdir(saveFile)
        saveFile = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_1-Cleaned.mat'));
        
        saveFailed = true;
        attempts = 0;
        
        while saveFailed
            try
                save(saveFile, 'sessionDir', 'ftData', 'trigTimes', 'trigLengths', 'channelInfo', 'trigShift', '-v7.3')
                saveFailed = false;
            catch ME
                saveFailed = true;
                attempts = attempts + 1;
                if attempts < 5
                    disp('Save Failed, Attempting Again...')
                    pause(5);
                else
                    disp('Save Failed, Giving up')
                    break
                end
            end
        end
        
        subparttime = toc(subparttime) / 60;
        disp(strcat(num2str(subparttime), ' Minutes have elapsed'))
        
    elseif queue{i,'Completed'} == 2
        disp('Loading Cleaned MAT Directly')
        
        %% Prep to load cleaned files if need to filter
        loadFile = fullfile(saveDirSess,'1-Cleaned');
        loadFile = fullfile(loadFile,strcat(queue{i,'File_Name'}{1},'_1-Cleaned.mat'));
        load(loadFile)
       
        if ~isempty(queue{i,'Example_Channel'}{1})
            split_channel = split(queue{i,'Example_Channel'}{1},';');
            channel_num = size(split_channel,1);
            
            for example_i = 1:channel_num
                for chn_i = 1:length(ftData.label)
                    if strcmp(ftData.label{chn_i,1},split_channel{example_i,1})
                        savename = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_1-Clean_',split_channel{example_i,1},'.png'));
                        plot_triggers(ftData,trigTimes,savename,chn_i);
                    end
                end
            end
        else
            savename = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_1-Clean_',ftData.label{1,1},'.png'));
            plot_triggers(ftData,trigTimes,savename,1); % if no example channel specified, plot the first channel in ftData.label
        end
    end
    
    %% FILTERING
    if queue{i,'Completed'} <= 2
        disp('Filtering Results')
        subparttime = tic;
        cfg = [];
        
        % Bandpass Filter
        cfg.bpfilter = 'yes';
        cfg.bpfiltord = 3;
        cfg.bpfreq = [bpfreqs];
        
        % NOTCH FILTER
        cfg.bsfilter       = 'yes';
        cfg.bsfiltord      = 3;
        line_noise_freqs   = basefreq * (1:harmonics);
        cfg.bsfreq         = [line_noise_freqs - 3; line_noise_freqs + 3].'; % use this for NY patients
        cfg.bsfreq = [cfg.bsfreq];
        
        ftData = ft_preprocessing(cfg,ftData); 
        subparttime = toc(subparttime) / 60;
        disp(strcat(num2str(subparttime), ' Minutes have elapsed'))
        
        %% Plot an example channel and save
        saveFile = fullfile(saveDirSess,'5-Figures');
        mkdir(saveFile)
        
        if ~isempty(queue{i,'Example_Channel'}{1})
            split_channel = split(queue{i,'Example_Channel'}{1},';');
            channel_num = size(split_channel,1);
            
            for example_i = 1:channel_num
                for chn_i = 1:length(ftData.label)
                    if strcmp(ftData.label{chn_i,1},split_channel{example_i,1})
                        savename = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_2-Filter_',split_channel{example_i,1},'.png'));
                        plot_triggers(ftData,trigTimes,savename,chn_i);
                    end
                end
            end
        else
            savename = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_2-Filter_',ftData.label{1,1},'.png'));
            plot_triggers(ftData,trigTimes,savename,1); % if no example channel specified, plot the first channel in ftData.label
        end
        
        %% SAVE FILTERED RESULTS
        disp('Saving Filtered MAT')
        subparttime = tic;
        
        saveFile = fullfile(saveDirSess,'2-Filtered');
        mkdir(saveFile)
        saveFile = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_2A-CleanFiltered.mat'));
        
        saveFailed = true;
        attempts = 0;
        
        while saveFailed
            try
                save(saveFile, 'sessionDir', 'ftData', 'trigTimes', 'trigLengths', 'channelInfo', 'trigShift', '-v7.3')
                saveFailed = false;
            catch ME
                saveFailed = true;
                attempts = attempts + 1;
                if attempts < 5
                    disp('Save Failed, Attempting Again...')
                    pause(5);
                else
                    disp('Save Failed, Giving up')
                    break
                end
            end
        end
        
        subparttime = toc(subparttime) / 60;
        disp(strcat(num2str(subparttime), ' Minutes have elapsed'))
    elseif queue{i,'Completed'} == 3
        disp('Loading Filtered MAT Directly')
        %% Prep to load cleaned files if need to filter
        loadFile = fullfile(saveDirSess,'2-Filtered');
        loadFile = fullfile(loadFile,strcat(queue{i,'File_Name'}{1},'_2A-CleanFiltered.mat'));
        load(loadFile)

        if ~isempty(queue{i,'Example_Channel'}{1})
            split_channel = split(queue{i,'Example_Channel'}{1},';');
            channel_num = size(split_channel,1);
            
            for example_i = 1:channel_num
                for chn_i = 1:length(ftData.label)
                    if strcmp(ftData.label{chn_i,1},split_channel{example_i,1})
                        savename = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_2-Filter_',split_channel{example_i,1},'.png'));
                        plot_triggers(ftData,trigTimes,savename,chn_i);
                    end
                end
            end
        else
            savename = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_2-Filter_',ftData.label{1,1},'.png'));
            plot_triggers(ftData,trigTimes,savename,1); % if no example channel specified, plot the first channel in ftData.label
        end
        
    end
    
    if queue{i,'Completed'} <= 3
        %% SUBSAMPLE
        cfg = [];
        cfg.resamplefs = 1000;
        ftData = ft_resampledata(cfg,ftData);
        
        
        %% Plot an example channel and save
        saveFile = fullfile(saveDirSess,'5-Figures');
        mkdir(saveFile)
 
        if ~isempty(queue{i,'Example_Channel'}{1})
            split_channel = split(queue{i,'Example_Channel'}{1},';');
            channel_num = size(split_channel,1);
            
            for example_i = 1:channel_num
                for chn_i = 1:length(ftData.label)
                    if strcmp(ftData.label{chn_i,1},split_channel{example_i,1})
                        savename = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_3-Sampled_',split_channel{example_i,1},'.png'));
                        plot_triggers(ftData,trigTimes,savename,chn_i);
                    end
                end
            end
        else
            savename = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_3-Sampled_',ftData.label{1,1},'.png'));
            plot_triggers(ftData,trigTimes,savename,1); % if no example channel specified, plot the first channel in ftData.label
        end
        
        
        %% SAVE SUBSAMPLED RESULTS
        saveFile = fullfile(saveDirSess,'3-Subsampled');
        mkdir(saveFile)
        saveFile = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_3A-CleanSubsampled.mat'));
        
        saveFailed = true;
        attempts = 0;
        
        while saveFailed
            try
                save(saveFile, 'sessionDir', 'ftData', 'trigTimes', 'trigLengths', 'channelInfo', 'trigShift', '-v7.3')
                saveFailed = false;
            catch ME
                saveFailed = true;
                attempts = attempts + 1;
                if attempts < 5
                    disp('Save Failed, Attempting Again...')
                    pause(5);
                else
                    disp('Save Failed, Giving up')
                    break
                end
            end
        end
        
        
    elseif queue{i,'Completed'} == 4
        disp('Loading Subsampled MAT Directly')
        %% Prep to load cleaned files if need to filter
        loadFile = fullfile(saveDirSess,'3-Subsampled');
        loadFile = fullfile(loadFile,strcat(queue{i,'File_Name'}{1},'_3A-CleanSubsampled.mat'));
        load(loadFile)
           
        if ~isempty(queue{i,'Example_Channel'}{1})
            split_channel = split(queue{i,'Example_Channel'}{1},';');
            channel_num = size(split_channel,1);
            
            for example_i = 1:channel_num
                for chn_i = 1:length(ftData.label)
                    if strcmp(ftData.label{chn_i,1},split_channel{example_i,1})
                        savename = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_3-Sampled_',split_channel{example_i,1},'.png'));
                        plot_triggers(ftData,trigTimes,savename,chn_i);
                    end
                end
            end
        else
            savename = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_3-Sampled_',ftData.label{1,1},'.png'));
            plot_triggers(ftData,trigTimes,savename,1); % if no example channel specified, plot the first channel in ftData.label
        end
        
    end
    
    if queue{i,'Completed'} <= 4
        disp('Epoching...')
        %% EPOCH DATA
        subparttime = tic;
        pre_epoch = 1.0; % Length to integrate beforehand
        length_epoch = 2.0; % Length of epoch in Seconds
        
        trigIdx = round(trigTimes * ftData.fsample);
        trialDef = zeros(length(trigIdx),3);
        trialDef(:,1) = round(trigIdx - ftData.fsample * pre_epoch); % Starting Index of Epoch
        trialDef(:,2) = round(trigIdx + ftData.fsample * length_epoch); % Ending Index
        trialDef(:,3) = -round(ftData.fsample * pre_epoch); % How far past the trigger the start of the epoch is (hence negative for starting before trigger)
        
        cfg = [];
        cfg.trl = trialDef;
        ftData_new = ft_redefinetrial(cfg,ftData);
%         ftData = ft_redefinetrial(cfg,ftData);
        
        %% SAVE EPOCHED RESULTS
        saveFile = fullfile(saveDirSess,'4-Epoched');
        mkdir(saveFile)
        saveFile = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_4A-CleanEpoched.mat'));
        
        saveFailed = true;
        attempts = 0;
        
        while saveFailed
            try
                save(saveFile, 'sessionDir', 'ftData', 'trigTimes', 'trigLengths', 'channelInfo', 'trigShift', '-v7.3')
                saveFailed = false;
            catch ME
                saveFailed = true;
                attempts = attempts + 1;
                if attempts < 5
                    disp('Save Failed, Attempting Again...')
                    pause(5);
                else
                    disp('Save Failed, Giving up')
                    break
                end
            end
        end
        
        subparttime = toc(subparttime) / 60;
        disp(strcat(num2str(subparttime), ' Minutes have elapsed'))
        
    end
    
    total_time = toc(total_time) / 60;
    disp(strcat('Patient finished in ', num2str(total_time), ' minutes'));
end
