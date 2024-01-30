%% START UP FILE FOR TMS-iEEG ANALYSIS
% DANNY HUANG, JEFF WANG, COREY KELLER
% 12/20/2022 (XIANQING

% REQUIRED TOOLBOXES
% FIELDTRIP
% NEURALYNX V6.0
% DANNY's IEP SCRIPTS (JEFF I SHARED THIS WITH YOU)
clear
close all

%% Set up directories
dataDir = 'C:\Users\xjl19\Desktop\TMS-iEEG\es-TT';
queueFName = 'C:\Users\xjl19\Desktop\TMS-iEEG\Processed_xlsx\xl_Process_Queue.xlsx';
saveDir =  'C:\Users\xjl19\Desktop\TMS-iEEG\xl_Data_Processed_esTT';

%% Parameters to use
% preBuffer = 0.01; % How much pre and post to exclude from trigger for artifact
% postBuffer = 0.05;
% preRecharge = -0.47;
% postRecharge = 0.53;
% basefreq = 60; % Base Frequency for Line Noise
% harmonics = 7; % How many harmonics of line noise to include in bandstop
% bpfreqs = [1,200]; %Bandpass filter frequency cutoffs

preBuffer = 0.008; % How much pre and post to exclude from trigger for artifact
postBuffer = 0.008;
basefreq = 60; % Base Frequency for Line Noise
harmonics = 7; % How many harmonics of line noise to include in bandstop
bpfreqs = [2,200]; %Bandpass filter frequency cutoffs

%% ESTABLISH DIRECTORY AND ADD PATH
cd(fileparts(which(mfilename))); %Navigate to the script's location
addpath(genpath('.\Jeff_Scripts'), ...
    genpath('.\dannyIEPScripts'), ...
    genpath('.\dannyGenScripts'), ...
    genpath('.\MatlabImportExport_v6.0.0'), ...
    genpath('.\util'));


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
    maxChannel = -1; % Initialize Max Channel as -1 to signify that it hasn't been identified yet
    
    
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
%             if trigTimes(1) < 0 % mark, error using fit
%                 trigTimes = trigTimes(2:end);
%                 trigLengths = trigLengths(2:end);
%             end
            trigDiffs = [2,diff(trigTimes)];
            trigTimes = trigTimes(trigDiffs > 1); 
            trigTimes = trigTimes(length(trigTimes)-59:end);
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
        
%         %% Shift Triggers if needed
%         trigShift = 0;
%         trigsChanged = false;
%         if ~exist('trigShift', 'var')
%             trigTimes = trigTimes + queue{i,'Trig_Shift'};
%             trigShift = queue{i,'Trig_Shift'};
%             trigsChanged = true;
%         elseif trigShift ~= queue{i,'Trig_Shift'}
%             trigTimes = trigTimes - trigShift + queue{i,'Trig_Shift'};
%             trigShift = queue{i,'Trig_Shift'};
%             trigsChanged = true;
%         end
        
        %% Plot an example channel and save
        
        saveFile = fullfile(saveDirSess,'5-Figures');
        mkdir(saveFile)
        saveFile = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_0-Raw.png'));
        
        mask = ftData.time{1} > trigTimes(1) & ftData.time{1} <= trigTimes(1) + 0.02;
        [~,maxChannel] = max(sum(abs(ftData.trial{1}(:,mask)),2));
        plot_triggers(ftData,trigTimes,saveFile,maxChannel)
        
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
        saveFile = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_0-Raw.png'));
        
        mask = ftData.time{1} > trigTimes(1) & ftData.time{1} <= trigTimes(1) + 0.02;
        [~,maxChannel] = max(sum(abs(ftData.trial{1}(:,mask)),2));
        plot_triggers(ftData,trigTimes,saveFile,maxChannel)
        
        
    end
    
    
    %% ARTIFACT REJECTION
    if queue{i,'Completed'} <= 1
%         %% Shift Triggers if needed
%         trigsChanged = false;
%         if ~exist('trigShift', 'var')
%             trigTimes = trigTimes + queue{i,'Trig_Shift'};
%             trigShift = queue{i,'Trig_Shift'};
%             trigsChanged = true;
%         elseif trigShift ~= queue{i,'Trig_Shift'}
%             trigTimes = trigTimes - trigShift + queue{i,'Trig_Shift'};
%             trigShift = queue{i,'Trig_Shift'};
%             trigsChanged = true;
%         end
%         
%         if trigsChanged
%             disp('Shifted Triggers')
%             saveFile = fullfile(saveDirSess,'0-Raw');
%             mkdir(saveFile)
%             saveFile = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_0-Raw.mat'));
%             save(saveFile, 'trigTimes', 'trigShift', '-append')
%         end
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
        
        % Corticals should use a common rereference
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
        %% Use Bipolar Montage for Depths
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
        
        %% DETREND
        disp('Detrending')
        ftData = detrend_decay(ftData,trigTimes,preBuffer,postBuffer);
        
        %% REJECT ARTIFACTS % mark
%         disp('Rejecting Artifacts')
%         subparttime = tic;
%         
%         % Save Unclean Version
%         ftData_orig = ftData;
%         
%         % Prefilter the data 
%         ftData = prefilter_data(ftData,trigTimes,preBuffer,postBuffer,basefreq,harmonics,bpfreqs,true);
%         
%         % Clean with Cubic interpolation % mark
%         ftData = dh_cleanartifact_interp(ftData,trigTimes, preBuffer, postBuffer,1); %Reject 25 ms after the trigger
%         ftData = dh_cleanartifact_interp(ftData,trigTimes, preRecharge, postRecharge,1);
%         
%         subparttime = toc(subparttime) / 60;
%         disp(strcat(num2str(subparttime), ' Minutes have elapsed During Artifact Rejection'))
        
        %% Plot an example channel and save
        saveFile = fullfile(saveDirSess,'5-Figures');
        mkdir(saveFile)
        saveFile = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_1-Cleaned.png'));

        maxChannel = 56; % mark

        plot_triggers(ftData,trigTimes,saveFile,maxChannel)
       
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
        
        %% Load uncleaned file as well
%         loadFile = fullfile(saveDirSess,'0-Raw');
%         loadFile = fullfile(loadFile,strcat(queue{i,'File_Name'}{1},'_0-Raw.mat'));
%         ftData_orig = load(loadFile,'ftData');
%         ftData_orig = ftData_orig.ftData;
        
        saveFile = fullfile(saveDirSess,'5-Figures');
        mkdir(saveFile)
        saveFile = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_1-Cleaned.png'));
        
        mask = ftData.time{1} > trigTimes(1) & ftData.time{1} <= trigTimes(1) + 0.02;
        [~,maxChannel] = max(sum(abs(ftData.trial{1}(:,mask)),2));
        plot_triggers(ftData,trigTimes,saveFile,maxChannel)
        
        if ~exist('trigShift', 'var')
            trigShift = queue{i,'Trig_Shift'};
        end
    end
    
    %% FILTERING
    if queue{i,'Completed'} <= 2
        disp('Filtering Results')

        % Save Unclean Version
        ftData_orig = ftData;
        
        subparttime = tic;
        cfg = [];
        
%         % Lowpass filter
%         cfg.lpfilter       = 'yes';
%         cfg.lpfiltord      = 3;
%         cfg.lpfreq       = bpfreqs(2);
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
%         ftData_orig = ft_preprocessing(cfg,ftData_orig);
        ftData = prefilter_data(ftData,trigTimes,preBuffer,postBuffer,basefreq,harmonics,bpfreqs,false);
%         ftData_orig = prefilter_data(ftData_orig,trigTimes,preBuffer,postBuffer + 0.01,basefreq,harmonics,bpfreqs,false);
        subparttime = toc(subparttime) / 60;
        disp(strcat(num2str(subparttime), ' Minutes have elapsed'))
        
        %% Plot an example channel and save
        saveFile = fullfile(saveDirSess,'5-Figures');
        mkdir(saveFile)
        saveFile = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_2-Filtered.png'));
        plot_triggers(ftData,trigTimes,saveFile,maxChannel)

        saveFile = fullfile(saveDirSess,'5-Figures');
        mkdir(saveFile)
        saveFile = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_2-Unfiltered.png'));
        plot_triggers(ftData_orig,trigTimes,saveFile,maxChannel)
        
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
        
        % Save uncfiltered version
        saveFile = fullfile(saveDirSess,'2-Filtered');
        mkdir(saveFile)
        saveFile = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_2B-Unfiltered.mat'));
        ftData_temp = ftData;
        ftData = ftData_orig;
        save(saveFile, 'sessionDir', 'ftData', 'trigTimes', 'trigLengths', 'channelInfo', '-v7.3')
        
        % Restore cleaned version
        ftData = ftData_temp;
        
        subparttime = toc(subparttime) / 60;
        disp(strcat(num2str(subparttime), ' Minutes have elapsed'))
    elseif queue{i,'Completed'} == 3
        disp('Loading Filtered MAT Directly')
        %% Prep to load cleaned files if need to filter
        loadFile = fullfile(saveDirSess,'2-Filtered');
        loadFile = fullfile(loadFile,strcat(queue{i,'File_Name'}{1},'_2A-CleanFiltered.mat'));
        load(loadFile)
        
        %% Load uncleaned file as well
        loadFile = fullfile(saveDirSess,'2-Filtered');
        loadFile = fullfile(loadFile,strcat(queue{i,'File_Name'}{1},'_2B-Unfiltered.mat'));
        ftData_orig = load(loadFile,'ftData');
        ftData_orig = ftData_orig.ftData;
        
        mask = ftData.time{1} > trigTimes(1) & ftData.time{1} <= trigTimes(1) + 0.02;
        [~,maxChannel] = max(sum(abs(ftData.trial{1}(:,mask)),2));
%         plot_triggers(ftData,trigTimes,saveFile,maxChannel)
    end
    
    if queue{i,'Completed'} <= 3
        %% SUBSAMPLE
        cfg = [];
        cfg.resamplefs = 1000;
        ftData = ft_resampledata(cfg,ftData);
        ftData_orig = ft_resampledata(cfg,ftData_orig);
        
        %% Plot an example channel and save
        saveFile = fullfile(saveDirSess,'5-Figures');
        mkdir(saveFile)
        saveFile = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_3-Subsampled.png'));
        plot_triggers(ftData,trigTimes,saveFile,maxChannel)

        saveFile = fullfile(saveDirSess,'5-Figures');
        mkdir(saveFile)
        saveFile = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_3-UnfilteredSubsampled.png'));
        plot_triggers(ftData_orig,trigTimes,saveFile,maxChannel)
        
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
        
        % Save Unfiltered Version
        saveFile = fullfile(saveDirSess,'3-Subsampled');
        mkdir(saveFile)
        saveFile = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_3B-UnfilterSubsampled.mat'));
        
        ftData_temp = ftData;
        ftData = ftData_orig;
        
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
        
        % Restore cleaned version
        ftData = ftData_temp;
        
    elseif queue{i,'Completed'} == 4
        disp('Loading Subsampled MAT Directly')
        %% Prep to load cleaned files if need to filter
        loadFile = fullfile(saveDirSess,'3-Subsampled');
        loadFile = fullfile(loadFile,strcat(queue{i,'File_Name'}{1},'_3A-CleanSubsampled.mat'));
        load(loadFile)
        
        %% Load uncleaned file as well
        loadFile = fullfile(saveDirSess,'3-Subsampled');
        mkdir(loadFile)
        loadFile = fullfile(loadFile,strcat(queue{i,'File_Name'}{1},'_3B-UnfilterSubsampled.mat'));
        ftData_orig = load(loadFile,'ftData');
        ftData_orig = ftData_orig.ftData;
        
        mask = ftData.time{1} > trigTimes(1) & ftData.time{1} <= trigTimes(1) + 0.02;
        [~,maxChannel] = max(sum(abs(ftData.trial{1}(:,mask)),2));
%         plot_triggers(ftData,trigTimes,saveFile,maxChannel)
    end
    
    if queue{i,'Completed'} <= 4
        disp('Epoching...')
        %% EPOCH DATA
        subparttime = tic;
        pre_epoch = 1; % Length to integrate beforehand
        length_epoch = 2; % Length of epoch in Seconds
        
        trigIdx = round(trigTimes * ftData.fsample);
        trialDef = zeros(length(trigIdx),3);
        trialDef(:,1) = round(trigIdx - ftData.fsample * pre_epoch); % Starting Index of Epoch
        trialDef(:,2) = round(trigIdx + ftData.fsample * length_epoch); % Ending Index
        trialDef(:,3) = -round(ftData.fsample * pre_epoch); % How far past the trigger the start of the epoch is (hence negative for starting before trigger)
        
        cfg = [];
        cfg.trl = trialDef;
        ftData = ft_redefinetrial(cfg,ftData);
        ftData_orig = ft_redefinetrial(cfg,ftData_orig);
        
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
        
        % Save Unfiltered Version
        saveFile = fullfile(saveDirSess,'4-Epoched');
        saveFile = fullfile(saveFile,strcat(queue{i,'File_Name'}{1},'_4B-UnfilterEpoched.mat'));
        
        ftData_temp = ftData;
        ftData = ftData_orig;
        
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
        
        % Restore cleaned version
        ftData = ftData_temp;
        
        subparttime = toc(subparttime) / 60;
        disp(strcat(num2str(subparttime), ' Minutes have elapsed'))
        
    end
    
    total_time = toc(total_time) / 60;
    disp(strcat('Patient finished in ', num2str(total_time), ' minutes'));
end
