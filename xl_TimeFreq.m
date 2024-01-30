%% TIME-FREQ ANALYSIS - for iES-iEEG
% Xianqing Bella Liu
% 09/17/2023

% REQUIRED TOOLBOXES
% FIELDTRIP
clear; clc; 
close all;

rootpath = '/Users/xianqliu/Library/CloudStorage/OneDrive-UniversityofIowa/1_TMS_iEEG';
cd(rootpath);
addpath /Users/xianqliu/Documents/MATLAB/fieldtrip-20230206
addpath(genpath('./0_Scripts'), ...
        genpath('./3_ProcessedData_Backup'));

StimType = 'esTT_L-DLPFC';
patient_list = [460 534 625 634];
chlinfo = readtable([rootpath '/3_ProcessedData_Backup/Selected_channel_AMY.xlsx'], 'Sheet', 'esTT L-DLPFC');

ft_defaults;

%% 1. Load Data for Each Patient
% Create a cell array to hold the data from each patient
all_patients_data = cell(1, length(patient_list)); % replace numPatients with the actual number of patients

% Loop through each patient's dataset
for i = 1:length(patient_list)
    
    % Load data (Replace this line with actual code to load data)
    data = load([rootpath '/3_ProcessedData_Backup/xl_Data_Processed/' StimType '/' ...
        num2str(patient_list(i)) '_' StimType '_0_5Hz_4B-UnfilterEpoched.mat']);
    
    % Store in cell array
    all_patients_data{i} = data;
end
clear data ftData trigLengths trigShift trigTimes sessionDir


%% 2. Select Channels and Trials, and Combine Channels Across Patients

% specify channels in amyg 
if strcmp(StimType, 'esTT_L-DLPFC')
    amygChannel{1} = {'LFPx105' 'LFPx106' 'LFPx107'};
    amygChannel{2} = {'LFPx228' 'LFPx229' 'LFPx230' 'LFPx231'};
    amygChannel{3} = {'LFPx132' 'LFPx133' 'LFPx134'};
    amygChannel{4} = {'LFPx57'};
end
amygChannelAll = [amygChannel{1} amygChannel{2} amygChannel{3} amygChannel{4}];
nChannel = length(amygChannelAll);

% Apply the channel selection to each patient's data
for i = 1:length(patient_list)
    cfg = [];
    cfg.channel = amygChannel{i}; % amyg channel each patient

    all_patients_data{i} = ft_selectdata(cfg, all_patients_data{i}.ftData);
end


% Append the data across patients
cfg = [];
combined_data = ft_appenddata(cfg, all_patients_data{:});
combined_data_org = combined_data;

removeTrials = cell(1,nChannel);
for i = 1:nChannel
    removeTrials{i} = str2double(strsplit(chlinfo{i,4}{1}, ';')); 
    removeTrials{i} = removeTrials{i}(~isnan(removeTrials{i}));
end

% Loop through each channel and set the data in the trials to be removed to NaN
for i = 1:nChannel % assuming patient_data is your FieldTrip data structure
    % Check if there are trials to remove for this channel
    if ~isnan(removeTrials{i})
        % Loop over the trials and set the data values to NaN
        for trl = removeTrials{i}
            patient_data.trial{trl}(i, :) = NaN;
        end
    end
end


% For each channel
singleChannelData = cell(1,nChannel);
for i = 1:nChannel
    cfg = [];
    cfg.channel = amygChannelAll{i}; % Select the specific channel
    cfg.trials = ones([1 length(combined_data.trial)]); 
    trl_deleted = str2double(strsplit(chlinfo{i,4}{1}, ';')); 
    trl_deleted = trl_deleted(~isnan(trl_deleted));
    cfg.trials(trl_deleted) = 0; 
    cfg.trials = logical(cfg.trials); % logical value of whether the trials should be kept
    singleChannelData{i} = ft_selectdata(cfg, combined_data);
    % Now singleChannelData will have only the desired trials for the specific channel.
end

% Append single channel data
cfg = [];
combined_data = ft_appenddata(cfg, singleChannelData{:});


% Medial and lateral division
cfg = [];
cfg.channel = {'LFPx105' 'LFPx228' 'LFPx229' 'LFPx230' 'LFPx231'};
medial_data = ft_selectdata(cfg, combined_data);

cfg = [];
cfg.channel = {'LFPx106' 'LFPx107' 'LFPx132' 'LFPx133' 'LFPx134' 'LFPx57'};
lateral_data = ft_selectdata(cfg, combined_data);


%% 3. Time-Frequency Analysis
freq = cell(1,nChannel);
freq_blc = cell(1,nChannel);
for i = 1:nChannel
    cfg = [];
    cfg.output     = 'pow';
    cfg.method     = 'wavelet';
    % cfg.method     = 'mtmfft';
    % cfg.method     = 'mtmconvol';
    % cfg.taper      = 'dpss';
    % cfg.tapsmofrq  = 2;  % Smoothing parameter for the multitaper
    % cfg.t_ftimwin  = 0.5 * ones(size(cfg.foi));  % length of time window, which typically scales with frequency
    cfg.width      = 7;
    cfg.foi        = 2:1:200;
    cfg.toi        = -0.2:0.005:1.5;
    cfg.keeptrials = 'yes';
    freq{i} = ft_freqanalysis(cfg, singleChannelData{i});
    
    % Baseline Correction
    cfg = [];
    cfg.baseline = [-0.2 -0.01];
    cfg.baselinetype = 'relative';
    freq_blc{i} = ft_freqbaseline(cfg, freq{i});
end



%% 4.1 Multi-Channel Power Spectrum Plot
for i = 1:nChannels
    subplot(ceil(sqrt(nChannels)), ceil(sqrt(nChannels)), i);
    
    % Extract Data for the channel
    data = squeeze(mean(freq_blc{i}.powspctrm(:, 1, :, :), 1));
    
    % Log-scale frequencies
    log_frequencies = log(freq_blc{i}.freq);
    
    % Plot
    imagesc(freq_blc{i}.time, log_frequencies, data);
    
    % Flip y-axis
    set(gca, 'YDir', 'normal');
    
    % Customize Axes
    yticks(log([2, 4, 10, 30, 70, 150]));
    yticklabels({'2', '4', '10', '30', '70', '150'});
    
    % Set color bar limits
    caxis([-0. 0.5]);
    
    title(['Channel: ' amygChannelAll{i}]);
end

colorbar;


%% 4.2 Average power spectrum across channels
% average_powspctrm - a matrix of size [N_Frequencies x N_Times]
average_powspctrm = nanmean(freq_blc{:}.powspctrm, 1);
average_powspctrm = nanmean(average_powspctrm, 2);
average_powspctrm = squeeze(average_powspctrm);

% Plot averaged power spectrum
% Log-scale frequencies
log_frequencies = log(freq_blc.freq);

figure;  % Create a new figure for the averaged power spectrum

imagesc(freq_blc.time, log_frequencies, average_powspctrm);

% Flip y-axis
set(gca, 'YDir', 'normal');

% Customize Axes
yticks(log([2, 4, 10, 30, 70, 150]));
yticklabels({'2', '4', '10', '30', '70', '150'});

% Set color bar limits
caxis([0 5]);

% Add labels and title
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Averaged Power Spectrum Across Channels');

% Add color bar
colorbar;



%% Plot power spectrum: Multi-Channel Overview
% layout = [];
% layout.label = combined_data.label';
% layout.pos = [
%     0.1 0.9; 0.5 0.9; 0.9 0.9; 1.3 0.9;  % Row 1
%     0.1 0.6; 0.5 0.6; 0.9 0.6; 1.3 0.6;  % Row 2
%     0.1 0.3; 0.5 0.3; 0.9 0.3];
% 
% cfg.showlabels = 'yes';
% cfg.fontsize = 8;
% cfg.layout = layout; 
% cfg.colorbar = 'yes';
% cfg.ylim = [2 150];
% cfg.zlim = [-1 5];
% cfg.baseline = [-0.1 -0.01];
% 
% ft_multiplotTFR(cfg, freq_blc);



%% Plot normalized data for each channel
% % Extract the baseline window data
% baseline_indices = (freq_blc.time >= -0.2) & (freq_blc.time <= -0.01);
% 
% for ch = 1:nChannels
%     subplot(ceil(sqrt(nChannels)), ceil(sqrt(nChannels)), ch);
%     
%     % Extract Data for the channel
%     channel_data = squeeze(freq_blc.powspctrm(:, ch, :, :));
%     
%     % Compute the baseline mean for this channel
%     baseline_data = channel_data(:, :, baseline_indices);
%     baseline_mean = nanmean(baseline_data, 3);
%     
%     % Normalize the data by the baseline mean
%     normalized_data = bsxfun(@rdivide, channel_data, baseline_mean);
%     
%     % Average over trials if necessary
%     avg_normalized_data = nanmean(normalized_data, 1);
%     
%     % Log-scale frequencies for visualization
%     log_frequencies = log(freq_blc.freq);
%     
%     % Plot the normalized data
%     imagesc(freq_blc.time, log_frequencies, squeeze(avg_normalized_data));
%     
%     % Flip y-axis
%     set(gca, 'YDir', 'normal');
%     
%     % Customize Axes
%     yticks(log([2, 10, 30, 70, 150]));
%     yticklabels({'2', '10', '30', '70', '150'});
%     
%     % Set color bar limits
%     caxis([-1 8]);  % or use a different scale depending on your normalized data
%     
%     title(['Channel: ' freq_blc.label{ch}]);
% end
% 
% colorbar;
% 
% %% Plot averaged channel after normalization
% baseline_data = freq_blc.powspctrm(:,:,:,baseline_indices);
% 
% % Compute the baseline mean for each channel and frequency
% baseline_mean = nanmean(baseline_data, 4);  % The 4th dimension is time
% 
% % Normalize each channel by its own baseline mean
% normalized_data = bsxfun(@rdivide, freq_blc.powspctrm, baseline_mean);
% 
% % Average across channels (normalized)
% average_powspctrm = nanmean(normalized_data, 1);
% average_powspctrm = nanmean(average_powspctrm, 2);
% average_powspctrm = squeeze(average_powspctrm);
% 
% 
% % Plot averaged power spectrum
% % Log-scale frequencies
% log_frequencies = log(freq_blc.freq);
% 
% figure;  % Create a new figure for the averaged power spectrum
% 
% imagesc(freq_blc.time, log_frequencies, average_powspctrm);
% 
% % Flip y-axis
% set(gca, 'YDir', 'normal');
% 
% % Customize Axes
% yticks(log([2, 4, 10, 30, 70, 150]));
% yticklabels({'2', '4', '10', '30', '70', '150'});
% 
% % Set color bar limits
% caxis([-1 5]);
% 
% % Add labels and title
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Averaged Power Spectrum Across Channels');
% 
% % Add color bar
% colorbar;






