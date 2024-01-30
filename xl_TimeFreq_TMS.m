%% TIME-FREQ ANALYSIS - for TMS-iEEG
% Xianqing Bella Liu
% 09/30/2023

% Toolbox: Fieidtrip

clear; clc; 
close all;

rootpath = '/Users/xianqliu/Library/CloudStorage/OneDrive-UniversityofIowa/1_TMS_iEEG';
cd(rootpath);
addpath /Users/xianqliu/Documents/MATLAB/fieldtrip-20230206
addpath(genpath('./0_Scripts'), ...
        genpath('./3_ProcessedData_Backup'));

StimType = 'TMS_L-DLPFC';
% patient_list = [429 430 460 483 518 524 534 593 634];
patient_list = [430 460 483 518 524 534 593 634 429];
chlinfo = readtable([rootpath '/3_ProcessedData_Backup/Selected_channel_AMY.xlsx'], 'Sheet', 'TMS L-DLPFC');

ft_defaults;

%% 1. Load Data for Each Patient
% Create a cell array to hold the data from each patient
all_patients_data = cell(1, length(patient_list)); % replace numPatients with the actual number of patients
all_patients_tms = cell(1, length(patient_list));
all_patients_sham = cell(1, length(patient_list));

% Loop through each patient's dataset
for i = 1:length(patient_list)
    
    % Load data (Replace this line with actual code to load data)
    data = load([rootpath '/3_ProcessedData_Backup/xl_Data_Processed/' StimType '/nodetrend_unfilter/' ...
        num2str(patient_list(i)) '_Comparison_L_DLPFC_0_5Hz_nodetrend_4B-UnfilterEpoched.mat']);
    
    % Store in cell array
    all_patients_data{i} = data;
end
clear data ftData trigLengths trigShift trigTimes sessionDir


%% 2. Select Channels and Trials, and Combine Channels Across Patients

% specify channels in amyg 
if strcmp(StimType, 'TMS_L-DLPFC')
    amygChannel{1} = {'LFPx43' 'LFPx44' 'LFPx8' 'LFPx9' 'LFPx10' 'LFPx11' 'LFPx12' 'LFPx13'}; % 430
    amygChannel{2} = {'LFPx105' 'LFPx106' 'LFPx107'}; % 460
    amygChannel{3} = {'LFPx185' 'LFPx186' 'LFPx187'}; % 483
    amygChannel{4} = {'LFPx27' 'LFPx29'}; % 518
    amygChannel{5} = {'LFPx149'}; % 524
    amygChannel{6} = {'LFPx231' 'LFPx229' 'LFPx230' 'LFPx228'}; % 534
    amygChannel{7} = {'LFPx233' 'LFPx186' 'LFPx187' 'LFPx234'}; % 593
    amygChannel{8} = {'LFPx57'}; % 634
    amygChannel{9} = {'LFPx120' 'LFPx121' 'LFPx122' 'LFPx123'}; % 429
end
amygChannelAll = [amygChannel{1} amygChannel{2} amygChannel{3} amygChannel{4} amygChannel{5} amygChannel{6} amygChannel{7} amygChannel{8} amygChannel{9}];
% amygChannelAll = [amygChannel{1} amygChannel{2} amygChannel{3} amygChannel{4} amygChannel{5} amygChannel{6} amygChannel{7} amygChannel{8}];
nChannel = length(amygChannelAll);

% Apply the channel selection to each patient's data
% TMS
for i = 1:length(patient_list)
    cfg = [];
    cfg.channel = amygChannel{i}; % amyg channel each patient

    all_patients_tms{i} = ft_selectdata(cfg, all_patients_data{i}.ftData_epoch_tms);
end

for i = 1:length(patient_list)
    cfg = [];
    cfg.channel = amygChannel{i}; % amyg channel each patient

    all_patients_sham{i} = ft_selectdata(cfg, all_patients_data{i}.ftData_epoch_sham);
end

% For each channel
singleChanneltms = cell(1,nChannel);
singleChannelsham = cell(1,nChannel);
allChannelIndex = 0;
for p = 1:length(patient_list)
    for n = 1:length(amygChannel{p})
        cfg = [];
        cfg.channel = amygChannel{1,p}{1,n}; % Select the specific channel
        singleChanneltms{1,n+allChannelIndex} = ft_selectdata(cfg, all_patients_tms{p});
    end
    allChannelIndex = allChannelIndex + n;
end

allChannelIndex = 0;
for p = 1:length(patient_list)
    for n = 1:length(amygChannel{p})
        cfg = [];
        cfg.channel = amygChannel{1,p}{1,n}; % Select the specific channel
        singleChannelsham{1,n+allChannelIndex} = ft_selectdata(cfg, all_patients_sham{p});
    end
    allChannelIndex = allChannelIndex + n;
end




% % Append single channel data
% cfg = [];
% combined_data = ft_appenddata(cfg, singleChannelData{:});
% 
% 
% % Medial and lateral division
% cfg = [];
% cfg.channel = {'LFPx105' 'LFPx228' 'LFPx229' 'LFPx230' 'LFPx231'};
% medial_data = ft_selectdata(cfg, combined_data);
% 
% cfg = [];
% cfg.channel = {'LFPx106' 'LFPx107' 'LFPx132' 'LFPx133' 'LFPx134' 'LFPx57'};
% lateral_data = ft_selectdata(cfg, combined_data);


%% 3. Time-Frequency Analysis
% TMS
freq_tms = cell(1,nChannel);
freq_tms_blc = cell(1,nChannel);
for i = 1:nChannel
    cfg = [];
    cfg.output     = 'pow';
    cfg.method     = 'wavelet';
    cfg.width      = 7;
    cfg.foi        = 2:1:200;
    cfg.toi        = -0.5:0.005:1.5;
    cfg.keeptrials = 'yes';
    freq_tms{i} = ft_freqanalysis(cfg, singleChanneltms{i});
    
    % Baseline Correction
    cfg = [];
    cfg.baseline = [-0.5 -0.05];
    cfg.baselinetype = 'relative';
    freq_tms_blc{i} = ft_freqbaseline(cfg, freq_tms{i});
end

% Sham
freq_sham = cell(1,nChannel);
freq_sham_blc = cell(1,nChannel);
for i = 1:nChannel
    cfg = [];
    cfg.output     = 'pow';
    cfg.method     = 'wavelet';
    cfg.width      = 7;
    cfg.foi        = 2:1:200;
    cfg.toi        = -0.5:0.005:1.5;
    cfg.keeptrials = 'yes';
    freq_sham{i} = ft_freqanalysis(cfg, singleChannelsham{i});
    
    % Baseline Correction
    cfg = [];
    cfg.baseline = [-0.5 -0.05];
    cfg.baselinetype = 'relative';
    freq_sham_blc{i} = ft_freqbaseline(cfg, freq_sham{i});
end




%% 4.1 Multi-Channel Power Spectrum Plot 
% TMS
for i = 1:nChannel
    subplot(ceil(sqrt(nChannel)), ceil(sqrt(nChannel)), i);
    
    % Extract Data for the channel
    data = squeeze(mean(freq_tms_blc{i}.powspctrm(:, 1, :, :), 1));
    
    % Log-scale frequencies
    log_frequencies = log(freq_tms_blc{i}.freq);
    
    % Plot
    imagesc(freq_tms_blc{i}.time, log_frequencies, data);
    
    % Flip y-axis
    set(gca, 'YDir', 'normal');
    
    % Customize Axes
    yticks(log([2, 4, 10, 30, 70, 150]));
    yticklabels({'2', '4', '10', '30', '70', '150'});
    
    % Set color bar limits
    caxis([0 3]);
    
    title(['TMS, Channel: ' amygChannelAll{i}]);
end

colorbar;

% Sham
for i = 1:nChannel
    subplot(ceil(sqrt(nChannel)), ceil(sqrt(nChannel)), i);
    
    % Extract Data for the channel
    data = squeeze(mean(freq_sham_blc{i}.powspctrm(:, 1, :, :), 1));
    
    % Log-scale frequencies
    log_frequencies = log(freq_sham_blc{i}.freq);
    
    % Plot
    imagesc(freq_sham_blc{i}.time, log_frequencies, data);
    
    % Flip y-axis
    set(gca, 'YDir', 'normal');
    
    % Customize Axes
    yticks(log([2, 4, 10, 30, 70, 150]));
    yticklabels({'2', '4', '10', '30', '70', '150'});
    
    % Set color bar limits
    caxis([0 3]);
    
    title(['Sham, Channel: ' amygChannelAll{i}]);
end

colorbar;


%% 4.2 Average power spectrum across channels
% TMS
% average_powspctrm - a matrix of size [N_Frequencies x N_Times]
powspctrm_tms = cell(1, nChannel);
for i = 1:nChannel
    powspctrm_tms{i} = squeeze(nanmean(freq_tms_blc{i}.powspctrm, 1));
end

% Get the averaged matrix 
[freq, bin] = size(powspctrm_tms{1});
average_powspctrm_tms = NaN(freq, bin, nChannel); % use NaN to initialize as there might be NaN values in the matrices
for i = 1:nChannel
    average_powspctrm_tms(:, :, i) = powspctrm_tms{i};
end
average_powspctrm_tms = nanmean(average_powspctrm_tms, 3);

% Plot averaged power spectrum
% Log-scale frequencies
log_frequencies = log(freq_tms_blc{1}.freq);

figure;  % Create a new figure for the averaged power spectrum

imagesc(freq_tms_blc{1}.time, log_frequencies, average_powspctrm_tms);

% Flip y-axis
set(gca, 'YDir', 'normal');

% Customize Axes
yticks(log([2, 4, 10, 30, 70, 150]));
yticklabels({'2', '4', '10', '30', '70', '150'});

% Set color bar limits
caxis([0 3]);

% Add labels and title
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Averaged Power Spectrum Across Amygdala Channels - TMS');

% Add color bar
colorbar;


% Sham
powspctrm_sham = cell(1, nChannel);
for i = 1:nChannel
    powspctrm_sham{i} = squeeze(nanmean(freq_sham_blc{i}.powspctrm, 1));
end

% Get the averaged matrix 
[freq, bin] = size(powspctrm_sham{1});
average_powspctrm_sham = NaN(freq, bin, nChannel); % use NaN to initialize as there might be NaN values in the matrices
for i = 1:nChannel
    average_powspctrm_sham(:, :, i) = powspctrm_sham{i};
end
average_powspctrm_sham = nanmean(average_powspctrm_sham, 3);

% Plot averaged power spectrum
% Log-scale frequencies
log_frequencies = log(freq_sham_blc{1}.freq);

figure;  % Create a new figure for the averaged power spectrum

imagesc(freq_sham_blc{1}.time, log_frequencies, average_powspctrm_sham);

% Flip y-axis
set(gca, 'YDir', 'normal');

% Customize Axes
yticks(log([2, 4, 10, 30, 70, 150]));
yticklabels({'2', '4', '10', '30', '70', '150'});

% Set color bar limits
caxis([0 3]);

% Add labels and title
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Averaged Power Spectrum Across Amygdala Channels - Sham');

% Add color bar
colorbar;

