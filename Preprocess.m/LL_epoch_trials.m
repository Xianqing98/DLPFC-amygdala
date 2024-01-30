%% START UP FILE FOR TMS-iEEG ANALYSIS - Epoching trials
% DANNY HUANG, JEFF WANG, COREY KELLER
% 3/23/2020

% REQUIRED TOOLBOXES
% FIELDTRIP
% NEURALYNX V6.0
% DANNY's IEP SCRIPTS (JEFF I SHARED THIS WITH YOU)

close all
clear
gcp % Start a parpool if it doesn't exist

%% Set up directories
dataDir =  'C:\Users\xjl19\Desktop\TMS-iEEG\xl_Data_Processed';
queueFName = 'C:\Users\xjl19\Desktop\TMS-iEEG\Processed_xlsx\xl_Epoch_Queue.xlsx';
saveDir = 'C:\Users\xjl19\Desktop\TMS-iEEG\xl_Data_Processed\Comparisons_TMS_sham_allSites';


%% ESTABLISH DIRECTORY AND ADD PATH
cd(fileparts(which(mfilename))); %Navigate to the script's location
addpath 'C:\Users\xjl19\Desktop\Matlab\fieldtrip-20230206'
addpath(genpath('.\Jeff_Scripts'), ...
    genpath('.\dannyIEPScripts'), ...
    genpath('.\dannyGenScripts'), ...
    genpath('.\MatlabImportExport_v6.0.0'), ...
    genpath('.\util'));
ft_defaults;

sham = true;

%% Import file names to be processed
queue = readtable(queueFName);

for i = 1:height(queue)
    close all
    disp("-------------------------------------" + newline + ...
        "Working on Patient " + num2str(queue{i,'Patient'}) + ...
        ", Condition " + queue{i,'Folder_Name'}{1})
    if queue{i,'Completed'} == 1
        disp('Patient already processed. Moving on.')
        continue;
    end
    tic;
    %Save Directory for everything
    dataDirSess = fullfile(dataDir, num2str(queue{i,'Patient'}), queue{i,'Folder_Name'}{1});


    %% IMPORT DATA
    tms_file_stem = strcat(num2str(queue{i,'Patient'}),'_TMS_',queue{i,'Folder_Name'}{1});
    sham_file_stem = strcat(num2str(queue{i,'Patient'}),'_Sham_',queue{i,'Folder_Name'}{1});
    tmsFName = fullfile(dataDirSess,'4-Epoched',strcat(tms_file_stem,'_4A-CleanEpoched.mat'));
    shamFName = fullfile(dataDirSess,'4-Epoched',strcat(sham_file_stem,'_4A-CleanEpoched.mat'));
    
    load(tmsFName);
    ftData_epoch_tms = ftData;
    
    load(shamFName);
    ftData_epoch_sham = ftData;
    
    disp('File Load Complete')
   
    
    %% Reselect Channels
    % Select only for channels that are connected
    validChannels = strtrim(cellstr(num2str(channelInfo.Channel)));
    validChannels = strcat('LFPx',validChannels);
    
    % Now select subset
    cfg = [];
    cfg.channel = validChannels;
    ftData = ft_selectdata(cfg, ftData);
    
    %% Apply Baseline Correction
    cfg=[];
    cfg.demean='yes';
    cfg.baselinewindow=[-.25, -.05];
    
    ftData_epoch_tms = ft_preprocessing(cfg,ftData_epoch_tms);
    if sham
        ftData_epoch_sham = ft_preprocessing(cfg,ftData_epoch_sham);
    end

    
    %% Timelock Analysis
    cfg = [];
    close all
    channel_names = cellfun(@(x) sscanf(x,'LFPx%f'),ftData.label);
    [~,sorted_idx] = sortrows(channel_names);
    
    ftData_timelocked_tms = ft_timelockanalysis(cfg, ftData_epoch_tms);
    ftData_timelocked_sham = ft_timelockanalysis(cfg, ftData_epoch_sham);


    %% Detrend
    if queue{i,'Detrend'} == 1
        mask = (ftData_epoch_tms.time{1} > 0.025) & (ftData_epoch_tms.time{1} < 0.5); % Limit fit to 25 to 500 ms post trigger
        time_fit = ftData_epoch_tms.time{1}(mask);
        % models_tms = struct([]);
        expFun = @(a,t) a(1) * exp(-abs(a(2)) * (t-a(5))) + a(3) * exp(-abs(a(4)) * (t-a(5)));
        start_fit = [200,14,200,10,0.03];
        parfevalOnAll(gcp(), @warning, 0, 'off');
        
        % progressbar = waitbar(0,'Detrending TMS Condition...');
        disp('Detrending TMS Condition...')
        trial_data = ftData_epoch_tms.trial;
        ftData_timelocked_tms.var(ftData_timelocked_tms.var == 0) = 1;
        all_weights = 1 ./ ftData_timelocked_tms.var;
        N_channels = length(ftData_epoch_tms.label);
        
        parfor tr = 1:length(trial_data)
            
            % waitbar((length(ftData_epoch_tms.label) - ch) / length(ftData_epoch_tms.label),progressbar,'Detrending TMS Condition...');
            for ch = 1:N_channels 
                weights = all_weights(ch,mask); % Compute weights as inverse of variance
                model_linear = fit(time_fit.',trial_data{tr}(ch,mask).','poly1','Weights',weights);
                try
                    model_exp = fitnlm(time_fit.',trial_data{tr}(ch,mask).',expFun,start_fit,'Weight',weights);
                catch exception
                    model_exp = [];
                end
                % Calculate Akaike Information Criterion (AIC)
                linear_aic = length(time_fit) * log(sum((model_linear(time_fit).' - trial_data{tr}(ch,mask)).^2 .* weights)) ...
                    + 2;
                
                if isempty(model_exp)
                    exp_aic = Inf;
                else
                    exp_aic = length(time_fit) * log(sum((model_exp.feval(time_fit) - trial_data{tr}(ch,mask)).^2 .* weights)) ...
                        + 2;
                end
                
                if linear_aic < exp_aic
                    model_model_choice = 0;
                    trial_data{tr}(ch,mask) = trial_data{tr}(ch,mask) - model_linear(time_fit).';
                else
                    model_model_choice = 1;
                    trial_data{tr}(ch,mask) = trial_data{tr}(ch,mask) - model_exp.feval(time_fit);
                end
            end
        end
        
        ftData_epoch_tms.trial = trial_data;
        
        
        %progressbar = waitbar(0,'Detrending Sham Condition...');
        disp('Detrending Sham Condition...')
        trial_data = ftData_epoch_sham.trial;
        ftData_timelocked_sham.var(ftData_timelocked_sham.var == 0) = 1;
        all_weights = 1 ./ ftData_timelocked_sham.var;
        N_channels = length(ftData_epoch_sham.label);
        
        parfor tr = 1:length(trial_data)
            % Compute weights as inverse of variance
            % waitbar((length(ftData_epoch_sham.label) - ch) / length(ftData_epoch_sham.label),progressbar,'Detrending Sham Condition...');
            for ch = 1:N_channels
                weights = all_weights(ch,mask);
                
                model_linear = fit(time_fit.',trial_data{tr}(ch,mask).','poly1','Weights',weights);
                try
                    model_exp = fitnlm(time_fit.',trial_data{tr}(ch,mask).',expFun,start_fit,'Weight',weights);
                catch exception
                    model_exp = [];
                end
                % Calculate Akaike Information Criterion (AIC)
                linear_aic = length(time_fit) * log(sum((model_linear(time_fit).' - trial_data{tr}(ch,mask)).^2 .* weights)) ...
                    + 2;
                
                if isempty(model_exp)
                    exp_aic = Inf;
                else
                    exp_aic = length(time_fit) * log(sum((model_exp.feval(time_fit) - trial_data{tr}(ch,mask)).^2 .* weights)) ...
                        + 2;
                end
                
                if linear_aic < exp_aic
                    model_model_choice = 0;
                    trial_data{tr}(ch,mask) = trial_data{tr}(ch,mask) - model_linear(time_fit).';
                else
                    model_model_choice = 1;
                    trial_data{tr}(ch,mask) = trial_data{tr}(ch,mask) - model_exp.feval(time_fit);
                end
            end
        end
        ftData_epoch_sham.trial = trial_data;
    end
    
    %% Rerun Timelock Analysis after Detrending
    cfg = [];
    
    ftData_timelocked_tms = ft_timelockanalysis(cfg, ftData_epoch_tms);
    ftData_timelocked_sham = ft_timelockanalysis(cfg, ftData_epoch_sham);

    
    %% Compute averages and SE's
    avg_clean = ftData_timelocked_tms.avg;
    avg_sham = ftData_timelocked_sham.avg;
    n_trials = length(ftData_epoch_tms.time);
    
    se_clean = sqrt(ftData_timelocked_tms.var / n_trials);
    se_sham = sqrt(ftData_timelocked_sham.var / n_trials);
    
    %% Calculate test-statistic
    t_stats = (avg_clean - avg_sham);
    t_stats = t_stats ./ sqrt(se_clean.^2 + se_sham.^2);
    dof = n_trials - 1;
    p_values = 2*tcdf(-abs(t_stats),dof);
   
    
    %% Save Results
    disp('Saving Epoched Data')
    mkdir(saveDir);
    
    save_file_stem = strcat(num2str(queue{i,'Patient'}),'_Comparison_',queue{i,'Folder_Name'}{1});
    if queue{i,'Detrend'} == 1 
        saveFName = fullfile(saveDir,strcat(save_file_stem,'_4A-CleanEpoched_detrend.mat'));
    else
        saveFName = fullfile(saveDir,strcat(save_file_stem,'_4A-CleanEpoched.mat'));
    end
    
    save(saveFName, 'sessionDir', 'ftData_epoch_tms', 'ftData_epoch_sham', 'ftData_timelocked_tms', 'ftData_timelocked_sham', 'channelInfo', 't_stats', 'p_values', '-v7.3')
    toc;


    %% plot data
    saveFile = fullfile(dataDir, num2str(queue{i,'Patient'}), queue{i,'Folder_Name'}{1},'5-Figures');
    mkdir(saveFile)

    ylabel_num = [-0.5,0,0.5,1,1.5,2];
    ylabel_name = {1,length(ylabel_num)};
    for y = 1:length(ylabel_num)
        ylabel_name{1,y} = [num2str(ylabel_num(y))];
    end
    ylabel_num = ylabel_num*1000+1000;

    disp('Ploting Epoched Data')
    for cond = ["TMS" "Sham"]
        if cond == "TMS"
            ftData = ftData_epoch_tms;
        elseif cond == "Sham"
            ftData = ftData_epoch_sham;
        end

        if ~isempty(queue{i,'Example_Channel'}{1})
            split_channel = split(queue{i,'Example_Channel'}{1},';');
            
            for example_i = 1:size(split_channel,1)
                for chn_i = 1:length(ftData.label)
                    if strcmp(ftData.label{chn_i,1},split_channel{example_i,1})
                        if queue{i,'Detrend'} == 1 
                            savename = fullfile(saveFile,strcat(num2str(queue{i,'Patient'}),'_',cond,'_', ...
                                queue{i,'Folder_Name'}{1},'_5-Detrended_',split_channel{example_i,1},'detrend.png'));
                        else
                            savename = fullfile(saveFile,strcat(num2str(queue{i,'Patient'}),'_',cond,'_', ...
                                queue{i,'Folder_Name'}{1},'_5-Detrended_',split_channel{example_i,1},'.png'));
                        end
                        disp(['Saving to: ' savename]); % Diagnostic print statement
                        
                        figure;
                        set(gcf, 'color', 'white');
                        set(gcf, 'Position', [100, 100, 3200, 800]);  % Set the figure size

                        for trial = 1:length(ftData.trial)
                            subplot(5,10,trial)
                            hold on
                            plot(ftData.trial{1,trial}(chn_i,:))
                            hold on
        
                            ylim([-100,100])
                            xticks(ylabel_num)
                            xticklabels(ylabel_name);
                            plot([990,990],[-100,100],'color',[0.5,0.5,0.5],'linewidth',1)
                            plot([1025,1025],[-100,100],'color',[0.5,0.5,0.5],'linewidth',1)
                            plot([0,3000],[0,0],'color',[0.5,0.5,0.5],'linewidth',0.5)
                        end
                        title([ftData.label{chn_i}]);
                        saveas(gcf, savename);
                        close(gcf); % Close the figure after saving
                    end
                end
            end
        end
    end   
end



