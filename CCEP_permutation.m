clear;
clc;

%% SETUP
rootpath = '/Users/xianqliu/Library/CloudStorage/OneDrive-UniversityofIowa/1_TMS_iEEG/3_ProcessedData_Backup';
cd(rootpath);

ROI_list = {'L-DLPFC', 'R-DLPFC', 'L-VLPFC', 'R-VLPFC'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ROI = 1;
n_permutations = 10000;  % Number of permutations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timepoints = 9:450;

CCEP_amy_avg = readtable(['CCEP_amy_avg_' ROI_list{ROI} '.csv']);
CCEP_amy_avg = readtable(['CCEP_amy_avg_' ROI_list{ROI} '_b200.csv']);
CCEP_amy_avg = CCEP_amy_avg(:,2:11);


%% PERMUTATION TEST
% Both
perm_stats_B = permutation_CCEP(CCEP_amy_avg, timepoints, n_permutations);
perm_stats_B_me = permutation_CCEP(CCEP_amy_avg(CCEP_amy_avg.is_lateral == 0, :), timepoints, n_permutations);
perm_stats_B_la = permutation_CCEP(CCEP_amy_avg(CCEP_amy_avg.is_lateral == 1, :), timepoints, n_permutations);


% % Left
% perm_stats_L = permutation_CCEP(CCEP_amy_avg(strcmp(CCEP_amy_avg.DKT, 'Left-Amygdala'), :), timepoints, n_permutations);
% perm_stats_L_la = permutation_CCEP(CCEP_amy_avg(strcmp(CCEP_amy_avg.DKT, 'Left-Amygdala') & CCEP_amy_avg.is_lateral == 1, :), timepoints, n_permutations);
% 
% % Right
% perm_stats_R_me = permutation_CCEP(CCEP_amy_avg(strcmp(CCEP_amy_avg.DKT, 'Right-Amygdala') & CCEP_amy_avg.is_lateral == 0, :), timepoints, n_permutations);


%% 
for ROI = 1:length(ROI_list)
    disp(['>>> ' ROI_list{ROI} ' <<<']); fprintf('\n');
    
    %% read CCEP amy avg file
    CCEP_amy_avg = readtable(['CCEP_amy_avg_' ROI_list{ROI} '.csv']);
    CCEP_amy_avg = readtable(['CCEP_amy_avg_' ROI_list{ROI} '_b200.csv']);
    CCEP_amy_avg = CCEP_amy_avg(:,2:11);


    %% Both amy 
    % 1.1 Both
    disp('>>>> Both Amygdala'); fprintf('\n');
    perm_stats_B_onet = perm_CCEP_vs_baseline(CCEP_amy_avg, timepoints, n_permutations, 'one-sampled');
    perm_stats_B_pairedt = perm_CCEP_vs_baseline(CCEP_amy_avg, timepoints, n_permutations, 'paired');
    

    % 1.2 Both medial (= right medial)
    disp('>>>> Both Medial Amygdala'); fprintf('\n');
    perm_stats_B_me_onet = perm_CCEP_vs_baseline(CCEP_amy_avg(CCEP_amy_avg.is_lateral == 0, :), timepoints, n_permutations, 'one-sampled');
    perm_stats_B_me_pairedt = perm_CCEP_vs_baseline(CCEP_amy_avg(CCEP_amy_avg.is_lateral == 0, :), timepoints, n_permutations, 'paired');
    

    % 1.2 Both lateral (= left lateral)
    disp('>>>> Both Lateral Amygdala'); fprintf('\n');
    perm_stats_B_la_onet = perm_CCEP_vs_baseline(CCEP_amy_avg(CCEP_amy_avg.is_lateral == 1, :), timepoints, n_permutations, 'one-sampled');
    perm_stats_B_la_pairedt = perm_CCEP_vs_baseline(CCEP_amy_avg(CCEP_amy_avg.is_lateral == 1, :), timepoints, n_permutations, 'paired');


%     %% Left amy
%     % 2.1 Left
%     disp('>>>> Left Amygdala'); fprintf('\n');
%     perm_stats_L_onet = perm_CCEP_vs_baseline(CCEP_amy_avg(strcmp(CCEP_amy_avg.DKT, 'Left-Amygdala'), :), timepoints, n_permutations, 'one-sampled');
%     perm_stats_L_pairedt = perm_CCEP_vs_baseline(CCEP_amy_avg(strcmp(CCEP_amy_avg.DKT, 'Left-Amygdala'), :), timepoints, n_permutations, 'paired');
%     
% 
%     % 2.2 Left medial 
%     disp('>>>> Left Medial Amygdala'); fprintf('\n');
%     perm_stats_L_me_onet = perm_CCEP_vs_baseline(CCEP_amy_avg(strcmp(CCEP_amy_avg.DKT, 'Left-Amygdala') & CCEP_amy_avg.is_lateral == 0, :), timepoints, n_permutations, 'one-sampled');
%     perm_stats_L_me_pairedt = perm_CCEP_vs_baseline(CCEP_amy_avg(strcmp(CCEP_amy_avg.DKT, 'Left-Amygdala') & CCEP_amy_avg.is_lateral == 0, :), timepoints, n_permutations, 'paired');

%     % 2.3 Left lateral 
%     disp('>>>> Left Lateral Amygdala'); fprintf('\n');
%     perm_stats_L_la_onet = perm_CCEP_vs_baseline(CCEP_amy_avg(strcmp(CCEP_amy_avg.DKT, 'Left-Amygdala') & CCEP_amy_avg.is_lateral == 1, :), timepoints, n_permutations, 'one-sampled');
%     perm_stats_L_la_pairedt = perm_CCEP_vs_baseline(CCEP_amy_avg(strcmp(CCEP_amy_avg.DKT, 'Left-Amygdala') & CCEP_amy_avg.is_lateral == 1, :), timepoints, n_permutations, 'paired');

  
    %% Right amy
%     % 3.1 Right-–_
%     disp('>>>> Right Amygdala'); fprintf('\n');
%     perm_stats_R_onet = perm_CCEP_vs_baseline(CCEP_amy_avg(strcmp(CCEP_amy_avg.DKT, 'Right-Amygdala'), :), timepoints, n_permutations, 'one-sampled');
%     perm_stats_R_pairedt = perm_CCEP_vs_baseline(CCEP_amy_avg(strcmp(CCEP_amy_avg.DKT, 'Right-Amygdala'), :), timepoints, n_permutations, 'paired');

%     % 3.2 Right medial
%     disp('>>>> Right Medial Amygdala'); fprintf('\n');
%     perm_stats_R_me_onet = perm_CCEP_vs_baseline(CCEP_amy_avg(strcmp(CCEP_amy_avg.DKT, 'Right-Amygdala') & CCEP_amy_avg.is_lateral == 0, :), timepoints, n_permutations, 'one-sampled');
%     perm_stats_R_me_pairedt = perm_CCEP_vs_baseline(CCEP_amy_avg(strcmp(CCEP_amy_avg.DKT, 'Right-Amygdala') & CCEP_amy_avg.is_lateral == 0, :), timepoints, n_permutations, 'paired');
% 
%     % 3.3 Right lateral 
%     disp('>>>> Right Lateral Amygdala'); fprintf('\n');
%     perm_stats_R_la_onet = perm_CCEP_vs_baseline(CCEP_amy_avg(strcmp(CCEP_amy_avg.DKT, 'Right-Amygdala') & CCEP_amy_avg.is_lateral == 1, :), timepoints, n_permutations, 'one-sampled');
%     perm_stats_R_la_pairedt = perm_CCEP_vs_baseline(CCEP_amy_avg(strcmp(CCEP_amy_avg.DKT, 'Right-Amygdala') & CCEP_amy_avg.is_lateral == 1, :), timepoints, n_permutations, 'paired');


    save(fullfile(rootpath, ['CCEP_amy_ttest_rand' ROI_list{ROI} '.mat']), '-v7.3');
    save(fullfile(rootpath, ['CCEP_amy_ttest_rand' ROI_list{ROI} '_b200.mat']), '-v7.3');
    save(fullfile(rootpath, ['CCEP_amy_perm_t_' ROI_list{ROI} '.mat']), '-v7.3');
    save(fullfile(rootpath, ['CCEP_amy_perm_t_' ROI_list{ROI} '_b200.mat']), '-v7.3');
    
end


%% Paired t test with random baseline distrubution
% Both
disp('>>>> Both Amygdala'); fprintf('\n');
perm_stats_B = perm_CCEP_vs_baseline2(CCEP_amy_avg, timepoints, n_permutations);
disp('>>>> Both Medial Amygdala'); fprintf('\n');
perm_stats_B_me = perm_CCEP_vs_baseline2(CCEP_amy_avg(CCEP_amy_avg.is_lateral == 0, :), timepoints, n_permutations);
disp('>>>> Both Lateral Amygdala'); fprintf('\n');
perm_stats_B_la = perm_CCEP_vs_baseline2(CCEP_amy_avg(CCEP_amy_avg.is_lateral == 1, :), timepoints, n_permutations);

% % Left
% disp('>>>> Left Amygdala'); fprintf('\n');
% perm_stats_L = perm_CCEP_vs_baseline2(CCEP_amy_avg(strcmp(CCEP_amy_avg.DKT, 'Left-Amygdala'), :), timepoints, n_permutations);
% 
% 
% % Right
% disp('>>>> Right Medial Amygdala'); fprintf('\n');
% perm_stats_R_me = perm_CCEP_vs_baseline2(CCEP_amy_avg(strcmp(CCEP_amy_avg.DKT, 'Right-Amygdala') & CCEP_amy_avg.is_lateral == 0, :), timepoints, n_permutations);


%% Create the histogram
% figure;
histogram(perm_stats_B.null_t, 'BinWidth', 0.2);
xlim([-10, 10]);

% Plot the vertical line
hold on;
line([perm_stats_B.CCEP_t(150,1), perm_stats_B_onet.CCEP_t(150,1)], ylim, 'LineWidth', 2, 'Color', 'r');

% Calculate 95% confidence interval
CI_lower = mean(perm_stats_B_onet.null_t) - 1.96 * std(perm_stats_B_onet.null_t);
CI_upper = mean(perm_stats_B_onet.null_t) + 1.96 * std(perm_stats_B_onet.null_t);

line([CI_lower, CI_lower], ylim, 'LineWidth', 2, 'Color', 'g');
line([CI_upper, CI_upper], ylim, 'LineWidth', 2, 'Color', 'g');


% Label the graph
legend('Null t distribition', 'Single timepoint t Value', 'Lower 95% CI', 'Upper 95% CI');


%% kstest
%% 

[perm_stats_B.swtest_h, perm_stats_B.swtest_p] = swtest(perm_stats_B.null_t);
qqplot(perm_stats_B.null_t);
std = std(perm_stats_B.null_t);
%% 

% Generate some normally distributed data
data = normrnd(0,0.1,1000,1);

% Test for normality
[h,p] = jbtest(data);
