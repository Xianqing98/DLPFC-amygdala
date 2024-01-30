% cluster-permutaion TMS L-DLPFC vs R-Parietal
clear all;

ies_DLPFC = readtable('/Users/xianqliu/Library/CloudStorage/OneDrive-UniversityofIowa/1_TMS_iEEG/3_ProcessedData_Backup/site_specificity_amy/TEP/archived/TEP_amy_ies_LDLPFC.csv');
ies_Parietal = readtable('/Users/xianqliu/Library/CloudStorage/OneDrive-UniversityofIowa/1_TMS_iEEG/3_ProcessedData_Backup/site_specificity_amy/TEP/archived/TEP_amy_ies_RParietal.csv');

condition = ies_DLPFC.Time > 25;  
ies_DLPFC = ies_DLPFC(condition, :);

condition = ies_Parietal.Time > 25;  
ies_Parietal = ies_Parietal(condition, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
laterality = 'R';   % B, L, R;
subdivision = 1;    % 2(whole), 0(medial), 1(lateral)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if laterality == 'L'
    condition = ismember(ies_DLPFC.DKT, ['Left-Amygdala']);
    ies_DLPFC = ies_DLPFC(condition, :);
    condition = ismember(ies_Parietal.DKT, ['Left-Amygdala']);
    ies_Parietal = ies_Parietal(condition, :);
    if subdivision == 0 % medial
        condition = ismember(ies_DLPFC.is_lateral, [0]);
        ies_DLPFC = ies_DLPFC(condition, :);
        condition = ismember(ies_Parietal.is_lateral, [0]);
        ies_Parietal = ies_Parietal(condition, :);
    elseif subdivision == 1 % lateral
        condition = ismember(ies_DLPFC.is_lateral, [1]);
        ies_DLPFC = ies_DLPFC(condition, :);
        condition = ismember(ies_Parietal.is_lateral, [1]);
        ies_Parietal = ies_Parietal(condition, :);
    end
elseif laterality == 'R'
    condition = ismember(ies_DLPFC.DKT, ['Right-Amygdala']);
    ies_DLPFC = ies_DLPFC(condition, :);
    condition = ismember(ies_Parietal.DKT, ['Right-Amygdala']);
    ies_Parietal = ies_Parietal(condition, :);
    if subdivision == 0 % medial
        condition = ismember(ies_DLPFC.is_lateral, [0]);
        ies_DLPFC = ies_DLPFC(condition, :);
        condition = ismember(ies_Parietal.is_lateral, [0]);
        ies_Parietal = ies_Parietal(condition, :);
    elseif subdivision == 1 % lateral
        condition = ismember(ies_DLPFC.is_lateral, [1]);
        ies_DLPFC = ies_DLPFC(condition, :);
        condition = ismember(ies_Parietal.is_lateral, [1]);
        ies_Parietal = ies_Parietal(condition, :);
    end
else
    if subdivision == 0 % medial
        condition = ismember(ies_DLPFC.is_lateral, [0]);
        ies_DLPFC = ies_DLPFC(condition, :);
        condition = ismember(ies_Parietal.is_lateral, [0]);
        ies_Parietal = ies_Parietal(condition, :);
    elseif subdivision == 1 % lateral
        condition = ismember(ies_DLPFC.is_lateral, [1]);
        ies_DLPFC = ies_DLPFC(condition, :);
        condition = ismember(ies_Parietal.is_lateral, [1]);
        ies_Parietal = ies_Parietal(condition, :);
    end
end

channel_DLPFC = unique(ies_DLPFC.Channel);
channel_Parietal = unique(ies_Parietal.Channel);


target_data = zeros(height(channel_DLPFC), 475);
for subject = 1:height(channel_DLPFC)
    condition = ismember(ies_DLPFC.Channel, channel_DLPFC(subject));
    data_temp = ies_DLPFC(condition, :);
    target_data(subject,:) = data_temp.iTEP_norm;
end

control_data = zeros(height(channel_Parietal), 475);
for subject = 1:height(channel_Parietal)
    condition = ismember(ies_Parietal.Channel, channel_Parietal(subject));
    data_temp = ies_Parietal(condition, :);
    control_data(subject,:) = data_temp.iTEP_norm;
end


%% cluster-based statistics

target_data;  % subject/electrodes * frequency/time
control_data;

time_num = size(target_data,2);

ttest_result = zeros(time_num,3);
for time_i = 1:time_num
    [h,p,ci,stats] = ttest2(target_data(:,time_i),control_data(:,time_i));
    ttest_result(time_i,1) = h;
    ttest_result(time_i,2) = p;
    ttest_result(time_i,3) = stats.tstat;
end

output_clusters = U_clusters_lzr(ttest_result(:,3),ttest_result(:,2),0.05);  % input: all statistics, p, threshold for cluster
cluster_real_raw = output_clusters;

temp_data = zeros(size(output_clusters,2),3);
for temp_i = 1:size(output_clusters,2)
    temp_data(temp_i,1) = sum(output_clusters{1,temp_i}.zs);  % cluster feature1: sum of statistics, most common
    temp_data(temp_i,2) = max(output_clusters{1,temp_i}.zs);  % cluster feature2: max of statistics
    temp_data(temp_i,3) = mean(output_clusters{1,temp_i}.zs); % cluster feature3: mean of statistics
end
cluster_real = temp_data;


%%
perm_num = 1000;   % >= 1000

tic
cluster_perm = zeros(perm_num,3);

for perm_i = 1:perm_num
    
    % shuffling method 1: channel level, shuffling the group labels among all sub/ch
    target_num = size(target_data,1);
    control_num = size(control_data,1);
    
    raw_data = [target_data; control_data];
    
    rand_info = randperm(target_num+control_num);
    temp_data1 = raw_data(find(rand_info<=target_num),:);
    temp_data2 = raw_data(find(rand_info>target_num),:);
    
    
    % or!!!!!
    % shuffling method 2: trial level, shuffling the trial labels for each sub/ch
    % for example: tms vs. sham 
    
    
    ttest_result = zeros(time_num,3);
    for time_i = 1:time_num
        [h,p,ci,stats] = ttest2(temp_data1(:,time_i),temp_data2(:,time_i));
        ttest_result(time_i,1) = h;
        ttest_result(time_i,2) = p;
        ttest_result(time_i,3) = stats.tstat;
    end
    
    output_clusters = U_clusters_lzr(ttest_result(:,3),ttest_result(:,2),0.05);
    
    if ~isempty(output_clusters)
        temp_data = zeros(size(output_clusters,2),3);
        for temp_i = 1:size(output_clusters,2)
            temp_data(temp_i,1) = sum(output_clusters{1,temp_i}.zs);
            temp_data(temp_i,2) = max(output_clusters{1,temp_i}.zs);
            temp_data(temp_i,3) = mean(output_clusters{1,temp_i}.zs);
        end
        if size(output_clusters,2) > 1
            cluster_perm(perm_i,:) = max(temp_data);  % only save the maximum of cluster features
        else
            cluster_perm(perm_i,:) = temp_data;
        end
    end
end
toc

%%
for i = 1:height(cluster_real)
    % Sort the distribution data
    sortedData = sort(cluster_perm(:,1));

    % Find the position of the cluster_real(i,1) in sorted data
    position = find(sortedData >= cluster_real(i,1), 1);

    % Calculate the percentile rank
    if isempty(position)
        percentile(i,1) = 100;  % If the value is greater than all elements
    else
        percentile(i,1) = (position / length(sortedData)) * 100;
    end
end









