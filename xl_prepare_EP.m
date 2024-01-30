close all; clc;
clear;
rootpath = '/Users/xianqliu/Library/CloudStorage/OneDrive-UniversityofIowa/Documents - JingJiang Lab/JiangLab/4_Projects/2_DLPFC_AMY_Bella/1_iES_iEEG/Preprocessed';
cd(rootpath);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimsite = "L_DLPFC";
filter = 2;               % 0, 1, 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datapath = strcat(rootpath, "/", stimsite);
cd(datapath);

%%
% if stimsite == "L_DLPFC"
%     patient_list = [585];
%     targetchannel_list{1} = [10 105:109];
% 
% %     patient_list = [460 534 625 634];
% %     targetchannel_list{1} = [105 106 107];
% %     targetchannel_list{2} = [228 229 230 231];
% %     targetchannel_list{3} = [132 133 134];
% %     targetchannel_list{4} = [57];
% 
% elseif stimsite == "R_DLPFC"
%     patient_list = [376 429 532 561 567 610];
%     targetchannel_list{1} = [51:53];
%     targetchannel_list{2} = [120:124];
%     targetchannel_list{3} = [99:102];
%     targetchannel_list{4} = [189 244 245 247];
%     targetchannel_list{5} = [73:78];
%     targetchannel_list{6} = [92:95 149:156];
% 
% elseif stimsite == "L_Parietal"
%     patient_list = [416 430 625 634];
%     targetchannel_list{1} = [15:17,84];
%     targetchannel_list{2} = [8:14 43 44];
%     targetchannel_list{3} = [132 133 134];
%     targetchannel_list{4} = [57];
% 
% elseif stimsite == "R_Parietal"
%     patient_list = [430 515 561];
%     targetchannel_list{1} = [8:14 43 44];
%     targetchannel_list{2} = [13:19];
%     targetchannel_list{3} = [189,244,245,247];
% 
% elseif stimsite == "L_IFG"
%     patient_list = [403 423 561 585];
%     targetchannel_list{1} = [58 59];
%     targetchannel_list{2} = [79:85 115 116];
%     targetchannel_list{3} = [189 244 245 247];
%     targetchannel_list{4} = [10 105:109];
% 
% elseif stimsite == "R_IFG"
%     patient_list = [376 413 418 429 532 567 610];
%     targetchannel_list{1} = [51:53];
%     targetchannel_list{2} = [49:51];
%     targetchannel_list{3} = [149:151];
%     targetchannel_list{4} = [120:123];
%     targetchannel_list{5} = [99:101];
%     targetchannel_list{6} = [73:78];
%     targetchannel_list{7} = [92:95 149:151];
% 
% elseif stimsite == "IFG" % Hipp
%     patient_list = [403 423 538 561 585];
%     targetchannel_list{1} = [54:57 112 119:123];
%     targetchannel_list{2} = [111:114];
%     targetchannel_list{3} = [7 195 198];
%     targetchannel_list{4} = [113 245];
%     targetchannel_list{5} = [87 88 97 17:19];
% 
% elseif stimsite == "DLPFC" % hipp
%     patient_list = [460 585 625 634];
%     targetchannel_list{1} = [83:86 102 103];
%     targetchannel_list{2} = [17:19 87 88 97];
%     targetchannel_list{3} = [137 146];
%     targetchannel_list{4} = [106];
% end

% L_DLPFC
session_list = {'460-088','534-077','585-076','625-123', '625-101','634-155','634-082','625-096','634-084'};
targetchannel_list{1} = [105 106 107]; %460
targetchannel_list{2} = [228 229 230 231]; %534
targetchannel_list{3} = [10 105:109]; %585
targetchannel_list{4} = [132 133 134]; %625
targetchannel_list{5} = [132 133 134]; 
targetchannel_list{6} = [57]; %634
targetchannel_list{7} = [57];
targetchannel_list{8} = [132 133 134]; 
targetchannel_list{9} = [57]; %634


session_list = {'625-126','634-097','634-159'};
targetchannel_list{1} = [132 133 134]; 
targetchannel_list{2} = [57]; %634
targetchannel_list{3} = [57]; %634

%L_VLPFC
session_list = {'403-067','423-125','423-127','423-133', '561-124','585-056','585-072','585-077','585-184'};
targetchannel_list{1} = [58 59]; %403
targetchannel_list{2} = [79:85 115 116]; %423
targetchannel_list{3} = [79:85 115 116]; %423
targetchannel_list{4} = [79:85 115 116]; %423
targetchannel_list{5} = [189 244 245 247]; %561
targetchannel_list{6} = [10 105:109]; %585
targetchannel_list{7} = [10 105:109]; %585
targetchannel_list{8} = [10 105:109]; %585
targetchannel_list{9} = [10 105:109]; %585

% patient
for p = 1:length(patient_list)
    fprintf(['Processing Patient ' num2str(patient_list(p)) '\n']);
    xl_get_EPtable(patient_list(p), stimsite, targetchannel_list{p},filter);
end

% session
for s = 1:length(session_list)
    fprintf(['Processing Session ' session_list{s} '\n']);
    xl_get_EPtable_sessions(session_list{s}, stimsite, targetchannel_list{s},filter);
end













