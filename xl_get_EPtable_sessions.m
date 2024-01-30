function EPtable = xl_get_EPtable_sessions(session, stimsite, targetchannel, filter)


%% Import data
saveDir = strcat('/Users/xianqliu/Library/CloudStorage/OneDrive-UniversityofIowa/1_TMS_iEEG/3_ProcessedData_Backup/allPatients_iEEGdata_csv/CCEP_csv_', ...
    stimsite);
mkdir(saveDir);

if filter == 0
    load(strcat(session, "_", stimsite, "_0_5Hz_4B-UnfilterEpoched.mat"));
    filetype = '_unfilter.csv';
elseif filter == 1
    load(strcat(session, "_", stimsite, "_0_5Hz_4A-CleanEpoched.mat"));
    filetype = '_1hp.csv';
elseif filter == 2 
    load(strcat(session, "_", stimsite, "_0_5Hz_4A-CleanEpoched.mat"));
    filetype = '_2hp.csv';
end

patient = session(1:3);


%% Creating EP table
% EPtable = {'Patient','TrialNumber','Channel','Time','Amplitude'};  
rowIndex = 2; 
nTrials = length(ftData.trial);

t_lim = [-0.5, 1]; 
t_idx = find((ftData.time{1,1}(1,:)>=t_lim(1))&(ftData.time{1,1}(1,:)<=t_lim(2))); 

chl_idx = [];
for i = 1:height(ftData.label)
    if ismember(str2double(ftData.label{i,1}(5:end)),targetchannel)
        chl_idx = [chl_idx, i];
    end
end

% for chl = 1:length(chl_idx)
%     for epc = 1:nTrials
%         subplot(10,7,epc)
%         hold on
%         title([num2str(ftData.label{chl_idx(chl),1}) ' trial ' num2str(epc)], 'Interpreter', 'none')
%         plot(ftData.time{epc},ftData.trial{epc}(chl_idx(chl),:));
%         ylim([-80,80])
%         hold off
%     end
%     saveFile = fullfile('C:\Users\xjl19\Desktop\TMS-iEEG\es-TT_Processed\epochs_plot', ...
%         strcat(num2str(patient),'_',num2str(ftData.label{chl_idx(chl),1}),'.png'));
%     set(gcf, 'PaperPosition', [0 0 12 10])
%     saveas(gcf,saveFile)
%     close
% end

EPtable = cell(nTrials*length(chl_idx)*length(t_idx), 5);
for trl = 1:nTrials
    fprintf(['>> Trial ' num2str(trl) ' <<\n']);
    for chl = 1:length(chl_idx) 
        for t = 1:length(t_idx)  % Timewindow loop
            EPtable{rowIndex,1} = patient;
            EPtable{rowIndex,2} = session;
            EPtable{rowIndex,3} = trl;
            EPtable{rowIndex,4} = ftData.label{chl_idx(chl),1};
            EPtable{rowIndex,5} = ftData.time{1,trl}(1,t_idx(t))*1000;
            EPtable{rowIndex,6} = ftData.trial{1,trl}(chl_idx(chl),t_idx(t));
            rowIndex = rowIndex + 1;
        end
    end
end

fprintf(['Saving Patient ' num2str(patient) ' <<\n']);
writecell(EPtable, strcat(saveDir, '/sessions/', session, '_EP_', stimsite, filetype));



















