clear;
clc;

%% Set up dataOrig
rootpath = '/Users/xianqliu/Library/CloudStorage/OneDrive-UniversityofIowa/1_TMS_iEEG/3_ProcessedData_Backup/TimeFrequency';
cd(rootpath);
addpath /Users/xianqliu/Documents/MATLAB/eeglab2023.1
addpath(genpath('/Users/xianqliu/Documents/MATLAB/eeglab2023.1'));

CCEP_amy_Orig.LDLPFC = readtable(['CCEP_amy_ChlxTimes_L-DLPFC.csv']);
ChannelInfo_LDLPFC = CCEP_amy_Orig.LDLPFC(:, 2:5);
CCEP_amy_Orig.LDLPFC = CCEP_amy_Orig.LDLPFC{:, 6:150005};



%% Set up data
dataOrig = CCEP_amy_Orig.LDLPFC; % ch x samples matrix;
sampling_rate = 1000; % srate;

timestochop = 1000:2500:1000 + (60 - 1) * 2500; % times of events;
timestochop = timestochop/1000;
MarkTimes = [-8 0 8]; % e.g. [-3.25 -3 0 5]*1000; % this would plot a marker at the x zero crossing. add any other times in secs to also plot these
chansToAnalyse = 1:size(dataOrig,1); % as an index (not relating to LFPx)
LFPxchansToAnalyse = []; % relating to LFPx channels - nice for subset analysis

baselineToRemove = [-500 -10]; % % baselineTime = [startTimeInMs endTimeInMs]; % in ms; removes this baseline prior to wavelet analysis 
timesout = 200; % temporal resolution of wavelet - divide entire duration (secs) of epoch by this to get resolution in seconds
alpha_val = 0.05; % significance threshold
preStimSecs = -1; % pre stim time in seconds (negative if before)
postStimSecs = 1.5; % post stim time in seconds

data = epoch( dataOrig, timestochop, [preStimSecs postStimSecs] ,'srate', sampling_rate);

%% Set up for wavelet analysis
eeglab
ElecNums = 1:size(data,1); % assuming you want to analyse all electrodes
NumTrials = size(data,3);
padratio = 4;
waveFreq = [2 150]; % frequency range
% plotting details:
maxersp = 4;
itcplot = 'off';

ersplim = [-maxersp maxersp];
trialBaseStatus = 'full'; % either 'on', 'off' or 'full'

%% don't clean or analyse zeroed channels
chansToAnalyse(mean(mean(data,3),2)==0) = [];
data(mean(mean(data,3),2)==0,:,:) = [];

%% load data into EEG lab and clean
EEG = pop_importdata( 'dataformat', 'array', 'data', 'data', 'setname', 'AllEpochs', 'srate',sampling_rate, 'pnts',0, 'xmin',preStimSecs, 'nbchan',0);
EEG = eeg_checkset( EEG );
EEG = pop_importepoch( EEG, ones(1,NumTrials),{'AllEpochs'}, 'latencyfields',{'AllEpochs'}, 'timeunit',1, 'headerlines',0);
EEG = eeg_checkset( EEG );

% EEG.xmin = EEG.xmin+preStimSecs;
% EEG.xmax = EEG.xmax+postStimSecs;

% you may or may not want these next two lines, depending on whether you
% want to remove a baseline period prior to wavelet analysis
EEG = pop_rmbase( EEG, baselineToRemove);
EEG = eeg_checkset( EEG );


%% Each channel, cycle through and perform wavelet analysis
for contact = chansToAnalyse

    subplot(4,4,contact)
    set(gcf,'units','normalized','outerposition',[0 0 1 1])

    elec = ElecNums(contact);
    elecidx = find(elec==ElecNums);

    currData = EEG;

    fprintf(['Electrode number(of ' num2str(length(chansToAnalyse)) '): %d\n'], elecidx)
    fprintf('\n\n');

    [ersp,itc,powbase,times,freqs,erspboot,itcboot,tfdata] = pop_newtimef(currData, ...
        1, elecidx, [EEG.xmin EEG.xmax]*1000,[2 0.7],'freqs',[waveFreq(1) waveFreq(2)],'padratio', padratio, ...
        'plotphase', 'off', 'timesout',timesout,'naccu', 200, 'alpha',alpha_val,'baseboot',1,'rmerp','off','freqscale','log', ...
        'plotersp','off', 'plotitc',itcplot,'marktimes',MarkTimes,'erspmax', maxersp,'trialbase',trialBaseStatus,'baseline',baselineToRemove);

    % plot using my function (it's a little nicer than EEGLAB)
    do_EEGLABPlot(ersp,erspboot,times,freqs,ersplim,MarkTimes,...
        ['Channel: ' num2str(contact)],1)
    
    clear in ersp itc powbase erspboot itcboot tfdata

end