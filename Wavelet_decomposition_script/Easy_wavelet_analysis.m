% Generic wavelet analysis
%% Main inputs
function Easy_wavelet_analysis(dataOrig,sampling_rate,timestochop,preStimSecs,postStimSecs,baselineTime,MarkTimes)

% dataOrig = ch x samples matrix;
% sampling_rate = srate;
% timestochop = times of events;
% preStimSecs = pre stim time in seconds (negative if before)
% postStimSecs = post stim time in seconds
% baselineTime = [startTimeInMs endTimeInMs];
% MarkTimes = e.g. [-3.25 -3 0 5]*1000; % this would plot a marker at the x zero crossing. add any other times in secs to also plot these
chansToAnalyse = 1:size(dataOrig,1); % as an index (not relating to LFPx)
LFPxchansToAnalyse = []; % relating to LFPx channels - nice for subset analysis
baselineToRemove = baselineTime; % in ms; removes this baseline prior to wavelet analysis
timesout = 200; % temporal resolution of wavelet - divide entire duration (secs) of epoch by this to get resolution in seconds
alpha_val = 0.01; % significance threshold

data = epoch( dataOrig, timestochop, [preStimSecs postStimSecs] ,'srate', sampling_rate);

%% Set up for wavelet analysis
eeglab
ElecNums = 1:size(data,1); % assuming you want to analyse all electrodes
NumTrials = size(data,3);
padratio = 4;
waveFreq = [2 150];
% plotting details:
maxersp = 4;
itcplot = 'off';

ersplim = [-maxersp maxersp];
trialBaseStatus = 'full'; % either 'on', 'off' or 'full'

%% don't clean or analyse zeroed channels
chansToAnalyse(mean(mean(data,3),2)==0) = [];
data(mean(mean(data,3),2)==0,:,:) = [];

%% load data into EEG lab and clean
EEG = pop_importdata( 'dataformat', 'array', 'data', 'data', 'setname', 'AllEpochs', 'srate',sampling_rate, 'pnts',0, 'xmin',0, 'nbchan',0);
EEG = eeg_checkset( EEG );
EEG = pop_importepoch( EEG, ones(1,NumTrials),{'AllEpochs'}, 'latencyfields',{'AllEpochs'}, 'timeunit',1, 'headerlines',0);
EEG = eeg_checkset( EEG );

EEG.xmin = EEG.xmin+preStimSecs;
EEG.xmax = EEG.xmax+postStimSecs;
% you may or may not want these next two lines, depending on whether you
% want to remove a baseline period prior to wavelet analysis
EEG = pop_rmbase( EEG, baselineToRemove);
EEG = eeg_checkset( EEG );

% go through cleaning
eventIDX = 1:NumTrials;
rejIndx = {};
for reject_cycle = 1:5
    EEG = pop_jointprob(EEG,1,chansToAnalyse,5,5,0,0);
    EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
    rejIndx{1}{reject_cycle} = EEG.reject.rejglobal;
    EEG = pop_rejepoch( EEG, find(EEG.reject.rejglobal), 0);
    EEG = eeg_checkset( EEG );
end
for reject_cycle = 1:5
    EEG = pop_rejkurt(EEG,1,chansToAnalyse,7,7,0,0);
    EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
    rejIndx{2}{reject_cycle} = EEG.reject.rejglobal;
    EEG = pop_rejepoch( EEG, find(EEG.reject.rejglobal), 0);
    EEG = eeg_checkset( EEG );
end

% Store the new event indices (removing the excluded trials) - in case you
% want them later
newIndices = eventIDX;
for stupidRejection = 1:5
    newIndices(rejIndx{1}{stupidRejection}) = [];
end
for stupidRejection = 1:5
    newIndices(rejIndx{2}{stupidRejection}) = [];
end

% specify particular LFPx channels (finds index for these)
if ~isempty(LFPxchansToAnalyse)
    indxchansToAnalyse = [];
    for allChans = 1:numel(LFPxchansToAnalyse)
        indxchansToAnalyse(allChans) = find(contains(chanNames,['LFPx' num2str(LFPxchansToAnalyse(allChans)) '_']));
    end
    chansToAnalyse = indxchansToAnalyse;
end
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
        1, elecidx, [EEG.xmin EEG.xmax]*1000,[3 0.5],'freqs',[waveFreq(1) waveFreq(2)],'padratio', padratio, ...
        'plotphase', 'off', 'timesout',timesout,'naccu', 200, 'alpha',alpha_val,'baseboot',1,'rmerp','off','freqscale','log', ...
        'plotersp','off', 'plotitc',itcplot,'marktimes',MarkTimes,'erspmax', maxersp,'trialbase',trialBaseStatus,'baseline',baselineToRemove);

    % plot using my function (it's a little nicer than EEGLAB)
    do_EEGLABPlot(ersp,erspboot,times,freqs,ersplim,MarkTimes,...
        ['Channel: ' num2str(contact)],1)

    clear in ersp itc powbase erspboot itcboot tfdata

end