% 27/01/2026

% TODO
% 1. Set up paths and session sessions
% 2. Compute particulars files for each session
% 3. Extract chunks for each condition
% 4. Superimpose spectrum for conditions for a single channel



%% 1. Set up paths and session sessions

% basepath to session folders
basepath = fullfile(getenv('HOME'), 'Local Data', 'EXELU');

% get session folder sessions
sessions = setdiff({dir(basepath).name}, {'.', '..', '.DS_Store'});

% load all sessions
% sessions = cell(1, numel(sessions));
% for i = 1:numel(sessions)
%     name = sessions{i};
%     fprintf('Loading session: %s\n', name);
%     sessions{i} = loadLFPSession(basepath, name);
% end

clear session % clear previous session variable
idx = 1; % index of session to process
session = loadLFPSession(basepath, sessions{idx});
fs = str2double(session.settings.parameters.fieldPotentials.lfpSamplingRate.Text);



%% 2. Load particulars files for each session and/or compute bad channels if necessary

% load existing particulars file and update bad channels if necessary
particulars = loadParticulars(fullfile(basepath, sessions{idx}, 'session_particulars.txt'));

% compute bad channels if not already present
if ~isfield(particulars,'bad_channels') || isempty(particulars.bad_channels)
    segmentDuration = 0.1; % seconds
    badChans = badChannels(session, fs, segmentDuration);
    display(['Computed bad channels: ', mat2str(badChans)]);
    particulars.bad_channels = badChans;

    if ~isempty(badChans)
        fid = fopen(fullfile(basepath, sessions{idx}, 'session_particulars.txt'), 'a');
        fprintf(fid, 'bad_channels= %s\n', strjoin(string(badChans), ' '));
        fclose(fid);
    end
end

session.qc.badChannels = str2double(split(particulars.bad_channels));
session.qc.goodMask = true(size(session.data,1),1);
session.qc.goodMask(session.qc.badChannels) = false;
session.qc.badMask = ~session.qc.goodMask;

data = session.data(session.qc.goodMask, :); % use only good channels
tWin = [5 10]; % [pre post] seconds
duration = size(data, 2) / fs;



%% 3.1 Extract chunks for condition 1
% Static Grating (on and off samples)

condition = 1;

% extract triggers
events = session.events.analogEvents{1, condition};
events = reshape(events, 2, [])'; % [start end] per trial
triggers = events(:,1) ./ fs; % trial onsets (seconds)

% clean triggers to ensure enough data for pre and post windows
triggers = filterValidTriggers(triggers, tWin, duration);
triggers = triggers(valid);

% parameters for multitaper spectral estimation
movingWin = [0.2 0.2*0.2]; % 200 ms window, 20% step size
params = struct();
params.tapers = [2 3];
params.Fs     = fs;
params.fpass  = [10 200];
params.pad    = 1;

% compute triggered spectra
spec = computeTriggeredSpectra( ...
    data, fs, triggers, tWin, movingWin, params);

session.responses.spectra{condition} = spec;



%% 3.2 Extract chunks for condition 2
% Full Screen Flickering @ 40 Hz (sample full screen flickers on then off)

condition = 2;
pulses = session.events.analogEvents{1, condition} ./ fs;
gapThresh = 5; % seconds
triggers = extractTrialTriggers(pulses, gapThresh);
triggers = filterValidTriggers(triggers, tWin, duration);

spec = computeTriggeredSpectra( ...
    data, fs, triggers, tWin, movingWin, params);

session.responses.spectra{condition} = spec;

%% 3.3 Extract chunks for condition 3
% Phase Reversing Grating @ 40 Hz (sample at which frame one turns on and off—starts with white)

condition = 3;

% trigger extraction (implicit trials via pulse gaps) 
pulses = session.events.analogEvents{1, condition} ./ fs;  % seconds
gapThresh = 5;  % seconds
triggers = extractTrialTriggers(pulses, gapThresh);
triggers = filterValidTriggers(triggers, tWin, duration);

%spectra computation (generalised, unchanged) ----
spec = computeTriggeredSpectra( ...
    data, fs, triggers, tWin, movingWin, params);

session.responses.spectra{condition} = spec;

%% 4. Plot spectra for across conditions

S = session.responses.spectra{condition}.S_mu;   % freq × channel
f = session.responses.spectra{condition}.f;

figure; hold on

nCh = size(S,2);

% depth-aware colormap: deep (ch=1) dark, superficial light
cmap = flipud(parula(nCh));

for ch = 1:nCh
    plot(f, S(:,ch), ...
        'Color', cmap(ch,:), ...
        'LineWidth', 1);
end

set(gca, 'XScale', 'log', 'YScale', 'log')

% physiologically sensible limits
xlim([min(f(f>1)) 100])
% ylim(prctile(S(:), [5 95]))

xlabel('Frequency (Hz)')
ylabel('Power')
title('Mean Triggered Spectra Across Depth')

colormap(cmap)
cb = colorbar;
cb.Label.String = 'Channel (deep → superficial)';
cb.Ticks = [0 1];
cb.TickLabels = {'Deep', 'Superficial'};

box off
grid on



% condition 3: Full Screen Flickering @ 40 Hz (sample full screen flickers on then off)
% condition 4: Phase Reversing Grating @ 40 Hz (sample at which frame two turns on and off —starts with black)
% condition 5: Running speed (sample at which pulse triggered; 600 pulses per rotation)
