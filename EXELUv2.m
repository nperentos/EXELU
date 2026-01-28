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



%% 2. Compute particulars files for each session

% identify bad channels
segmentDuration = 100; % ms
badChans = badChannels(session, fs, segmentDuration);

% load existing particulars file and update bad channels if necessary
particulars = loadParticulars(fullfile(basepath, sessions{idx}, 'session_particulars.txt'));

if ~isfield(particulars, 'bad_channels') || ...
        ~isequal(particulars.bad_channels, badChans)
    particulars.bad_channels = badChans;
    fid = fopen(fullfile(basepath, sessions{idx}, 'session_particulars.txt'), 'a');
    fprintf(fid, 'bad_channels= %s\n', strjoin(string(badChans), ' '));
    fclose(fid);
end



%% 3. Extract chunks for each condition
% condition 1: Static Grating (on and off samples)
% condition 2: Phase Reversing Grating @ 40 Hz (sample at which frame one turns on and off—starts with white)
% condition 3: Full Screen Flickering @ 40 Hz (sample full screen flickers on then off)
% condition 4: Phase Reversing Grating @ 40 Hz (sample at which frame two turns on and off —starts with black)
% condition 5: Running speed (sample at which pulse triggered; 600 pulses per rotation)

condition = 1;

nCh = size(session.data, 1);
duration = size(session.data, 2) / fs; 

events = session.events.analogEvents{1, condition}; % convert to seconds
events = reshape(events, 2, []); % each row is [start, end] of event
triggers = events(1, :)./30e3; % take start times as triggers IN SECONDS

tWin = [5, 10]; % time window around event in seconds [pre, post]
movingWin = [0.2 0.2*0.2]; % moving window for chunking [window size, step size] in seconds

valid = triggers + tWin(2) <= duration; % ensure triggers fit within data duration
triggers = triggers(valid); 

% parameters for PSD computation
params.tapers = [2 3]; % time-bandwidth product 3, 5 tapers
params.Fs = fs;
params.fpass = [10 200]; % frequency range
params.pad = 1 ; % padding

% preallocate cell containers (good hygiene)
S      = cell(nCh,1);
S_bsl  = cell(nCh,1);
S_trig = cell(nCh,1);
S_diff = cell(nCh,1);
S_mu   = cell(nCh,1);
S_sem  = cell(nCh,1);

for ch = 1:nCh
    
    [S{ch}, t, f] = mtspecgramtrigc(session.data(ch, :), triggers, tWin, movingWin, params);

    % S is in dimensions: times x frequencies x trials

    baseline_mask = (t>=4.5 & t<5); % baseline period before event
    triggered_mask = (t>5); % period after event

    S_bsl{ch} = sq(mean(S{ch}(baseline_mask, :, :), 1)); % average over baseline time
    S_trig{ch} = sq(mean(S{ch}(triggered_mask, :, :), 1)); % average over triggered time
    S_diff{ch} = S_trig{ch} - S_bsl{ch}; % difference between triggered and baseline
    S_mu{ch} = sq(mean(S_diff{ch}, 2)); % mean over trials
    S_sem{ch} = sq(std(S_diff{ch}, 0, 2)) ./ sqrt(size(S{ch}, 3)); % SEM over trials

end

nF = size(S_mu{1}, 1);

S_mu_all  = zeros(nF, nCh);
S_sem_all = zeros(nF, nCh);

for ch = 1:nCh
    S_mu_all(:, ch)  = S_mu{ch};
    S_sem_all(:, ch) = S_sem{ch};
end

session.responses.spectra{condition} = struct(...
    'S_mu', S_mu_all, ...
    'S_sem', S_sem_all, ...
    'f', f ...
);



%% 4. Superimpose spectrum over channels for a single condition

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