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
fs_data  = str2double(session.settings.parameters.fieldPotentials.lfpSamplingRate.Text);

%% 2. Load particulars files for each session and/or compute bad channels if necessary

% load existing particulars file and update bad channels if necessary
particulars = loadParticulars(fullfile(basepath, sessions{idx}, 'session_particulars.txt'));

% compute bad channels if not already present
if ~isfield(particulars,'bad_channels') || isempty(particulars.bad_channels)
    segmentDuration = 0.1; % seconds
    badChans = badChannels(session, fs_data, segmentDuration);
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
duration = size(data, 2) / fs_data; 

% parameters for multitaper spectral estimation
tWin = [5 10]; % [pre post] seconds
movingWin = [0.2 0.2*0.2]; % 200 ms window, 20% step size
params = struct();
params.tapers = [2 3];
params.Fs     = fs_data; 
params.fpass  = [10 200];
params.pad    = 1;

%% 3.1 Extract chunks for condition 1
% Static Grating (on and off samples)
 
fs_events = 30e3; % event sampling rate (Hz)

condition = 1;

% extract triggers
events = session.events.analogEvents{1, condition};
events = reshape(events, 2, [])'; % [start end] per trial
triggers = events(:,1) ./ fs_events;  % trial onsets (seconds)

% clean triggers to ensure enough data for pre and post windows
validTriggers = filterValidTriggers(triggers, tWin, duration);

% compute triggered spectra
spec = computeTriggeredSpectra( ...
    data, fs_data,  validTriggers, tWin, movingWin, params ...
);
session.responses.spectra{condition} = spec;

%% 3.2 Extract chunks for condition 2
% Full Screen Flickering @ 40 Hz (sample full screen flickers on then off)

condition = 2;
pulses = session.events.analogEvents{1, condition} ./ fs_events; 
gapThresh = 5; % seconds
triggers = extractTrialTriggers(pulses, gapThresh);
validTriggers = filterValidTriggers(triggers, tWin, duration);

spec = computeTriggeredSpectra( ...
    data, fs_data, validTriggers, tWin, movingWin, params ...
);
session.responses.spectra{condition} = spec;

%% 3.3 Extract chunks for condition 3
% Phase Reversing Grating @ 40 Hz (sample at which frame one turns on and off—starts with white)

condition = 3;

% trigger extraction (implicit trials via pulse gaps) 
pulses = session.events.analogEvents{1, condition} ./ fs_events;   % seconds
gapThresh = 5;  % seconds
triggers = extractTrialTriggers(pulses, gapThresh);
validTriggers = filterValidTriggers(triggers, tWin, duration);

% spectra computation (generalised, unchanged) 
spec = computeTriggeredSpectra( ...
    data, fs_data, validTriggers, tWin, movingWin, params ...
);

session.responses.spectra{condition} = spec;

%% 4.1 Plot depth x frequency heatmap per condition

nConditions = numel(session.responses.spectra);

f = session.responses.spectra{1}.f;
fMask = f >= 5 & f <= 150;
xLimits = [f(find(fMask,1,'first')) f(find(fMask,1,'last'))];

colorlim = [-5, 5]; % dB color limits

test = session.responses.spectra{1}.S_mu;
summary(test)

figure('Color', 'w', 'Position', [100 100 1200 1200]);
tiledlayout(nConditions,1, 'Padding', 'compact', 'TileSpacing', 'compact');

for c = 1:nConditions
    nexttile;

    S = session.responses.spectra{c}.S_mu; % freq × channel
    f = session.responses.spectra{c}.f;

    img = S(fMask, :)'; % apply frequency mask and transpose
    imagesc(f(fMask), 1:size(S,2), 10*log10(abs(img).^2));

    axis xy;
    % set(gca);
    xlim(xLimits);
    clim(colorlim);
    colormap(gca, parula);

    xline(40, '--w', 'LineWidth', 1.5); % 60 Hz line

    xlabel('Frequency (Hz)');
    ylabel('Channel (deep → superficial)');
    title(sprintf('Condition %d: Mean Triggered Spectra', c));
    
    set(gca, 'FontSize', 12, 'LineWidth', 1);

end

cb = colorbar;
cb.Label.String = 'Power (dB)';
cb.Layout.Tile = 'east';

%% 4.2 Channel 22 across conditions

ch = 22;

f = session.responses.spectra{1}.f;
fMask = f >= 5 & f <= 150;

figure('Color','w'); hold on;

for c = 1:numel(session.responses.spectra)
    spec = session.responses.spectra{c};

    mu  =  10*log10(abs(spec.S_mu(fMask, ch)).^2);
    sem =  10*log10(abs(spec.S_sem(fMask, ch)).^2);

    plot(f(fMask), mu, 'LineWidth', 3);

    % optional SEM
    fill([f(fMask) fliplr(f(fMask))], ...
         [mu-sem; flipud(mu+sem)]', ...
         'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
end

set(gca, 'XScale', 'log');
xline(40,'--k','LineWidth',1.5);

xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title(sprintf('Single-channel spectra (ch %d)', ch));

legend({'Cond 1','SEM','Cond 2','SEM','Cond 3','SEM'}, 'Location','best');
set(gca,'FontSize',12,'LineWidth',1);

%% 4.3 Gamma power across conditions
ch = 22;
conds = numel(session.responses.spectra);

f = session.responses.spectra{1}.f;
gammaMask = f >= 30 & f <= 50;
[~, f40] = min(abs(f - 40));

figure('Color','w'); hold on;

for c = 1:conds
    spec = session.responses.spectra{c};

    % extract all trial-wise data
    pAll = squeeze(spec.S_diff(:, ch, :));   % [freq × trial]
    pAll_db = 10*log10(abs(pAll));

    % --- non-gamma frequencies (light grey) ---
    nongamma = find(~gammaMask);
    for k = nongamma'
        x = c + 0.15*randn(1, size(pAll_db,2));
        scatter(x, pAll_db(k,:), 8, ...
            'MarkerEdgeColor', [0.8 0.8 0.8], ...
            'MarkerFaceAlpha', 0.1);
    end

    % --- gamma frequencies except 40 Hz (dark grey) ---
    gammaOther = find(gammaMask);
    gammaOther(gammaOther == f40) = [];
    for k = gammaOther'
        x = c + 0.15*randn(1, size(pAll_db,2));
        scatter(x, pAll_db(k,:), 10, ...
            'MarkerEdgeColor', [1 0 1], ...
            'MarkerFaceAlpha', 0.5);
    end

    % --- 40 Hz (red) ---
    x = c + 0.15*randn(1, size(pAll_db,2));
    scatter(x, pAll_db(f40,:), 18, ...
        'MarkerEdgeColor', [1 0 0], ...
        'LineWidth', 1.2);

    % --- band-averaged mean ± SEM (gamma band) ---
    pGamma = squeeze(mean(spec.S_diff(gammaMask, ch, :),1));
    pGamma_db = 10*log10(abs(pGamma));

    mu  = mean(pGamma_db);
    sem = std(pGamma_db) / sqrt(numel(pGamma_db));

    errorbar(c, mu, sem, 'ko', ...
        'MarkerFaceColor','k', ...
        'LineWidth',2);
end

xlim([0.5 conds+0.5]);
xticks(1:conds);
xticklabels({'Cond 1','Cond 2','Cond 3'});
ylabel('Power change (dB)');
ylim([0 40]);
title(sprintf('Spectral power by condition (ch %d)', ch));
set(gca,'FontSize',12,'LineWidth',1);

%% NOTES

% condition 3: Full Screen Flickering @ 40 Hz (sample full screen flickers on then off)
% condition 4: Phase Reversing Grating @ 40 Hz (sample at which frame two turns on and off —starts with black)
% condition 5: Running speed (sample at which pulse triggered; 600 pulses per rotation)
