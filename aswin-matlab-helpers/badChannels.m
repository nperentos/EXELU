
function badChans = badChannels(session, fs, segmentDuration)
% badChannels - Identify bad channels in an LFP session based on signal quality.
% Syntax: badChans = badChannels(session)
% Inputs:
%    session - An LFP session structure containing channel data.
%    fs      - Sampling frequency of the LFP data.
% Outputs:          
%    badChans - A vector of channel indices identified as bad.

if nargin < 2
    fs = session.settings.parameters.acquisitionSystem.samplingRate;
    segmentDuration = 100; % ms
end

% select samples in first segmentDuration and normalize by median
segment = session.data(:, 1:segmentDuration*fs)'; % transpose to have time x channels
segment = segment - median(segment, 1); % remove DC offset

% compute correlation between adjacent channels of linear probe
C = size(segment, 2);
rho = zeros(C-1, 1);
for i = 1:C-1
    r = corr(segment(:, i), segment(:, i+1), 'type', 'Spearman');
    rho(i) = r;
end

% compute score for each channel based on correlation with neighbors
score = nan(C, 1);
score(2:C-1) = (rho(1:end-1) + rho(2:end))/2;
z = (score - median(score, 'omitnan')) / 1.4826*mad(score, 1);
badChans = find(z < -3); % threshold at z < -3

% compute segment PSD for each channel
movingwindow = [1.5 1.5*0.2]; % 1.5 ms window, 20% overlap
params.tapers = [2 3] ; % time-bandwidth product 2, 3 tapers
params.Fs = fs;
params.fpass = [0 200]; % frequency range 0-200 Hz
params.pad = 0 ; % no padding

[S_seg, t_seg, f_seg] = mtspecgramc(segment, movingwindow, params);
S_seg_mean = sq(mean(S_seg, 1)); % average over time segments

lineHz = 50; % line noise frequencies
bw = 1 ; % bandwidth around line noise to consider
lineBand = f_seg >= lineHz - bw & f_seg <= lineHz + bw;
linePower = sum(S_seg_mean(lineBand, :), 1);
lfpBand = f_seg > 1 & f_seg < 200 & ~lineBand;
lfpPower = sum(S_seg_mean(lfpBand, :), 1);
lineRatio = linePower ./ lfpPower;

% combine metrics
score = sqrt(max(z, 0)) .* max(zscore(lineRatio)', 0);
badChans = find(score > 3); % threshold at score > 3
disp(['Bad channels identified: ' mat2str(badChans)]);
end

