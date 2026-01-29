
function badChans = badChannels(session, fs, segmentDuration, W)
% badChannels - Identify bad channels in an LFP session based on signal quality.
% Syntax: badChans = badChannels(session)
% Inputs:
%    session - An LFP session structure containing channel data.
%    fs      - Sampling frequency of the LFP data.
%    segmentDuration - Duration (in seconds) of the segment to analyze for bad channels.
%    W       - Window size for local z-score computation (default: 5).
% Outputs:          
%    badChans - A vector of channel indices identified as bad.

if nargin < 2 || isempty(fs)
    fs = session.settings.parameters.fieldPotentials.lfpSamplingRate;
end
if nargin < 3 || isempty(segmentDuration)
    segmentDuration = 0.1; 
end
if nargin < 4 || isempty(W)
    W = 5;
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

% compute z-score for each channel based on correlation with neighbors
score = nan(C, 1);
score(2:C-1) = (rho(1:end-1) + rho(2:end))/2;

% local z-score (within groups of ±W channels)
z_local_corr = nan(size(score)); 
for ch = 1:numel(score)
    lo = max(1, ch-W); % lower bound
    hi = min(numel(score), ch+W); % upper bound
    local = score(lo:hi); 
    z_local_corr(ch) = (score(ch) - median(local,'omitnan')) / (1.4826 * mad(local,1));
end

bad_corr_idx = z_local_corr < -3; % threshold at z < -3

% compute segment PSD for each channel
movingwindow = [0.2 0.04]; % 200 ms window, 20% step size
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
lineRatio = (linePower ./ lfpPower)';

% compute z-score for line noise ratio
z_line = (lineRatio - median(lineRatio)) / (1.4826*mad(lineRatio, 1));
bad_line_idx = z_line > 3; % threshold at z > 3

% combine metrics
% score = sqrt(max(z, 0)) .* max(zscore(lineRatio)', 0);

badChans = find(bad_corr_idx | bad_line_idx);

fprintf('Bad channels identified: %s\n', mat2str(badChans));
end

