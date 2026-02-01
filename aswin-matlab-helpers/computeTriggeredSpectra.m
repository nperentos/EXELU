function spec = computeTriggeredSpectra(data, fs, triggers, tWin, movingWin, params)
% computeTriggeredSpectra
%
% Pure function: given data + triggers, compute baseline-corrected spectra.
%
% Inputs:
%   data       : channels × samples
%   fs         : sampling rate (Hz)
%   triggers   : vector of trigger times (seconds), one per trial
%   tWin       : [pre post] window in seconds (e.g. [5 10])
%   movingWin  : [window step] in seconds
%   params     : Chronux params struct (Fs, tapers, fpass, pad)
%
% Output:
%   spec : struct with fields
%       .S_mu   (freq × channel)
%       .S_sem  (freq × channel)
%       .f      frequency axis (Hz)

nCh = size(data, 1);

% preallocate per-channel containers
S_diff = cell(nCh,1);
S_mu  = cell(nCh,1);
S_sem = cell(nCh,1);

for ch = 1:nCh

    [S, t, f] = mtspecgramtrigc(data(ch,:), triggers, tWin, movingWin, params);
    % S: time × freq × trials

    baseline_mask  = (t >= 4.5 & t < 5);
    triggered_mask = (t > 5);

    S_bsl  = sq(mean(S(baseline_mask,:,:), 1));
    S_trig = sq(mean(S(triggered_mask,:,:), 1));
    S_diff{ch} = S_trig - S_bsl;

    S_mu{ch}  = sq(mean(S_diff{ch}, 2));
    S_sem{ch} = sq(std(S_diff{ch}, 0, 2)) ./ sqrt(size(S,3));
end

% pack into arrays (freq × channel)
nF = size(S_mu{1},1);
nTr = size(S_diff{1},2);

S_diff_all = zeros(nF, nCh, nTr);
S_mu_all  = zeros(nF, nCh);
S_sem_all = zeros(nF, nCh);

for ch = 1:nCh
    S_diff_all(:,ch,:) = S_diff{ch};
    S_mu_all(:,ch)  = S_mu{ch};
    S_sem_all(:,ch) = S_sem{ch};
end

spec = struct( ...
    'S_diff', S_diff_all, ...
    'S_mu',  S_mu_all, ...
    'S_sem', S_sem_all, ...
    'f',     f ...
);
end