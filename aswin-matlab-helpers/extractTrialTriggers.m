function triggers = extractTrialTriggers(pulses, gapThresh)
% pulses    : vector of event times (seconds)
% gapThresh : minimum gap (seconds) that separates trials
% triggers  : one trigger per trial (seconds)

    pulses = pulses(:);                 % enforce column
    gapIdx = find(diff(pulses) > gapThresh);
    startIdx = [1; gapIdx + 1];
    triggers = pulses(startIdx);
end