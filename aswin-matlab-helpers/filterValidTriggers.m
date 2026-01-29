function triggers = filterValidTriggers(triggers, tWin, duration)
    valid = triggers + tWin(2) <= duration;
    triggers = triggers(valid);
end