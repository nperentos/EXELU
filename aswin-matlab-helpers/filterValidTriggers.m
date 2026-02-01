function validTriggers = filterValidTriggers(triggers, tWin, duration)
    valid = triggers + tWin(2) <= duration;
    validTriggers = triggers(valid);
end