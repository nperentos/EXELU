function sessionStruct = loadLFPSession(basepath, session, ch)
% loadSession 
%
% sessionStruct.name
% sessionStruct.chanmap
% sessionStruct.events
% sessionStruct.data
% sessionStruct.settings
% sessionStruct.tScale
%
% basepath : root directory containing session folders
% session  : session folder name (char or string)
% ch       : optional subset of channels

    sessionPath = fullfile(basepath, session);
    sessionStruct.name = session;

    % metadata local to session folder
    sessionStruct.chanmap = load(fullfile(sessionPath, session + ".chanmap.mat"));
    sessionStruct.events  = load(fullfile(sessionPath, "analog_events.mat"));

    % LFP via original getLFP function
    if nargin >= 3 && ~isempty(ch)
        [sessionStruct.data, sessionStruct.settings, sessionStruct.tScale] = ...
            getLFP(session, basepath, ch);
    else
        [sessionStruct.data, sessionStruct.settings, sessionStruct.tScale] = ...
            getLFP(session, basepath);
    end
end
