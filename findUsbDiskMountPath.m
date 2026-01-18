function mountPath = findUsbDiskMountPath(targetLabel)
%FINDUSBDISKMOUNTPATH Return mount path for a disk with given volume label.
% Works on Windows, Linux, macOS. Returns '' if not found (or not mounted).

if nargin < 1 || strlength(string(targetLabel)) == 0
    error('Provide a target label, e.g. "EXELU_SSD1".');
end
targetLabel = char(string(targetLabel));

mountPath = '';

if ispc
    % Use PowerShell for robust label->drive-letter mapping
    % Outputs CSV: DriveLetter,FileSystemLabel
    ps = [
        'powershell -NoProfile -Command "'
        '$v=Get-Volume | Select-Object DriveLetter,FileSystemLabel; '
        '$v | ConvertTo-Csv -NoTypeInformation"'
    ];
    [status, out] = system(ps);
    if status ~= 0 || isempty(strtrim(out)), return; end

    lines = splitlines(string(strtrim(out)));
    % First line is header; subsequent lines are CSV rows: "C","OS"
    for i = 2:numel(lines)
        line = strtrim(lines(i));
        if line == "", continue; end
        parts = parseCsvLine(line);
        if numel(parts) < 2, continue; end

        drive = strip(parts{1}, '"');
        label = strip(parts{2}, '"');

        if strcmpi(char(label), targetLabel) && strlength(drive) > 0
            mountPath = char(drive + ":\");
            return;
        end
    end

elseif ismac
    % Most macOS external volumes mount under /Volumes/<Label>
    candidate = fullfile('/Volumes', targetLabel);
    if isfolder(candidate)
        mountPath = candidate;
        return;
    end

    % Fallback: query diskutil for mounted volumes and labels
    cmd = 'diskutil list -plist';
    [status, ~] = system(cmd);
    if status ~= 0
        return; % if diskutil plist not available, keep simple
    end

    % Practical fallback: scan /Volumes for case-insensitive match
    d = dir('/Volumes');
    for k = 1:numel(d)
        if d(k).isdir && ~startsWith(d(k).name, '.')
            if strcmpi(d(k).name, targetLabel)
                mountPath = fullfile('/Volumes', d(k).name);
                return;
            end
        end
    end

elseif isunix
    % Linux: get LABEL + MOUNTPOINT for mounted block devices
    % -n: no header, -r: raw, -p: full device paths (optional)
    cmd = 'lsblk -o LABEL,MOUNTPOINT -nr';
    [status, out] = system(cmd);
    if status ~= 0 || isempty(strtrim(out)), return; end

    lines = splitlines(string(strtrim(out)));
    for i = 1:numel(lines)
        line = strtrim(lines(i));
        if line == "", continue; end

        % Split into: LABEL MOUNTPOINT (mountpoint may be missing)
        parts = regexp(char(line), '\s+', 'split');
        if isempty(parts), continue; end

        label = parts{1};
        mnt = '';
        if numel(parts) >= 2
            % mountpoint itself can contain spaces rarely; lsblk prints it as one field,
            % but regexp split would separate it. So re-join the rest:
            mnt = strtrim(strjoin(parts(2:end), ' '));
        end

        if strcmpi(label, targetLabel) && ~isempty(mnt)
            mountPath = mnt;
            return;
        end
    end
else
    error('Unsupported OS.');
end
end

function parts = parseCsvLine(line)
% Minimal CSV parser for PowerShell ConvertTo-Csv output lines like: "E","EXELU_SSD1"
% Handles quoted fields and commas.
line = char(line);
parts = {};
field = '';
inQuotes = false;
for i = 1:numel(line)
    c = line(i);
    if c == '"'
        inQuotes = ~inQuotes;
        field(end+1) = c; %#ok<AGROW>
    elseif c == ',' && ~inQuotes
        parts{end+1} = strtrim(field); %#ok<AGROW>
        field = '';
    else
        field(end+1) = c; %#ok<AGROW>
    end
end
parts{end+1} = strtrim(field);
end
