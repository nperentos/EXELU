function vars = loadParticulars(filename)
% loadVars  Load name=value pairs from a text file into a struct.
%
%   vars = loadVars('vars.txt')
%
% The file should contain lines like:
%   CA1_chan = 55
%   theta = 67
%   day = tues
%
% Numeric values will be converted automatically. Text is kept as string.

    % Read table, but don't attempt automatic type detection
    T = readtable(filename, ...
        'Delimiter', '=', ...
        'ReadVariableNames', false, ...
        'TextType', 'string', ...   % keep as string when possible
        'Format', '%s%s');          % force both columns to be text

    % Convert the columns to strings no matter what
    names  = string(T.Var1);
    values = string(T.Var2);

    % Trim whitespace
    names  = strtrim(names);
    values = strtrim(values);

    vars = struct();

    for i = 1:numel(names)
        v = values(i);

        % Attempt numeric conversion
        num = str2double(v);

        if ~isnan(num)
            vars.(names(i)) = num;
        else
            vars.(names(i)) = v;
        end
    end
end