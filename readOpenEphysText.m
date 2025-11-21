function messages = readOpenEphysText(filename)
% readOpenEphysText Loads Open Ephys message center text.npy files into MATLAB
%
%   messages = readOpenEphysText(filename)
%
% REQUIREMENTS:
%   - MATLAB with Python configured (`pyenv` must point to a Python installation with NumPy).
%   - Works with Open Ephys `text.npy` files (pickled object arrays).
%
% OUTPUT:
%   - messages: cell array of MATLAB character vectors

    if ~isfile(filename)
        error('File not found: %s', filename);
    end

    % Ensure Python is available and has NumPy
    try
        py.importlib.import_module('numpy');
    catch
        error('Python with NumPy is required. Install NumPy and ensure MATLAB is using the correct Python (see "pyenv").');
    end

    % Load using numpy with allow_pickle=True
    npData = py.numpy.load(filename, pyargs('allow_pickle', true));
    npList = npData.tolist();        % Convert NumPy array to Python list
    messages = cellfun(@char, cell(npList), 'UniformOutput', false); % Convert to MATLAB cell array of strings
end