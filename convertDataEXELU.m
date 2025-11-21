function [] = convertDataEXELU(fileBase,basepath)
%% PRE
if nargin < 2
    basepath = 'C:\Users\MDRTC5707A\Documents\Open Ephys\';
end
pth = fullfile(basepath,fileBase);

% check if folder exists
if ~exist(pth, 'dir') == 7
    error('recording folder not found. Please check');
end

%% RAW DATA BACKUP
% first, backup the original onto an external SSD EXELU_SSD2 for example
% Check for external drive named 'EXELU'
[~, result] = system('wmic logicaldisk get VolumeName,DeviceID');

% Look for the line containing 'EXELU'
driveLetter = '';
lines = strsplit(result, newline);
for i = 1:length(lines)
    if contains(lines{i}, 'EXELU_SSD2')
        parts = strsplit(strtrim(lines{i}));
        driveLetter = parts{1}; % e.g. 'E:'
        break;
    end
end

% Use the result
if ~isempty(driveLetter)
    externalPath = [driveLetter '\']; % e.g. 'E:\'
    fprintf('External EXELU_SSD2 found at %s\n', externalPath);
    sourceFolder = pth;
    destinationFolder = fullfile(externalPath, fileBase);
    if ~exist(destinationFolder, 'dir')
        copyfile(sourceFolder, destinationFolder);
        disp('Folder copied successfully.');
    else
        disp('Folder already exists on external SSD.');
    end
else
    disp('External SSD EXELU not found.');
end

%% Create LFP and dat files
% load metadata
xml_pth = fullfile(pth,'Record Node 102','experiment1','recording1','structure.oebin');  
%  assuming the folder structure wont change!
[sample_rate, num_channels, channels_description] = metadata_load(xml_pth);

% load timestamps
timestamps_pth = fullfile(pth,'Record Node 102\experiment1\recording1\continuous\Acquisition_Board-100.acquisition_board','timestamps.npy');
timestamps = readNPY(timestamps_pth);

% check if any samples are missings
% tmp = unique(diff(timestamps));
% exc = find(tmp>1/sample_rate);
% if (max(tmp) - min(tmp)) >= 1/sample_rate
%     error('discontinuous samples detected')
% end

% make an xml for neuroscope
cd(pth)
%neuroscope_xml_creator(fileBase, sample_rate, num_channels,1000); % the standard one
neuroscope_xml_creator_EXELU(fileBase, sample_rate, num_channels,1000,channels_description);% the EXELU specific one

% % hard link to the bin file at the filebase 
% need to ensure that an actual copy is generated when the raw data are
% saved elsewhere to the processing data
link = [fileBase,'.dat'];
target = fullfile(pth,'Record Node 102\experiment1\recording1\continuous\Acquisition_Board-100.acquisition_board','continuous.dat');
cmd = sprintf('fsutil hardlink create "%s" "%s"',link,target);
[status,result] = system(cmd);
% for now we avoid making a new dat file with bit2volts taken into account
% as there is a possibility ot accuracy loss (unless we use float32 which will increase file sizes substantially
% . We just need to remember in
% later stages to do the bit to volt conversions

% make .lfp
dat2lfp([fileBase,'.dat']);

% extract events by loading ADCs only
ev_chs = find([channels_description.type]==2);
fid = fopen([fileBase,'.dat']);
fseek(fid,0,'eof');
len = ftell(fid); % filesize in bytes
len = len/(num_channels*2);
m = memmapfile(link,'Format',{'int16',[num_channels len],'m'},'writable',false);

d = m.Data.m;
ADCs = d(ev_chs,:);
for i = 1:size(ADCs,1)
    display('')
    test = ADCs(i,:);
    tmp1 = find(test(1:end-1)>max(test)*.75); % if this occasionaly doesnt work we have to switch to prctile
    tmp2 = find(test(2:end)>max(test)*.75);
    [both, fall, rise] = setxor(tmp1,tmp2);
    analogEvents{i} = rise; % add names of event channels as per arduino connections and correspondece to the py code conditions
end
analog_ch_names = {'static','flashing','reversing','encoder'};

% get any user text messages 
pth_txt_msg = fullfile(basepath,fileBase,'Record Node 102\experiment1\recording1\events\MessageCenter');
ts_user_msgs = readNPY(fullfile(pth_txt_msg,'timestamps.npy'));
user_msgs = readOpenEphysText(fullfile(pth_txt_msg,'text.npy')); % doesnt work

save('analog_events.mat','analogEvents','analog_ch_names','user_msgs','ts_user_msgs');


%% SPIKESORTING

% make a map file specific to H3(64) + EEG(6) + ADC(5)
ephys_chs = find([channels_description.type]==0);
Nchannels = num_channels;%length(ephys_chs);
chanMap = 1:Nchannels;
connected = false(Nchannels,1);
connected(ephys_chs) = true; % if we always have 6 EEGs
xcoords = zeros(Nchannels,1);
ycoords = (0:20:(20*(Nchannels-1)))';
kcoords = zeros(Nchannels,1);
fs = sample_rate;
save([fileBase,'.chanmap.mat'],'chanMap','connected','xcoords','ycoords','kcoords','fs');
% figure; % scatter(xcoords,ycoords,[],connected); % set(gca, 'YDir','reverse');

%vconfigure python for kilosort
py_exe = 'C:\Users\MDRTC5707A\anaconda3\envs\kilosort\python.exe'; % Python inside the conda env
script = 'C:\Users\MDRTC5707A\Documents\MATLAB\KSpy\run_KS_wrapper_NP.py';
SAVE_PATH = fullfile(pth,[fileBase,'.dat']);
N_CHAN_BIN = num_channels;
MAP = fullfile(pth,[fileBase,'.chanmap.mat']);

% === BUILD AND RUN COMMAND ===
logFile = fullfile(pth,'kilosort_log.txt');

cmd = sprintf('"%s" "%s" "%s" %d "%s"', py_exe, script, SAVE_PATH, N_CHAN_BIN, MAP);
%status = system(cmd);
status = system([cmd ' > "' logFile '" 2>&1']);

if status == 0
    disp('Kilosort completed successfully.');
else
    warning('ks algo script exited with error code %d.', status);
end


script2 = 'C:\Users\MDRTC5707A\Documents\MATLAB\KSpy\generate_KS_report.py';
command = sprintf('"%s" "%s" "%s"', py_exe, script2, SAVE_PATH);
status = system(command);
if status == 0
    disp('Kilosort report generated successfully.');
else
    warning('ks report script exited with error code %d.', status);
end

%% POST
% once everything is finished,copy the processed data to EXELU_SSD1
% (without the original OE files)
% Check for external drive named 'EXELU'
[~, result] = system('wmic logicaldisk get VolumeName,DeviceID');

% Look for the line containing 'EXELU'
driveLetter = '';
lines = strsplit(result, newline);
for i = 1:length(lines)
    if contains(lines{i}, 'EXELU_SSD1')
        parts = strsplit(strtrim(lines{i}));
        driveLetter = parts{1}; % e.g. 'E:'
        break;
    end
end

% Use the result
if ~isempty(driveLetter)
    externalPath = [driveLetter '\']; % e.g. 'E:\'
    fprintf('External EXELU_SSD1 found at %s\n', externalPath);
    sourceFolder = pth;
    destinationFolder = fullfile(externalPath, fileBase);
    excludeFolder  = 'Record Node 102';
    if ~exist(destinationFolder, 'dir')
        mkdir(destinationFolder);
        items = dir(sourceFolder);
        disp('copying please wait...');
        for i = 1:length(items)
            name = items(i).name;
        
            % Skip '.' and '..'
            if strcmp(name, '.') || strcmp(name, '..')
                continue;
            end
        
            % Full paths
            srcPath = fullfile(sourceFolder, name);
            destPath = fullfile(destinationFolder, name);
        
            % Skip the excluded folder
            if items(i).isdir && strcmp(name, excludeFolder)
                fprintf('Skipping folder: %s\n', srcPath);
                continue;
            end
        
            % Copy files or folders
            if items(i).isdir
                copyfile(srcPath, destPath); % Copies entire folder recursively
            else
                copyfile(srcPath, destPath);
            end

        end
        disp('Folder copied successfully.');
    else
        disp('Folder already exists on external SSD.');
    end
else
    disp('External SSD_EXELU1 not found.');
end

end