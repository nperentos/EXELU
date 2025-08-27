function [outputArg1,outputArg2] = convertDataEXELU(fileBase,basepath)
% Create LFP and dat files
%   Detailed explanation goes here
if nargin < 2
    basepath = 'C:\Users\MDRTC5707A\Documents\Open Ephys\';
end
pth = fullfile(basepath,fileBase);

% check if folder exists
if ~exist(pth, 'dir') == 7
    error('recording folder not found. Please check');
end

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
% as there is a possibility ot accuracy loss. We just need to remember in
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

save('analog_events.mat','analogEvents','analog_ch_names');


% spikesorting
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
% figure;
% scatter(xcoords,ycoords,[],connected);
% set(gca, 'YDir','reverse');

% % we postpone trying to run KS from matlab for now 
% % invoke the KS py env and spike sort with default vars
% %pyenv('Version','C:\Users\MDRTC5707A\anaconda3\envs\kilosort\python.exe'); % run only once to set matlab default python
% py.print("hello from KS env!");
% py.runpy.run_path('C:\Users\MDRTC5707A\Documents\MATLAB\KSpy\KSpy.py');


outputArg1 = 0;
outputArg2 = 0;

end