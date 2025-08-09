function [out] = load_OE_data_EXELU(session_path)
%% loads data into matlab workspace and creates an xml for neuroscope
% to do: generate a .lfp file 
%        see how to incoroporate into cellexplorer to take advantage of that pipeline

% check if folder is present
if ~exist(session_path, 'dir') == 7
    error('recording folder not found. Please check');
end

try
    % load metadata
    pth = fullfile(session_path,'structure.oebin');    
    [sample_rate, num_channels, channels_description] = metadata_load(pth);
    
    % load timestamps
    pth = fullfile(session_path,'continuous','Acquisition_Board-100.acquisition_board','timestamps.npy');
    timestamps = readNPY(pth);

    % load actual data
    pth = fullfile(session_path,'continuous','Acquisition_Board-100.acquisition_board','continuous.dat');
    fle = fopen(pth, 'r');
    raw_data = fread(fle,[num_channels,Inf], 'uint16');
    fclose(fle);
    raw_data = (raw_data'.*[channels_description.bit_volts])';
catch
    error('there was an error. Either there are files missing or files are corrupt?');
end

out.t = timestamps;
out.d = raw_data;
out.fs = sample_rate;
out.nCh = num_channels;
out.ch_description = channels_description;

neuroscope_xml_creator(fullfile(session_path,'continuous','Acquisition_Board-100.acquisition_board','settings'), sample_rate, num_channels,1000);

% % this is how we checked for the pulses to verify that they are correct
% test = raw_data(71,:);
% tmp1 = find(test(1:end-1)>max(test)*.75);
% tmp2 = find(test(2:end)>max(test)*.75);
% [both, fall, rise] = setxor(tmp1,tmp2);
% a = diff(both');
% a(a>550)=[];
% figure;
% hist(a,40);
