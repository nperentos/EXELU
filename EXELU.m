% 27/09/2025
% first script for basic data processing of EXELU data post dat and lfp creation
% KS has already ran at this stage by another script that resides on the recording PC
% the script is called convertDataEXELU.m

addpath(genpath('C:\Users\npere\MATLAB'));

dataPath = 'I:\';
fileBase = '2025-09-20_16-51-22';
session_path = fullfile(dataPath,fileBase);
load(fullfile(dataPath,fileBase,[fileBase,'.chanmap.mat']));
load(fullfile(dataPath,fileBase,'analog_events.mat'));
num_channels = length(chanMap);
[data,settings,tScale] = getLFP(fileBase, dataPath);

SR = str2num(settings.parameters.fieldPotentials.lfpSamplingRate.Text);

particulars = loadParticulars(fullfile(session_path,'session_particulars.txt'));

% keep the probe only the following parts
x = data(1:128,:); 
clear data; % just to keep some space

% theta
    movingwin = [1.5 1.5*0.2];
    params.tapers = [3 5];    % Or [2 3] for slightly smoother spectra
    params.Fs     = SR;     % Sampling rate (Hz)
    params.fpass  = [0 60];   % Focus on theta band
    params.pad    = 0;        % No additional frequency padding
    % params.win    = [1 0.5];  % 1 s window, 0.5 s step (50% overlap)
    theta_ch = particulars.CA1_chan-7;
    tic;
    [S,t,f]=mtspecgramc(data',movingwin,params); %(theta_ch,:) %28 seconds!
    toc
    % figure;
    % for i = 1:num_channels
    %     imagesc(t,f,S(:,:,i)');
    %     axis xy; title(['ch ',num2str(i)]);
    %     pause(0.2); cla;
    % end


% gamma 

    movingwin = [0.2 0.2*0.2]; % 0.2 s x 20% overlap step
    params.tapers = [2 3];    % Or [2 3] for slightly smoother spectra
    params.Fs     = SR;     % Sampling rate (Hz)
    params.fpass  = [10 200];   % Focus on gamma band
    params.pad    = 1;        % No additional frequency padding
    % params.win    = [1 0.5];  % 1 s window, 0.5 s step (50% overlap)
    theta_ch = particulars.CA1_chan-7;
    
    % whitening is desirable
    whitenorder = 2;
    whitenwindow = 20 * SR; % 20sec window

    if size(x,1)<size(x,2); x = x'; end;
    if whitenwindow>size(x,1); whitenwindow = size(x,1); end; % In case signal is less than 20 sec

    for c=1:size(x,2);% cycle channels
        [A ,nv] = arburg(double(x(1:whitenwindow,c)),whitenorder);
        x(:,c) = filter(A,1,double(x(:,c)));
        options.NoiseVariance(c) = nv;
        options.WhitenModel(c,:) = A;
    end


    tic;
    [S,t,f]=mtspecgramc(x,movingwin,params); %(theta_ch,:) %28 seconds!
    toc
    % figure;
    % for i = 1:num_channels
    %     imagesc(t,f,S(:,:,i)');
    %     axis xy; title(['ch ',num2str(i)]);
    %     pause(0.2); cla;
    % end


% gather the relevant chunks for each condition. The events are in dat timescale 
% lets convert them to  

%static
    test = analogEvents{1,1}./30e3;
    test = reshape(test,2,[]);
    % t indices in the spectra matrix demarcating the start of the static exp.
    %idx = arrayfun(@(x) find(abs(t - x) == min(abs(t - x)), 1), test);
    idx = interp1(t, 1:numel(t), test, 'nearest');
    SegLen = median(diff(idx,1)); % should be about 10 s
    M = idx(1,:)' + (0:SegLen);
    M_static = vc(M');
    % figure;
    % for i = 1:num_channels
    %     imagesc(t(M_static),f,S(M_static,:,i)');
    %     axis xy; title(['ch ',num2str(i)]);
    %     pause(0.2); cla;
    % end





% flashing
    % find diffs in pulse times that exceed 5s
    test = analogEvents{1,2}./30e3;
    
    idx_end = [find(diff(test) > 5) ];
    idx_strt= [1 idx_end(1:end-1)+1];

    idx = interp1(t, 1:numel(t), test(idx_strt), 'nearest');
    M = idx(1,:)' + (0:SegLen);
    M_flashing = vc(M');
    




% reversing
    % find diffs in pulse times that exceed 5s
    test = analogEvents{1,3}./30e3;
    
    idx_end = [find(diff(test) > 5) ];
    idx_strt= [1 idx_end(1:end-1)+1];

    idx = interp1(t, 1:numel(t), test(idx_strt), 'nearest');
    M = idx(1,:)' + (0:SegLen);
    M_reversing = vc(M');
    


% compare same channel all spectra v static
    
    figure; 
    subplot(311); imagesc(t(M_static),f   ,S(M_static,:,55)'); axis xy; title('static');colorbar;
    subplot(312); imagesc(t(M_flashing),f ,S(M_flashing,:,55)'); axis xy; title('flashing');colorbar;
    subplot(313); imagesc(t(M_reversing),f ,S(M_reversing,:,55)'); axis xy; title('reversing');colorbar;
    mod_all_plots('clim([0 200])');




 



    %out = specmt_low_mem(data([theta_ch:theta_ch+2],:)',{'defaults','gamma','blocksize',2^10}); %data(1:2:end,:) 
out = specmt_low_mem(data([1:3],:)',{'defaults','gamma','blocksize',2^10}); %data(1:2:end,:) 
%%




test = analogEvents{1};
test =reshape(test,2,floor(length(test)/2));
durations = test(2,:)-test(1,:);











% % we need to run some corrections on the event detection since the script used had 
% % an incorrect event indexing method
% fbs = {'2025-09-18_09-59-40','2025-09-18_20-10-12','2025-09-20_10-49-09', ...
% '2025-09-20_13-33-23','2025-09-20_16-51-22'};
% dataPath = 'I:\';
% 
% for j = 1:length(fbs)
%     fileBase = fbs{j};
%     load(fullfile(dataPath,fileBase,'analog_events.mat'));
%     load(fullfile(dataPath,fileBase,[fileBase,'.chanmap.mat']));
%     num_channels = length(chanMap);
%     fle = fullfile(dataPath,fileBase,[fileBase,'.dat']);
%     fid = fopen(fle);
%     fseek(fid,0,'eof');
%     len = ftell(fid);
%     len = len/(num_channels*2);
%     m = memmapfile(fle,'Format',{'int16',[num_channels len],'m'},'writable',false);
%     d = m.Data.m;
%     ev_chans = [135:139];
%     ADCs = d(ev_chans,:);
% 
%     for i = 1:size(ADCs,1)
%         test = double(ADCs(i,:));
%         tmp1 = find(test(1:end-1)>max(test)*.75); % if this occasionaly doesnt work we have to switch to prctile
%         tmp2 = find(test(2:end)>max(test)*.75);
%         [both, fall, rise] = setxor(tmp1,tmp2);
%         analogEvents{i} = tmp1(rise); % add names of event channels as per arduino connections and correspondece to the py code conditions
%     end
%     analog_ch_names{5} = 'run_encoder';
%     % overwrite the analog_events.mat
%     save(fullfile(dataPath,fileBase,'analog_events.mat'),'analogEvents','analog_ch_names','user_msgs','ts_user_msgs');
% end
% 
% 
% 
% 
% 
% % we should make a session variable that holds
% 
% session.events =
% session.num_channels = 
% session.SR = 
% session.ADCs = 
% session.speed = 
% session.CA1
% session.V4