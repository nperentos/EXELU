% 27/09/2025
% first script for basic data processing of EXELU data post dat and lfp creation
% KS has already ran at this stage by another script that resides on the recording PC
% the script is called convertDataEXELU.m

% TODO
% check that the spectra t-indices for the different conditions do not intersect
% save the spectra in each relevant folder so we dont have to recompute each time
% add the speed, theta and ripple for filtering epochs
% compute speed and resample at the spectrogram values
% add ripples to remove ripples


if ispc 
    addpath(genpath('C:\Users\npere\MATLAB'));
elseif isunix
    addpath(genpath('/home/nikolas/Documents/MATLAB'))
else
    error('cannot mount code folders -> please resolve...');
end

p = findUsbDiskMountPath("EXELU_SSD1");

if p == ""
    disp("EXELU_SSD1 not found or not mounted.");
else
    fprintf("Data disk found at: %s\n", p);
end



dataPath = p;
fileBase = '2025-09-20_16-51-22';
session_path = fullfile(dataPath,fileBase);
load(fullfile(dataPath,fileBase,[fileBase,'.chanmap.mat']));
load(fullfile(dataPath,fileBase,'analog_events.mat'));
num_channels = length(chanMap);
[data,settings,tScale] = getLFP(fileBase, dataPath);

SR = str2num(settings.parameters.fieldPotentials.lfpSamplingRate.Text);

particulars = loadParticulars(fullfile(session_path,'session_particulars.txt'));

% keep the probe only the following parts

x = single(data(1:128,:)); 
clear data; % just to keep some space

%% Find bad channels - works for linear probe, not tested for multishank!
if ~isfield(particulars,'bad_channels')
    seg = x(:,1:100*SR)';
    seg = seg - median(seg,1);
    C = size(seg,2);
    rho = zeros(C-1,1);
    for i = 1:C-1
        r = corr(seg(:,i), seg(:,i+1), 'type','Spearman'); % Spearman is robust
        rho(i) = r;
    end
    score = nan(C,1);
    score(2:C-1) = (1-rho(1:C-2)) + (1-rho(2:C-1));  % high = suspicious
    z = (score - median(score,'omitnan')) / (1.4826*mad(score,1));
    idx_bad = find(z > 10);
    % compute short segment psd to use for bad channel detection
    movingwin = [1.5 1.5*0.2];
    params.tapers = [2 3];    % Or [2 3] or [3 5] for less smoother spectra
    params.Fs     = SR;       % Sampling rate (Hz)
    params.fpass  = [0 200];   % Focus on theta band
    params.pad    = 0;        % No additional frequency padding
    % params.win    = [1 0.5];  % 1 s window, 0.5 s step (50% overlap)
    [S_seg,t_seg,f_seg]=mtspecgramc(seg,movingwin,params);
    S_seg = sq(mean(S_seg,1));
    lineHz = 50;        % change to 60 if needed
    bw = 1;             % +/- 1 Hz
    lineBand  = f_seg > lineHz-bw & f_seg < lineHz+bw;
    lfpBand   = f_seg > 1 & f_seg < 200;
    linePower = sum(S_seg(lineBand,:),1);
    lfpPower  = sum(S_seg(lfpBand,:),1);
    lineRatio = linePower ./ lfpPower;
    %combine the two metrics
    score = sqrt(max(z,0).* max(zscore(lineRatio)',0));
    bad_chans = find(score > 3);
    disp(['identified bad channels: ', mat2str(bad_chans)]);
    particulars.bad_channels = bad_chans;
    % append the bad channels list to the session_particulars.txt file
    fid = fopen(fullfile(dataPath,fileBase,'session_particulars.txt'),'a');
    %fprintf(fid, ['bad_channels= ',mat2str(bad_chans), '\n']);
    fprintf(fid, 'bad_channels= %s\n', strjoin(string(bad_chans), ' '));
    fclose(fid);
else
    bad_chans = particulars.bad_channels;
end

%% SPECTRA
% theta spectrogram

    movingwin = [1.5 1.5*0.2];
    params.tapers = [3 5];    % Or [2 3] or [3 5] for less smoother spectra
    params.Fs     = SR;       % Sampling rate (Hz)
    params.fpass  = [0 60];   % Focus on theta band
    params.pad    = 0;        % No additional frequency padding
    % params.win    = [1 0.5];  % 1 s window, 0.5 s step (50% overlap)
    theta_ch = particulars.CA1_chan-7;
    tic;
    [S_theta,t_theta,f_theta]=mtspecgramc(x',movingwin,params); %(theta_ch,:) %28 seconds!
    toc
    % figure;
    % for i = 1:num_channels
    %     imagesc(t,f,S(:,:,i)');
    %     axis xy; title(['ch ',num2str(i)]);
    %     pause(0.2); cla;
    % end


% gamma spectrogram
    tic;
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
    toc

    tic;
    [S,t,f]=mtspecgramc(x,movingwin,params); %(theta_ch,:) %28 seconds!
    S = single(S);
    toc %  200 seconds
    % save the spectra to save 200s next time
    % tic;save(fullfile(session_path,'gamma_spectra.mat'),'f','S','t','-v7.3'); toc
% running speed
    runPulses = analogEvents{1,5}./30e3;
    bins = 0:movingwin(1):max(runPulses)+0.1;
    cnts = histcounts(runPulses,bins).*20*pi*9/600; % cm/s
    speed_tmp = smooth(cnts,10); 
    speed = interp1(bins(2:end),speed_tmp,t, 'nearest');
    figure; 
    subplot(3,1,[1,2]); imagesc(t,f,S(:,:,20)'); axis xy; hold on;
    clim([0 150]);
    plot(t,speed.*max(f)./max(speed),'.w','MarkerSize',5,'LineWidth',2);
 


% lets use the mt triggered spectrogram instead
nCh = size(x,2);
% trigs must be in seconds 
test = analogEvents{1,1}./30e3; % events in seconds
test = reshape(test,2,[]);
trigs = test(1,:);
twin = [5, 10];
movingwin = [0.2 0.2*0.2]; % 0.2 s x 20% overlap step
params.tapers = [2 3];    % Or [2 3] for slightly smoother spectra
params.Fs     = SR;     % Sampling rate (Hz)
params.fpass  = [10 200];   % Focus on gamma band
params.pad    = 1;        % No additional frequency padding

% compute triggered spectrograms
for ch = 1:nCh
    [S{ch}, t, f] = mtspecgramtrigc(x(:,ch), trigs, twin, movingwin, params);
    % S is in: t x f x tr
    bsl_mask = (t>=0) & (t<5);
    dur_mask = (t>5);
    % get averages of bsl ts and dur ts. maintain trials
    S_bsl{ch} = sq(mean(S{ch}(bsl_mask,:,:),1));% average of epochs, not trials
    S_dur{ch} = sq(mean(S{ch}(dur_mask,:,:),1));% average of epochs, not trials
    S_diff{ch} = S_dur{ch}-S_bsl{ch};% baselinecorrected trials
    S_mu{ch} = mean(S_diff{ch},2); % mean across trials
    S_sem{ch} = std(S_diff{ch},0,2)./sqrt(length(trigs));% std across trials
end

% plot the RESIDUAL spectra after bsl correction
nRows = 16;
nCols = 8;   % 16 × 8 = 128

figure('Color','w');

t = tiledlayout(nRows, nCols);
t.Padding = 'none';
t.TileSpacing = 'none';

for ch = 1:128
    ax = nexttile;
    % plot(f,mean(S_bsl{ch},2),f,mean(S_dur{ch},2)); hold on;
    
    plot(f,S_mu{ch},'k'); hold on;
    plot(f,S_mu{ch} - S_sem{ch},'color',[0.5 0.5 0.5]);
    plot(f,S_mu{ch} + S_sem{ch},'color',[0.5 0.5 0.5]);
    axis tight
    plot([40 40],ylim,'r--');
    plot(xlim,[0 0],'r--');
    ax.XTick = [];
    ax.YTick = [];
end
title(t,'static grating');











% trialSum = squeeze(sum(S_tmp,[1 2]));
% plot(zscore(trialSum));

figure; 
for i = 1:size(S_tmp,3)
    clf;

    subplot(1,3,[1 2]);
    imagesc(t-5,f,abs(S_tmp(:,:,i)'));
    caxis([0 max(vc(S_tmp))*0.1]);
    hold on;
    plot([0,0],ylim,'k--','LineWidth',2);
    title(['static, trial: ',num2str(i)]);

    subplot(1,3,[3]);
    
    pause;
end

% gather the relevant chunks for each condition. The events are in dat timescale 
% (samples) so lets convert them to seconds

%static
    test = analogEvents{1,1}./30e3; % events in seconds
    test = reshape(test,2,[]);
    % t indices in the spectra matrix demarcating the start of the static exp.
    %idx = arrayfun(@(x) find(abs(t - x) == min(abs(t - x)), 1), test);
    idx = interp1(t, 1:numel(t), test, 'nearest');
    SegLen = median(diff(idx,1)); % should be about 10 s
    M = idx(1,:)' + (0:SegLen);
    %M_static = vc(M');

    % the 5 seconds before
    idx = interp1(t, 1:numel(t), test(1,:), 'nearest');
    nback = length([0:t(2)-t(1):5]);
    N_static = idx' - (0:nback-1);

    % % the 5 seconds before
    % N_static = vc((idx(1,1)'-5/median(diff(t)) + (0:5/median(diff(t))))');
    % % the 10s after
    % N_static = vc((idx(2,:)'+10 + (0:SegLen-20))');


% flashing
    % find diffs in pulse times that exceed 5s
    test = analogEvents{1,2}./30e3;
    
    idx_end = [find(diff(test) > 5) ];
    idx_strt= [1 idx_end(1:end-1)+1];

    idx = interp1(t, 1:numel(t), test(idx_strt), 'nearest');
    M = idx(1,:)' + (0:SegLen);
    M_flashing = vc(M');
    
    % the 10s after
    idx = interp1(t, 1:numel(t), test(idx_end), 'nearest');
    N_flashing = vc((idx'+10 + (0:SegLen-20))');


% reversing
    % find diffs in pulse times that exceed 5s
    test = analogEvents{1,4}./30e3;
    
    idx_end = [find(diff(test) > 5) ];
    idx_strt= [1 idx_end(1:end-1)+1];

    idx = interp1(t, 1:numel(t), test(idx_strt), 'nearest');
    M = idx(1,:)' + (0:SegLen);
    M_reversing = vc(M');

    % the 10s after
    idx = interp1(t, 1:numel(t), test(idx_end), 'nearest');
    N_reversing = vc((idx'+10 + (0:SegLen-20))');

% package the above into a single structure for ease of use
    cond_ind{1,1} = M_static;cond_ind{1,2} = M_flashing;cond_ind{1,3} = M_reversing;
    cond_ind{2,1} = N_static;cond_ind{2,2} = N_flashing;cond_ind{2,3} = N_reversing;

% visual compare same channel different conditions    
    figure; 
    for ch = 22%1:2:128
        subplot(321); imagesc(t(M_static),f   ,S(M_static,:,ch)'); axis xy; title('static grating ');colorbar;
        subplot(323); imagesc(t(M_flashing),f ,S(M_flashing,:,ch)'); axis xy; title('full field flashing at 40 Hz');colorbar;
        subplot(325); imagesc(t(M_reversing),f ,S(M_reversing,:,ch)'); axis xy; title('reversing gratings');colorbar;
        subplot(322); imagesc(t(N_static),f ,S(N_static,:,ch)'); axis xy; title('no stim');colorbar;
        subplot(324); imagesc(t(N_flashing),f ,S(N_flashing,:,ch)'); axis xy; title('no stim');colorbar;
        subplot(326); imagesc(t(N_reversing),f ,S(N_reversing,:,ch)'); axis xy; title('no stim');colorbar;
        mod_all_plots('clim([0 150])');
        pause;%(0.2);clf
    end



%% lets average and get spectra rather than spectrograms for each state and channel.
% figure; 
avgs = [];
for ch = 1:size(S,3) % n of chans

% Static
    S_tmp = S(M_static,:,ch)';
    ttt = sum(S_tmp,1);
    thr = prctile(ttt,95);
    keep = find(ttt < thr);
    S_av = mean(S_tmp(:,keep),2);
    S_sd = std(S_tmp(:,keep),0,2);
    avgs(:,1,ch) = S_av;
    stds(:,1,ch) = S_sd;
    if intersect(ch,bad_chans);  disp(['cj',num2str(ch),' match...']); avgs(:,1,ch) = nan; end
    % plot(f,S_av,'LineWidth',2); hold on;
    
% Flashing
    S_tmp = S(M_flashing,:,ch)';
    ttt = sum(S_tmp,1);
    thr = prctile(ttt,95);
    keep = find(ttt < thr);
    S_av = mean(S_tmp(:,keep),2);
    S_sd = std(S_tmp(:,keep),0,2);
    avgs(:,2,ch) = S_av;
    stds(:,2,ch) = S_sd;
    if intersect(ch,bad_chans);  disp(['cj',num2str(ch),' match...']); avgs(:,2,ch) = nan; end
    % plot(f,S_av,'LineWidth',2); 
    
% Reversing
    S_tmp = S(M_reversing,:,ch)';
    ttt = sum(S_tmp,1);
    thr = prctile(ttt,95);
    keep = find(ttt < thr);
    S_av = mean(S_tmp(:,keep),2); 
    S_sd = std(S_tmp(:,keep),0,2);
    avgs(:,3,ch) = S_av;
    stds(:,3,ch) = S_sd;
    if intersect(ch,bad_chans);  disp(['cj',num2str(ch),' match...']); avgs(:,3,ch) = nan; end
    % plot(f,S_av,'LineWidth',2); 

% No Stimulation
    S_tmp = S([N_static; N_flashing; N_reversing],:,ch)';
    ttt = sum(S_tmp,1);
    thr = prctile(ttt,95);
    keep = find(ttt < thr);
    S_av = mean(S_tmp(:,keep),2); 
    S_sd = std(S_tmp(:,keep),0,2);
    avgs(:,4,ch) = S_av;
    stds(:,4,ch) = S_sd;
    if intersect(ch,bad_chans);  disp(['cj',num2str(ch),' match...']); avgs(:,4,ch) = nan; end
    % plot(f,S_av,'LineWidth',2);    

end
%% PLOTS
% plots depth(channel) vs freq.
figure; 
ttl = {'static','flashing','reversing','no stim'};
for i = 1:4
    plt_data = zscore(sq(avgs(:,i,:)),0,1)';
    subplot(2,2,i); h=imagesc(f,[1:size(S,3)],plt_data);title(ttl{i});xlabel('f (Hz)');ylabel('depth');
    set(gca,'Color',[1 1 1]);  % gray background
    set(h,'AlphaData',~isnan(plt_data));  
    colorbar;
end
mod_all_plots('clim([-5 5])');
% set(gcf,'color','white');
%% sad 
% use the avgs variable to contrast average spectra for V1L4
    ch = 12;
    figure;
    ch = [1,5,12,22,40, 60];
    co = get(groot,'defaultAxesColorOrder');co = [[0,0,0];co];
    colors4 = co([4 5 7 1],:);
    for j = 1:length(ch)
        subplot(2,3,j)
        for i = 1:4
            % plot(f, avgs(:,i,22),'linewidth',1.8); hold on;
            h(i) = shadedErrorBar(f,avgs(:,i,ch(j)),stds(:,i,ch(j))./40,{'linewidth',2,'color',colors4(i,:)},1);
            hold on;
        end
        legend([h(1).mainLine h(2).mainLine h(3).mainLine h(4).mainLine],ttl);
        hold on; 
        xlabel('frequency (Hz)'); ylabel('amplitude (a.u.)');
        set(gcf,'color','default');
        set(gca,'fontsize',16);
        ylim([0 30]);
        plot([40 40], ylim, '--r','LineWidth',2);
        text(max(xlim)/2,max(ylim)*0.9,['ch ',num2str(ch(j))],'FontWeight','bold','FontSize',14);
    end

%%













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