function parameters=neuroscope_xml_creator_EXELU(filebase,SR,nChannels,lfpSR, mapping)
%s = xml2struct([filebase '.settings.xml']);
%SR = s.Settings.AcquisitionSettings.Attributes.SR;
%nChannels = str2num(s.Settings.AcquisitionSettings.Attributes.nChannels);
%lfpSR = s.Settings.ProcessingSettings.Attributes.lfpSR;

% working off of the assumption that are using single shank probes always
multicol = 0;

parameters.Attributes.creator='neuroscope-1.3.3';
parameters.Attributes.version='1.0';
parameters.acquisitionSystem.nBits = '16';
parameters.acquisitionSystem.nChannels =num2str(nChannels);
parameters.acquisitionSystem.samplingRate = SR;
parameters.acquisitionSystem.voltageRange = '10';
parameters.acquisitionSystem.amplification = '1';
parameters.acquisitionSystem.offset = '0';
parameters.fieldPotentials.lfpSamplingRate = lfpSR;
parameters.files.file.samplingRate = lfpSR;
parameters.files.file.extension = 'lfp';

% lets use mapping to define anatomical groups bit2volts colors etc

tpe = [mapping.type];
bit2volt = [mapping.bit_volts];
grps = unique(tpe);

for i = 1: length(grps)
    ch_in_grp = find(tpe == grps(i));
    for j = 1:length(ch_in_grp)
        ch = ch_in_grp(j);
        parameters.anatomicalDescription.channelGroups.group{i}.channel{j}.Text = num2str(ch-1);  % Channel index
        parameters.anatomicalDescription.channelGroups.group{i}.channel{j}.Attributes.skip = '0';  % Optional
        parameters.anatomicalDescription.channelGroups.group{i}.channel{j}.ADBitVolts.Text = num2str(bit2volt(ch), '%.12f');
    end
end

% if nargin>4 & ~isempty(mapping);
%     switch mapping
%         case 'Buzsaki64';
%             multicol=1;
%             for c=0:7
%                 for ch=1:8;
%                     parameters.anatomicalDescription.channelGroups.group{c+1}.channel{ch}.Attributes.skip = '0';
%                     parameters.anatomicalDescription.channelGroups.group{c+1}.channel{ch}.Text = num2str(c*8+ch-1);
%                 end
%             end
%         case 'Buzsaki64sp';
%             multicol=1;
%             for c=0:7
%                 for ch=1:8;
%                     parameters.anatomicalDescription.channelGroups.group{c+1}.channel{ch}.Attributes.skip = '0';
%                     parameters.anatomicalDescription.channelGroups.group{c+1}.channel{ch}.Text = num2str(c*8+ch-1);
%                 end
%             end
%         case 'A8x8';
%             multicol=1;
%             for c=0:7
%                 for ch=1:8;
%                     parameters.anatomicalDescription.channelGroups.group{c+1}.channel{ch}.Attributes.skip = '0';
%                     parameters.anatomicalDescription.channelGroups.group{c+1}.channel{ch}.Text = num2str(c*8+ch-1);
%                 end
%             end
%         case 'A4x8';
%             multicol=1;
%             for c=0:3
%                 for ch=1:8;
%                     parameters.anatomicalDescription.channelGroups.group{c+1}.channel{ch}.Attributes.skip = '0';
%                     parameters.anatomicalDescription.channelGroups.group{c+1}.channel{ch}.Text = num2str(c*8+ch-1);
%                 end
%             end
%         otherwise
%             multicol=0;          
%     end
% end
% if nargin<5 || isempty(mapping) || multicol==0;
%     for k=0:(nChannels-1)
%         parameters.anatomicalDescription.channelGroups.group.channel{k+1}.Attributes.skip = '0';
%         parameters.anatomicalDescription.channelGroups.group.channel{k+1}.Text = num2str(k);
%     end
% end

parameters.spikeDetection='';
parameters.neuroscope.Attributes.version='1.3.3';
parameters.neuroscope.miscellaneous.screenGain='0.01';
parameters.neuroscope.miscellaneous.traceBackgroundImage='';
parameters.neuroscope.video.rotate='0';
parameters.neuroscope.video.flip='0';
parameters.neuroscope.video.videoImage='0';
parameters.neuroscope.video.positionsBackground='0';
parameters.neuroscope.spikes.nSamples='32';
parameters.neuroscope.spikes.peakSampleIndex='16';
for k=0:(nChannels-1)
    parameters.neuroscope.channels.channelColors{k+1}.channel=num2str(k);
    parameters.neuroscope.channels.channelColors{k+1}.color='#0080ff';
    parameters.neuroscope.channels.channelColors{k+1}.anatomyColor='#0080ff';
    parameters.neuroscope.channels.channelColors{k+1}.spikeColor='#0080ff';
    parameters.neuroscope.channels.channelOffset{k+1}.channel=num2str(k);
    parameters.neuroscope.channels.channelOffset{k+1}.defaultOffset='0';
end
paramsxml.parameters=parameters;
struct2xml(paramsxml,[filebase '.xml']);