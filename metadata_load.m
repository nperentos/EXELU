function [sample_rate, num_channels, channels_description] = metadata_load(full_oebin_path)

% load the .oebin JSON metadata file that comes with each recording
    json_data = jsondecode(fileread(full_oebin_path));
    
% access the continuous stream
    stream_info = json_data.continuous;
    
% see all field names
    field_names = fieldnames(stream_info);
    % for i = 1:length(field_names)
    %     fprintf(field_names{i})
    %     fprintf("\n")
    % end
    
% extract important information
    sample_rate = stream_info.sample_rate;
    num_channels = stream_info.num_channels;
    channels_description = stream_info.channels;

end