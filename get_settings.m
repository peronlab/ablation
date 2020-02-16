function settings = get_settings
%
% This should be edited to point to the data and where intermediate files should go
%
%   settings.data_path: where data lives (top-level directory)
%   settings.processed_data_path: summary output files that are generated from the "raw" data so that things go faster on subsequent runs
%
    % data_path: where the data lives (downloaded from CRCNS)    
    data_path = '/Volumes/proc_ablate/crcns_upload_for_manuscript/';
    if (~exist(data_path,'dir'))
        error([data_path ' is not a valid directory; see get_settings.m for this']);
    end
    settings.data_path = data_path;

    % processed_data_path
    processed_data_path = '~/Desktop/ablation_data';
    if (~exist(processed_data_path,'dir'))
        disp([processed_data_path ' is not a valid directory; see get_settings.m to change this, or hit any key to create']);
        pause;
        mkdir(processed_data_path);
    end
    settings.processed_data_path = processed_data_path;

    % print_dir and print_enabled -- if print_enabled, it will print any output figures to print_dir
    print_dir = '~/Desktop/ablation_figs';
    settings.print_dir = print_dir;
    settings.print_on = 0;


