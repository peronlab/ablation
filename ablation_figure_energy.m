function ablation_figure_energy
%
% Effectiveness of ablation vs depth and energy
%
    force_redo = 0; % reload everything from rawish data
    settings = get_settings;
    processed_data_path = [settings.processed_data_path filesep 'energy'];
    if (~exist(processed_data_path)) ; mkdir(processed_data_path) ; end

    abl_types = {};
    abl_type = 'touch';
    anims = get_anims(abl_type);
    for a=1:length(anims) ; abl_types{end+1} = abl_type ;end
    abl_type = 'silent';
    anims_t = get_anims(abl_type);
    for a=1:length(anims_t) ; abl_types{end+1} = abl_type ; anims{end+1} = anims_t{a}; end
    abl_type = 'whisking';
    anims_t = get_anims(abl_type);
    for a=1:length(anims_t) ; abl_types{end+1} = abl_type ; anims{end+1} = anims_t{a}; end

    % Gather data
    depths = [];
    net_energy_joule = [];
    max_power_watt = [];
    success = [];
    num_attempts_per_anim = zeros(1,length(anims));
    num_success_per_anim = zeros(1,length(anims));
    depth_bins = 0:25:400;

    for a=1:length(anims)

        % save locally 
        data_file = sprintf('%s%c%s_%s_ablation_energy_data.mat', processed_data_path, filesep, anims{a}, abl_types{a});
        if (exist(data_file, 'file') & ~force_redo)
            disp(['Loading ' data_file '; to force regeneration, set force_redo to 1']);
            load(data_file);
        else
            [depth_cell success_single max_power_cell_watt net_energy_cell_joule] = get_single_animal_data(anims{a}, abl_types{a});
            save(data_file, 'depth_cell','success_single','max_power_cell_watt','net_energy_cell_joule');
        end

        depths = [depths depth_cell];
        success = [success success_single];
        max_power_watt = [max_power_watt max_power_cell_watt];
        net_energy_joule = [net_energy_joule net_energy_cell_joule];

        num_attempts_per_anim(a) = length(depth_cell);
        num_success_per_anim(a) = sum(success_single);
    end

    % Who got how many?  sanity check against cell count table
    for a=1:length(anims)
        disp(sprintf('ID: %s %s success/attempt: %03d/%03d', anims{a}, abl_types{a}, num_success_per_anim(a), num_attempts_per_anim(a)));
    end

    successi = find(success);
    faili = find(~success);

    frac_success = nan*depth_bins;
    max_power = nan*depth_bins;
    net_energy = nan*depth_bins;
    attempt_count = nan*depth_bins;
    for d=1:length(depth_bins)-1
        depthi = find(depths >= depth_bins(d) & depths < depth_bins(d+1));
        frac_success(d) = sum(success(depthi))/length(depthi);
        attempt_count(d) = length(depthi);

        depthi = intersect(successi,depthi);
        max_power(d) = nanmean(max_power_watt(depthi));
        net_energy(d) = nanmean(net_energy_joule(depthi));
    end

    ms_small = 7;
    ms_large = 15;
    raw_color = [1 1 1]*0.5;

    % --- 1) peak power v. depth
    fh(1) = figure('Position',[0 0 400 400]);
    plot(depths(successi), max_power_watt(successi), 'o', 'Color', [0 0 0], 'MarkerSize', ms_small, 'MarkerFaceColor', raw_color);
    hold on;
    plot(depth_bins+mode(diff(depth_bins)), max_power, 'o', 'Color', [0 0 0], 'MarkerSize', ms_large, 'MarkerFaceColor', [0 0 0]);
    axis([0 depth_bins(end) 0 1.5]);

    % prettify
    ylabel('Peak power (W)');
    xlabel('Depth \mum');
    set(gca,'TickDir','out','FontSize',15, 'XTick',0:100:400, 'YTick',[0 0.5 1 1.5]);
    abl_print_fig (fh(1), 'method_peak_power_v_depth');

    % --- 2) net energy v. depth
    fh(2) = figure('Position',[0 0 400 400]);
    plot(depths(successi), net_energy_joule(successi), 'o', 'Color', [0 0 0], 'MarkerSize', ms_small, 'MarkerFaceColor',raw_color);
    hold on;
    plot(depth_bins+mode(diff(depth_bins)), net_energy, 'o', 'Color', [0 0 0], 'MarkerSize', ms_large, 'MarkerFaceColor', [0 0 0]);
    axis([0 depth_bins(end) 0 10]);

    % prettify
    ylabel('Total energy (J)');
    xlabel('Depth \mum');
    set(gca,'TickDir','out','FontSize',15, 'XTick',0:100:400, 'YTick',[0 5 10]);
    abl_print_fig (fh(2), 'method_net_energy_v_depth');

    % --- 3) success v. depth
    fh(3) = figure('Position',[0 0 400 400]);
    plot(depth_bins+mode(diff(depth_bins)), frac_success, 'o', 'Color', [0 0 0], 'MarkerSize', ms_large, 'MarkerFaceColor', [0 0 0]);
    for d=1:length(depth_bins)-1
        if(~isnan(frac_success(d)))
            text(depth_bins(d)+mode(diff(depth_bins)), frac_success(d)-0.1, sprintf('%d',attempt_count(d)));
        end
    end
    axis([0 depth_bins(end) 0 1.1]);

    % prettify
    ylabel('Fraction success');
    xlabel('Depth \mum');
    title(sprintf('Total attempts: %d ; success: %d ; fail: %d', length(successi)+length(faili), length(successi), length(faili)));
    set(gca,'TickDir','out','FontSize',15, 'XTick',0:100:400, 'YTick',[0 0.5 1]);
    abl_print_fig (fh(3), 'method_frac_success');




function [depth_cell success max_power_cell_watt net_energy_cell_joule] = get_single_animal_data(anim, abl_type)

    z0_um = 128;
    max_y = 600;

    data_type = 'jrc';
    nyu_anims = {'276013','276015','275801','275798','274577','274424','272761','278288','278759'};
    if (ismember(anim, nyu_anims)) ; data_type = 'nyu'; end

    [pre_path post_path abl_curation_file_path] = get_data_path(anim, abl_type);
    if (iscell(abl_curation_file_path) & length(abl_curation_file_path) > 1) ; abl_curation_file_path = abl_curation_file_path{1} ; end 
    [data dc] = get_ablation_data(anim,abl_type,0);
    cell_id_data = data(:,find(strcmp(dc,'cell_id')));
    y_data = data(:,find(strcmp(dc,'y_um')));
    z_offs = 15*(y_data/max_y); 
    depth_data = data(:,find(strcmp(dc,'z_um'))) + z0_um + z_offs;

    %% collect data -- all ablations, on a per-cell basis
    switch data_type
        case 'nyu'
            max_power_watt = 1; % 1 W max
            abl_path = fileparts_crossplatform(abl_curation_file_path);
            fl = dir([abl_path filesep 'session*attempt*mat']);
            cell_id = zeros(1,length(fl));
            net_energy_attempt_joule = zeros(1,length(fl));
            max_power_attempt_watt = zeros(1,length(fl));

            % per attempt data gather
            for f=1:length(fl)
                uidx = find(fl(f).name == '_');
                cell_id(f) = str2num(fl(f).name(uidx(4)+1:uidx(5)-1));

                load([abl_path filesep fl(f).name]);
                max_time = max(ablationLog.rawFluoTimeVecSec);
                if (isempty(max_time)) ; max_time = 0 ;end

                maxi = max(find(ablationLog.powerAndShutterTimeVecSec <= max_time ));

                max_power_attempt_watt(f) = max(ablationLog.ablationPowerVec(1:maxi)/100)*max_power_watt;

                vali = find(ablationLog.ablationPowerVec > 0);
                vali = vali(find(vali <= maxi));

                dt_sec = mode(diff(ablationLog.powerAndShutterTimeVecSec));
                net_time = length(vali)*dt_sec;
                net_energy_attempt_joule(f) = mean(ablationLog.ablationPowerVec(vali)/100)*max_power_watt*net_time;
            end

            md = load(abl_curation_file_path);
            successful_id = unique(md.ablationTargetIds(find(md.ablationSuccessful)));
            
        case 'jrc'
            max_power_watt = 1; % 1 W max

            % per attempt
            abl_path = fileparts_crossplatform(abl_curation_file_path);
            subname = fileparts_crossplatform(abl_curation_file_path,2);
            uidx = find(subname == '_');
            abl_path = [abl_path filesep 'ablation_logs' filesep subname(uidx(1)+1:end)];

            % per attempt data gather
            fl = dir([abl_path filesep 'ablation_log*mat']);
            disp(['Getting log files from : ' abl_path filesep 'ablation_log*mat']);
            cell_id = nan*zeros(1,length(fl));
            net_energy_attempt_joule = nan*zeros(1,length(fl));
            max_power_attempt_watt = nan*zeros(1,length(fl));
            for f=1:length(fl)
                disp(['Loading ' abl_path filesep fl(f).name]);
                load([abl_path filesep fl(f).name]);
                if (isfield(la, 'ablatedId'))
                    cell_id(f) = la.ablatedId;
                    max_time = la.pmtTimeVec(la.endIdx);
                    maxi = max(find(la.beamTime <= max_time ));

                    max_power_attempt_watt(f) = max(la.beamBuffer(1:maxi))*max_power_watt;

                    vali = find(la.beamBuffer > la.beamBuffer(1)+.01);
                    vali = vali(find(vali <= maxi));

                    dt_sec = mode(diff(la.beamTime));
                    net_time = length(vali)*dt_sec;
                    net_energy_attempt_joule(f) = mean(la.beamBuffer(vali))*max_power_watt*net_time  ; 
                else
                    disp('  No field ablatedId; this could be a problem -- check the log file loaded before this statement');
                end
            end

            vali = find(~isnan(cell_id));
            cell_id = cell_id(vali);
            net_energy_attempt_joule = net_energy_attempt_joule(vali);
            max_power_attempt_watt = max_power_attempt_watt(vali);

            md = load(abl_curation_file_path);
            successful_id = md.ablated_ids;
           
    end

    %% convert to per-cell data
    u_cell_id = unique(cell_id);
    net_energy_cell_joule = 0*u_cell_id;
    max_power_cell_watt = 0*u_cell_id;
    success = 0*u_cell_id;
    depth_cell = 0*u_cell_id;
    for c=1:length(u_cell_id)
        fi = find(cell_id == u_cell_id(c));

        net_energy_cell_joule(c) = sum(net_energy_attempt_joule(fi));
        max_power_cell_watt(c) = max(max_power_attempt_watt(fi));
        
        % find in data
        datai = find(cell_id_data == u_cell_id(c)); 
        depth_cell(c) = depth_data(datai);

        if (ismember(u_cell_id(c), successful_id)) ; success(c) = 1; end
    end


