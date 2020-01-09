function ablation_figure_correlation
%
% This generates the correlation figure (3)
%

    % overall summary panels
    fh(1) = figure('Position', [0 0 600 300]);
    ah(1) = axes('Parent', fh(1), 'Position', [.1 .2 .35 .7]);
    ah(2) = axes('Parent', fh(1), 'Position', [.6 .2 .35 .7]);

    abl_type = 'touch';
    anims = get_anims(abl_type);

    summary_plot(fh(1), ah(1), 't', abl_type, anims);
    summary_plot(fh(1), ah(2), 'w', abl_type, anims);
    abl_print_fig (fh(1), ['correlation_plot_abl_type_' abl_type]);


    fh(2) = figure('Position', [0 0 600 300]);
    ah(3) = axes('Parent', fh(2), 'Position', [.1 .2 .35 .7]);
    ah(4) = axes('Parent', fh(2), 'Position', [.6 .2 .35 .7]);

    abl_type = 'whisking';
    anims = get_anims(abl_type);

    summary_plot(fh(2), ah(3), 't', abl_type, anims);
    summary_plot(fh(2), ah(4), 'w', abl_type, anims);
    abl_print_fig (fh(2), ['correlation_plot_abl_type_' abl_type]);

    % Example animal
    fh(3) = figure('Position', [0 0 600 600]);
    ah(5) = axes('Parent', fh(3), 'Position', [.1 .2 .7 .7]);
    anim = '272761';
    example_plot (fh(3), ah(5), 't','touch',anim);
    title(ah(5), ['Touch ablation, ' anim]);
    abl_print_fig (fh(3), ['correlation_plot_example_' anim]);


function example_plot (fh, ah, resp_type, abl_type, anim);
%
% for the example animal -- you can change the anim/abl_type to get plots for whomever you want
%

    % setup figure
    figure(fh);
    axes(ah);
    hold on;

    % main calls
    force_redo = 0;
    score_bins = -1:.1:1;
   
    % actual data
    plot_single_animal (anim, abl_type, resp_type, 0, ah, 'all_cells', force_redo, score_bins);
    plot_single_animal (anim, abl_type, resp_type, 0, ah, 'average_line_thick', force_redo, score_bins);

    % prettify
    figure(fh);
    axes(ah);
    axis([-1 1 -1 1]);
    axis([-1 1 -1 1]);
    ylabel(['Change in ' resp_type ' score']);
    xlabel('Correlation to ablated');
    title([abl_type ' ablation']);
    set(gca,'TickDir','out','FontSize',15, 'XTick',[-1 0 1], 'YTick',[-1 0 1]);
    plot([-1 1],[0 0],'k:');
    plot([0 0],[-1 1],'k:');


function summary_plot(fh, ah, resp_type, abl_type, anims);
%
% summary plot across multiple animals for single resp/abl type (resp = response we look at - touch or whisking ; abl = who was ablated t,w, or silent)
%
    min_bin_count = 3; % how many values we must have in a bin to take average

    % setup figure
    figure(fh);
    axes(ah);
    hold on;

    % main calls
    force_redo = 0;
    score_bins = -1:.1:1;

    % definitions
    if (resp_type == 't')
        plot_color = [0 .5 1];
        sim_plot_color = [1 1 1]*.75;
    elseif (resp_type == 'w')
        plot_color = [0 .8 0];
        sim_plot_color = [1 1 1]*.75;
    end 

    % actual data
    delta_score_abl = nan*zeros(length(anims),length(score_bins)-1);
    slopes= nan*zeros(length(anims),1);
    n_cells = 0;
    for a=1:length(anims)
        [delta_score_abl(a,:) n_cells_anim slopes(a)] = plot_single_animal (anims{a}, abl_type, resp_type, 0, ah, 'average_line', force_redo, score_bins);
        n_cells = n_cells+n_cells_anim;
    end

    xvals = score_bins(1:end-1)+mode(diff(score_bins))/2;

    nvals = sum(~isnan(delta_score_abl));
    vali = find(nvals >= min_bin_count); % must have @ least n values
    
    plot(xvals(vali), nanmean(delta_score_abl(:,vali)), '-', 'Color', plot_color, 'LineWidth',5);
    [corr_r corr_p] = nancorr(xvals(vali), nanmean(delta_score_abl(:,vali)));
    text(-0.75, -0.25, sprintf('Pearson R: %0.2f p-val: %0.3f n_cells: %d ; pval sign test slope: %0.3f', corr_r, corr_p, n_cells, signtest(slopes)));

    % prettify
    figure(fh);
    axes(ah);
    axis([-1 1 -0.5 0.5]);
    ylabel(['Change in ' resp_type ' score']);
    xlabel('Correlation to ablated');
    title([abl_type ' ablation']);
    set(gca,'TickDir','out','FontSize',15, 'XTick',[-1 0 1], 'YTick',[-0.5 0 0.5]);
    plot([-1 1],[0 0],'k:');
    plot([0 0],[-1 1],'k:');



function [delta_score_per_bin n_cells slope] = plot_single_animal (anim, abl_type, resp_type, simulate_touch_ablation, ah, plot_mode, force_redo, score_bins)
% 
% Generates single animal's figure (or just the data)
%
%   anim: animal id
%   abl_type: type of ablation
%   resp_type: t or w ; which response category to look at?
%   simulate_touch_ablation: if 1, it will simulate touch ablation based on parameters defined here and plot that instead; abl cannot be touch, resp_type must be t
%   ah: axes handle where things are plotted
%   plot_mode: 'average_line', 'all_cells'
%   force_redo: 1 if you want to regeneration correlation matrices (will not apply to get_ablation_data)
%   score_bins: how to divvy up score axis
%
% Will return binned data for this animal so you can do grand mean
%

    %% definitions
    thresh = 0.25;
    thresh_meth = 1;
    min_cells_per_bin = 3;

    %% pull data // generate correlation matrices 
    [data data_columns] = get_ablation_data(anim, abl_type, 0); 
    [score_pre score_post vali ablati] = get_restricted_neurons (data, data_columns, resp_type, thresh, thresh_meth);
    [corr_abl_pre corr_abl_post simulated_ablati] = get_single_animal_data(anim,abl_type, simulate_touch_ablation, data, data_columns,force_redo);

    %% now the core stuff
    delta_score_per_bin = nan*zeros(1,length(score_bins)-1);
    n_cells = 0;
    for b=1:length(score_bins)-1
        subvali = intersect(vali,find(corr_abl_pre > score_bins(b) & corr_abl_pre <= score_bins(b+1)));
        if (length(subvali) > min_cells_per_bin)
            delta_score_per_bin(b) = nanmean(score_post(subvali)-score_pre(subvali));
            n_cells = n_cells+length(subvali);
        end
    end

    %% plot
    t_color = [0 .5 1];
    w_color = [0 .8 0];
    if (resp_type == 'w')
        plot_color = w_color;
        if (simulate_touch_ablation)
            plot_color = [1 1 1]*.75;
        end 
    elseif (resp_type == 't')
        plot_color = t_color;
        if (simulate_touch_ablation)
            plot_color = [1 1 1]*.75;
        end
    end

    slope = nan;
    switch plot_mode
        case 'average_line'
            plot(ah, score_bins(1:end-1)+mode(diff(score_bins))/2, delta_score_per_bin, '-', 'Color', plot_color, 'LineWidth',1);
            xvec= score_bins(1:end-1)+mode(diff(score_bins))/2;
            yvec = delta_score_per_bin;
            vali = find(~isnan(xvec) & ~isnan(yvec));
            pfit = polyfit(xvec(vali), yvec(vali), 1);
            slope = pfit(1);

        case 'average_line_thick'
            plot(ah, score_bins(1:end-1)+mode(diff(score_bins))/2, delta_score_per_bin, '-', 'Color', plot_color, 'LineWidth',5);
            [corr_r corr_p] = nancorr( score_bins(1:end-1)+mode(diff(score_bins))/2, delta_score_per_bin);

            xvec= score_bins(1:end-1)+mode(diff(score_bins))/2;
            yvec = delta_score_per_bin;
            vali = find(~isnan(xvec) & ~isnan(yvec));
            pfit = polyfit(xvec(vali), yvec(vali), 1);
            slope = pfit(1);

            text(-0.75, -0.25, sprintf('Pearson R: %0.2f p-val: %0.3f ; n_cells: %d slope: %0.3f', corr_r, corr_p, n_cells, slope));

        case 'all_cells'
            plot(corr_abl_pre(vali), score_post(vali)-score_pre(vali),'o', 'MarkerFaceColor', t_color, 'Color', [1 1 1]*0.75, 'MarkerSize',8,'LineWidth',1);
    end


function [corr_abl_pre corr_abl_post simulated_ablati] = get_single_animal_data(anim, abl_type, simulate_touch_ablation, data, data_columns, force_redo);

    settings = get_settings;
    processed_data_path = [settings.processed_data_path filesep 'corr'];
    if (~exist(processed_data_path)) ; mkdir(processed_data_path) ; end
    data_file = sprintf('%s%c%s_%s_correlation_data.mat', processed_data_path, filesep, anim, abl_type);
    if (simulate_touch_ablation) ; data_file = strrep(data_file, 'correlation', 'simulated_ablation_correlation') ; end


    if (exist(data_file) & ~force_redo)
        disp(['Loading: ' data_file]);
        load(data_file);
    else
        [pre_path post_path abl_curation_file_path] = get_data_path(anim, abl_type);

        corr_abl_pre = nan*zeros(size(data,1),1);
        corr_abl_post = nan*zeros(size(data,1),1);

        % compile the vectors from all the session files
        [cell_id_pre trial_avg_dff_pre] = pull_session_objects(pre_path, anim);
        [cell_id_post trial_avg_dff_post] = pull_session_objects(post_path, anim);

        % sort the above to make them consistent with "data" matrix
        ids = data(:, find(strcmp(data_columns,'cell_id')));
        ni_pre = 0*ids;
        ni_post = 0*ids;
        for i=1:length(ids)
            ni_pre(i) = find(cell_id_pre == ids(i));
            ni_post(i) = find(cell_id_post == ids(i));
        end
        trial_avg_dff_pre = trial_avg_dff_pre(:,ni_pre);
        trial_avg_dff_post = trial_avg_dff_post(:,ni_post);

        % generate the ablated mean dff -- or simulated
        if (simulate_touch_ablation) % simulated: 20 of top 50 in top 2 subvolumes
            t_pre = data(:,find(strcmp(data_columns,'t_pre')));
            t_pre(find(isnan(t_pre))) = -2;
            subvol = data(:,find(strcmp(data_columns,'subvol')));
            [irr idxi] = sort(t_pre, 'descend');
            [irr ia ib] = intersect(idxi, find(subvol < 3));
            idxi = idxi(ib);
            t_pre(idxi(1:20));
            randi = randperm(50);
            simulated_ablati = idxi(randi(1:20));
            ablati = simulated_ablati;
        else
            ablati = find(data(:,find(strcmp(data_columns,'is_ablated'))));
            simulated_ablati = [];
        end
        trial_avg_ablated = nanmean(trial_avg_dff_pre(:,ablati)')';
        
        % compute correlations to this vector
        for i=1:length(ids)
            corr_abl_pre(i) = nancorr(trial_avg_ablated, trial_avg_dff_pre(:,i));
            corr_abl_post(i) = nancorr(trial_avg_ablated, trial_avg_dff_post(:,i));
        end

        save(data_file, 'corr_abl_pre', 'corr_abl_post','anim','abl_type', 'simulated_ablati');
    end
   
function [cell_id trial_avg_dff] = pull_session_objects(dat_path, anim);

    cd (dat_path);
    fl = dir('an*sess.mat'); 
    ntp_single = 70;
    ntp = ntp_single*2;

    cell_id = [];
    trial_avg_dff = [];

    % loop thru files
    for f=1:length(fl) ; 
        load(fl(f).name) ; 
    
        cell_id = [cell_id s.caTSA.ids];

        % sanity
        if (length( s.caTSA.cellFeatures.get('dffBased_HitLTrialAvg')) == 0)
            gcfParamsBasic.caDffTSA = 'caTSA.dffTimeSeriesArray';
            gcfParamsBasic.rootTag = 'dffBased_';
            gcfParamsBasic.analyses = {'travg'};
            s.generateCellFeaturesHash('local', gcfParamsBasic);       
        end

        % compile trial_avg_dff
        this_matrix =  nan*zeros(ntp, length(s.caTSA.ids));
        vv = s.caTSA.cellFeatures.get('dffBased_HitLTrialAvg');
        this_matrix(1:ntp_single,:) = vv(1:ntp_single,:);
        vv = s.caTSA.cellFeatures.get('dffBased_HitRTrialAvg');
        this_matrix(1:ntp_single,:) = vv(1:ntp_single,:);

        trial_avg_dff = [trial_avg_dff this_matrix];
    end 
