function ablation_figure_beh (force_redo)
%
% Looks at performance and kinematic parameters
%
    if (nargin < 1) ; force_redo = 0 ;end

    types = {'touch','whisking','silent'};

    for t=1:length(types)
        abl_type = types{t};
        anims = get_anims(abl_type);

        frac_corr_overall = [];
        touch_per_trial_overall = [];
        for a=1:length(anims)
            [frac_corr touch_per_trial] = plot_single_animal (anims{a}, abl_type, force_redo);
            frac_corr_overall = [frac_corr_overall ; frac_corr];
            touch_per_trial_overall = [touch_per_trial_overall ; touch_per_trial];
        end
        plot_overall(frac_corr_overall, touch_per_trial_overall, abl_type);
    end

function plot_overall(frac_corr_overall, touch_per_trial_overall, abl_type)
    n_anims = size(frac_corr_overall,1);

    switch abl_type
        case 'touch'
            plot_col = [0 .5 1];

        case 'whisking'
            plot_col = [0 .8 0];

        case 'silent'
            plot_col = 0*[ 1 1 1];
    end

    % plot de stuff
    fh = figure('Position',[0 0 300 300]);

    ax(1) = axes('Position', [.15 .1 .3 .8]);
    hold on;
    for a=1:n_anims
        plot(frac_corr_overall(a,:), 'o-', 'Color', plot_col, 'MarkerSize', 15);
    end
    plot(nanmean(frac_corr_overall), 'o-', 'Color', [0 0 0 ], 'MarkerSize', 10,'MarkerFaceColor',[0 0 0], 'LineWidth',2);
    title(sprintf('%s p=%0.3f', abl_type, ranksum(frac_corr_overall(:,1), frac_corr_overall(:,2))));
    ylabel('Fraction correct');
    axis([0.5 2.5 0 1.1]);
    set(gca,'TickDir','out','FontSize',15, 'XTick',[1 2], 'XTickLabel', {'Pre', 'Post'}, 'YTick', [0 1]);

    ax(2) = axes('Position', [.6 .1 .3 .8]);
    hold on;
    for a=1:n_anims
        plot(touch_per_trial_overall(a,:), 'o-', 'Color', plot_col, 'MarkerSize', 15);
    end
    plot(nanmean(touch_per_trial_overall), 'o-', 'Color', [0 0 0 ], 'MarkerSize', 10,'MarkerFaceColor',[0 0 0], 'LineWidth',2);
    title(sprintf('%s p=%0.3f', abl_type, ranksum(touch_per_trial_overall(:,1), touch_per_trial_overall(:,2))));
    ylabel('Touch per trial');
    maxv = ceil(max(touch_per_trial_overall(:)));
    maxv = 5;
    axis([0.5 2.5 0 maxv]);
    set(gca,'TickDir','out','FontSize',15, 'XTick',[1 2], 'XTickLabel', {'Pre', 'Post'}, 'YTick', [0:maxv]);

    abl_print_fig (fh, ['behavior_overall_abl_type_' abl_type]);


function [frac_corr touch_per_trial] = plot_single_animal(anim, abl_type, force_redo);
    % get data
    settings = get_settings;
    processed_data_path = [settings.processed_data_path filesep 'beh'];
    if (~exist(processed_data_path)) ; mkdir(processed_data_path) ; end
    data_file = sprintf('%s%c%s_%s_behavior_data.mat', processed_data_path, filesep, anim, abl_type);

    if (exist(data_file) & ~force_redo)
        disp(['Loading: ' data_file]);
        load(data_file);
    else
        [pre_path post_path abl_curation_file_path] = get_data_path(anim, abl_type);

        % compile the vectors from all the session files
        [net_dkappa_at_touch_pre theta_at_touch_pre theta_pre frac_correct_pre touches_per_trial_pre] = pull_session_objects(pre_path, anim);
        [net_dkappa_at_touch_post theta_at_touch_post theta_post frac_correct_post touches_per_trial_post] = pull_session_objects(post_path, anim);

        save(data_file, 'net_dkappa_at_touch_pre','theta_at_touch_pre','theta_pre','frac_correct_pre','touches_per_trial_pre', ...
                        'net_dkappa_at_touch_post','theta_at_touch_post','theta_post','frac_correct_post','touches_per_trial_post');
    end

    frac_corr = [frac_correct_pre frac_correct_post];
    touch_per_trial = [touches_per_trial_pre touches_per_trial_post];

    switch abl_type
        case 'touch'
            plot_col_pre = [0 .5 1];
            plot_col_post = [0.5 .75 1];

        case 'whisking'
            plot_col_pre = [0 .8 0];
            plot_col_post = [0.5 1 0.5];

        case 'silent'
            plot_col_pre = 0*[ 1 1 1];
            plot_col_post = 0.5*[1 1 1];
    end

    % plot de stuff
    plot_enabled = 1;
    disp('Plotting of individual variables enabled ; to stop a million plots, find this in the ablation_figure_beh.m and set plot_enabled = 0 above this line');
    if (plot_enabled)
        fh = figure('Position',[0 0 900 300]);
     
        ax(1) = axes('Position', [.1 .2 .2 .6]);
        [n_pre x] = hist(theta_pre, -45:90);
        plot(x, n_pre/sum(n_pre), '-', 'Color', plot_col_pre, 'LineWidth', 2);
        hold on;
        [n_post x] = hist(theta_post, -45:90);
        plot(x, n_post/sum(n_post), '-', 'Color', plot_col_post, 'LineWidth', 2);
        [h pval] = kstest2(n_pre/sum(n_pre), n_post/sum(n_post));
        title(sprintf('%s KS p=%0.3f', anim, pval));
        xlabel('\theta (deg.)');
        aa = axis;
        axis([-45 90 0 aa(4)]);
        set(gca,'TickDir','out','FontSize',15, 'XTick',[-45 0 90], 'YTick', [0 aa(4)]);

        ax(2) = axes('Position', [.4 .2 .2 .6]);
        [n_pre x] = hist(theta_at_touch_pre, -45:90);
        plot(x, n_pre/sum(n_pre), '-', 'Color', plot_col_pre, 'LineWidth', 2);
        hold on;
        [n_post x] = hist(theta_at_touch_post, -45:90);
        plot(x, n_post/sum(n_post), '-', 'Color', plot_col_post, 'LineWidth', 2);
        [h pval] = kstest2(n_pre/sum(n_pre), n_post/sum(n_post));
        title(sprintf('%s KS p=%0.3f', anim, pval));
        xlabel('\theta at touch (deg.)');
        aa = axis;
        axis([-45 90 0 aa(4)]);
        set(gca,'TickDir','out','FontSize',15, 'XTick',[-45 0 90], 'YTick', [0 aa(4)]);

        ax(3) = axes('Position', [.7 .2 .2 .6]);
        [n_pre x] = hist(net_dkappa_at_touch_pre, -.025:.001:.025);
        plot(x, n_pre/sum(n_pre), '-', 'Color', plot_col_pre, 'LineWidth', 2);
        hold on;
        [n_post x] = hist(net_dkappa_at_touch_post, -.025:.001:.025);
        plot(x, n_post/sum(n_post), '-', 'Color', plot_col_post, 'LineWidth', 2);
        [h pval] = kstest2(n_pre/sum(n_pre), n_post/sum(n_post));
        title(sprintf('%s KS p=%0.3f', anim, pval));
        xlabel('Net \Delta\kappa at touch (mm^-^1)');
        aa = axis;
        axis([-.025 .025 0 aa(4)]);
        set(gca,'TickDir','out','FontSize',15, 'XTick',[-.025 0 .025], 'YTick', [0 aa(4)]);

        abl_print_fig (fh, ['behavior_anim_summ_' anim '_abl_type_' abl_type]);
    end


function [net_dkappa_at_touch theta_at_touch theta frac_correct touches_per_trial] = pull_session_objects(dat_path, anim);

    cd (dat_path);
    fl = dir('an*sess.mat'); 

    net_dkappa_at_touch = []; 
    theta_at_touch = [];
    theta = [];

    total_trials = 0;
    total_correct = 0;
    total_touches = 0;

    % loop thru files
    for f=1:length(fl) ; 
        load(fl(f).name) ; 

        theta = [theta s.whiskerAngleTSA.valueMatrix(1,:)];

        net_dkappa_at_touch_tmp = s.whiskerBarContactESA.esa{1}.eventPropertiesHash.get('kappaMaxAbsOverTouch');
        net_dkappa_at_touch = [net_dkappa_at_touch net_dkappa_at_touch_tmp(1:2:end)];

        theta_at_touch_tmp = s.whiskerBarContactESA.esa{1}.eventPropertiesHash.get('thetaAtTouchOnset');
        theta_at_touch = [theta_at_touch theta_at_touch_tmp(1:2:end)];

        tt = sum(s.trialTypeMat');
        total_trials = total_trials + sum(tt(1:4));
        total_correct = total_correct + sum(tt(1:2));

        % intertouch interval of 50 ms 
        [tch_times tch_trials] = s.whiskerBarContactESA.esa{1}.getBurstTimes(pldo.timeSeries.convertTime(50, pldo.timeSeries.millisecond, s.whiskerBarContactESA.esa{1}.timeUnit));
        total_touches = total_touches + length(tch_times);
    end  

    frac_correct = total_correct/total_trials;
    touches_per_trial = total_touches/total_trials;
    
