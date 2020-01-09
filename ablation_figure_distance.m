function ablation_figure_distance
%
% Distance-to-closest-ablated analyses
%

    t_color = [0 0.5 1];
    w_color = [0 0.8 0];

    abl_type = 'touch';
    anims = get_anims(abl_type);

    resp_type = 't'; % touch cells
    fh(1) = single_plot(resp_type, t_color, abl_type, anims, [-0.25 0.15]);
    print_fig (fh(1), ['distance_abl_' abl_type '_resp_' resp_type]);

    resp_type = 'w'; % whisking cells
    fh(2) = single_plot(resp_type, w_color, abl_type, anims, [-0.2 0.2]);
    print_fig (fh(2), ['distance_abl_' abl_type '_resp_' resp_type]);

    abl_type = 'silent'; % this is its own thing -- different
    anims = get_anims(abl_type);
    fh(3) = event_rate_plot([1 1 1]*0, abl_type, anims, [-0.4 0.4],75, [10 25]);
    print_fig (fh(3), ['distance_abl_' abl_type '_event_rate']);

function fh = event_rate_plot(color, abl_type, anims, yrange, max_dist, adjacent_compare_dist)

    if (nargin < 5) ; max_dist = 150 ; end
    if (nargin < 6) ; adjacent_compare_dist = [] ;end

    fh = figure ('Position', [0 0 600 600]); % per-animal delta event rate
    hold on;
    linestyles = {'-',':','--','-.'};

    %% pull data
    dbin = 10;
    dist_bins = 0:dbin:max_dist;
    distances = dist_bins(1:end-1)+dbin/2;
    d_er_by_bin = nan*zeros(length(anims),length(dist_bins)-1);
    all_d_ablated = [];
    all_delta_er = [];
    lii = 1;
    adjacent_data = nan*zeros(length(anims),2);
    n_cells = 0;
    for a=1:length(anims)
        anim = anims{a};
        
        [data data_columns] = get_ablation_data (anim, abl_type, 0);

        ablati_and_proxi = find(data(:,find(strcmp(data_columns, 'proximity_exclude'))));
        d_ablated = data(:,find(strcmp(data_columns, 'd_ablated_um')));
        er_pre = data(:,find(strcmp(data_columns, 'er_pre')));
        er_post = data(:,find(strcmp(data_columns, 'er_post')));
        delta_er = er_post - er_pre;

        all_d_ablated = [d_ablated' all_d_ablated];
        all_delta_er = [delta_er' all_delta_er];

        vali = setdiff(1:size(data,1), ablati_and_proxi);

        % bins
        for d=1:length(dist_bins)-1
            dvali = find(d_ablated >= dist_bins(d) & d_ablated < dist_bins(d+1));
            dvali = intersect(vali,dvali);
            d_er_by_bin(a,d) = nanmedian(delta_er(dvali));
        end
        plot(distances, d_er_by_bin(a,:), linestyles{lii}, 'Color', [1 1 1]*0.75, 'LineWidth', 1);
        lii = lii+1;
        if (lii > length(linestyles)) ; lii = 1; end

        % plot raw data
        plot(d_ablated(vali),delta_er(vali), 'o','Color',[1 1 1]*0,'MarkerFaceColor', [1 1 1]*.85, 'MarkerSize',5);

        % adjacentmal cells
        adjacenti = find(d_ablated > adjacent_compare_dist(1) & d_ablated < adjacent_compare_dist(2));
        adjacenti = setdiff(adjacenti, ablati_and_proxi);
        adjacent_data(a,1) = nanmedian(er_pre(adjacenti));
        adjacent_data(a,2) = nanmedian(er_post(adjacenti));
        n_cells = n_cells + length(adjacenti);
    end
    plot(distances, nanmean(d_er_by_bin), 'Color', color, 'LineWidth', 3);

    % statistical comparison of pre to post for adjacentmal
    if (length(adjacent_compare_dist) == 2)
        adjacenti = find(all_d_ablated > adjacent_compare_dist(1) & all_d_ablated < adjacent_compare_dist(2));
        distali = find(all_d_ablated > adjacent_compare_dist(2) & all_d_ablated < max_dist);

        % this number not in paper, pre/post is given instead
        title(sprintf('distal vs adjacent, p ranksum: %0.3f', ranksum(all_delta_er(adjacenti), all_delta_er(distali))));

        disp(sprintf('abl type: silent ; event rate, adjacent %d - %d', adjacent_compare_dist(1), adjacent_compare_dist(2)));
        disp(sprintf('  adjacent: grand median +/- MAD pre %0.5f +/- %0.3f post %0.5f +/- %0.3f', ...
            nanmedian(adjacent_data(:,1)), 1.4826*mad(adjacent_data(:,1)), nanmedian(adjacent_data(:,2)), 1.4826*mad(adjacent_data(:,2))));
        disp(sprintf('  pval pre v post: %0.3f ; n_anims = %d ; ncells = %d', signrank(adjacent_data(:,1), adjacent_data(:,2)), size(adjacent_data,1), n_cells));
    end

    %% prettify 
    axis([0 max(dist_bins) yrange(1) yrange(2)]);
    plot([0 max(dist_bins)], [0 0],'k-');
    ylabel('Ca^2^+ event rate change (Hz)');
    xlabel('Distance to ablated (\mum)');
    set(gca,'TickDir','out','FontSize',15, 'XTick', [0 max(dist_bins)/2 max(dist_bins)], 'YTick', [yrange(1) 0  yrange(2)]);



function fh = single_plot(resp_type, color, abl_type, anims, yrange)

    thresh = 0.1; % stricter than for core figure

    proxi_dist = [15 35];
    distal_dist = [115 135];

    fh = figure; % per-animal delta
    hold on;
    linestyles = {'-',':','--','-.'};

    %% pull data
    dbin = 25;
    dist_bins = 0:dbin:150;
    distances = dist_bins(1:end-1)+dbin/2;
    d_score_by_bin = nan*zeros(length(anims),length(dist_bins)-1);

    d_score_proximal = [];
    d_score_distal = [];

    lii = 1;
    d_proxi_ani = [];
    d_distal_ani = [];
    for a=1:length(anims)
        anim = anims{a};
        
        [data data_columns] = get_ablation_data (anim, abl_type, 0);
        [score_pre score_post vali ablati] = get_restricted_neurons (data, data_columns, resp_type, thresh, 1);

        d_ablated = data(vali,find(strcmp(data_columns, 'd_ablated_um')));
        d_score = score_post(vali)-score_pre(vali);

        % bins
        for d=1:length(dist_bins)-1
            dvali = find(d_ablated >= dist_bins(d) & d_ablated < dist_bins(d+1));
            d_score_by_bin(a,d) = nanmean(d_score(dvali));
        end
        plot(distances, d_score_by_bin(a,:), linestyles{lii}, 'Color', [1 1 1]*0.75, 'LineWidth', 1);
        lii = lii+1;
        if (lii > length(linestyles)) ; lii = 1; end
    
        d_score_proximal = [d_score_proximal d_score(find(d_ablated > proxi_dist(1) & d_ablated < proxi_dist(2)))'];
        d_score_distal = [d_score_distal d_score(find(d_ablated > distal_dist(1) & d_ablated < distal_dist(2)))'];
    
        d_proxi_ani(a) = nanmean(d_score(find(d_ablated > proxi_dist(1) & d_ablated < proxi_dist(2))));
        d_distal_ani(a) = nanmean(d_score(find(d_ablated > distal_dist(1) & d_ablated < distal_dist(2))));
    end

    plot(distances, nanmean(d_score_by_bin), 'Color', color, 'LineWidth', 3);
    f = fit(distances',nanmean(d_score_by_bin)','exp1');
    cif = confint(f);
    fit_y = f.a*exp(f.b*distances);
    plot(distances, fit_y, '--', 'Color', [0 0 0], 'LineWidth', 2);
    text(50, -0.15, ['\lambda = ' sprintf('%d', round(abs(1/f.b))) '\mum' sprintf(' 95CI %0.2f %0.2f', round((1/cif(1,2))), round((1/cif(2,2))))], 'FontSize', 15);
    disp(sprintf('abl type: %s ; pop: %s, proxi %d - %d ; distal %d - %d', abl_type, resp_type, proxi_dist(1), proxi_dist(2), distal_dist(1), distal_dist(2)));
    disp(sprintf('  proxi: med +/- MAD %0.3f +/- %0.3f distal: %0.3f +/- %0.3f', nanmedian(d_score_proximal), 1.4826*mad(d_score_proximal), nanmedian(d_score_distal), 1.4826*mad(d_score_distal)) );
    disp(sprintf('  proxi vs. distal pval signrank: %0.3f N=%d animals', signrank(d_proxi_ani, d_distal_ani), length(d_proxi_ani)));

    %% prettify 
    axis([0 max(dist_bins) yrange(1) yrange(2)]);
    plot([0 max(dist_bins)], [0 0],'k-');
    if (resp_type == 't')
        ylabel('Touch score change');
    else
        ylabel('Whisking score change');
    end
    xlabel('Distance to ablated (\mum)');
    set(gca,'TickDir','out','FontSize',15, 'XTick', [0 max(dist_bins)/2 max(dist_bins)], 'YTick', [yrange(1) 0  yrange(2)]);



