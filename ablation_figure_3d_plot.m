% will print pre, ablated (orange), and post
function ablation_figure_3d_plot (anim, ablation_type)
    t_color = [0 0.5 1];
    w_color = [0 0.8 0];

    % touch ablation
    abl_type = 'touch';
    val_range = [0.1 0.4];
    anim = '271211';
    resp_type = 't';
    fh(1) = plot_single_ablation(anim, abl_type, resp_type, t_color, val_range);
    abl_print_fig (fh(1), ['3d_plot_' resp_type '_map_' abl_type '_ablation_anim_' anim]);
    resp_type = 'w';
    fh(2) = plot_single_ablation(anim, abl_type, resp_type, w_color, val_range);
    abl_print_fig (fh(2), ['3d_plot_' resp_type '_map_' abl_type '_ablation_anim_' anim]);

    % w ablation
    abl_type = 'whisking';
    val_range = [0.05 0.35];
    anim = '281915';
    resp_type = 't';
    fh(1) = plot_single_ablation(anim, abl_type, resp_type, t_color, val_range);
    abl_print_fig (fh(1), ['3d_plot_' resp_type '_map_' abl_type '_ablation_anim_' anim]);
    resp_type = 'w';
    fh(2) = plot_single_ablation(anim, abl_type, resp_type, w_color, val_range);
    abl_print_fig (fh(2), ['3d_plot_' resp_type '_map_' abl_type '_ablation_anim_' anim]);

    % si ablation
    abl_type = 'silent';
    val_range = [0.05 0.3];
    anim = '278939';
    resp_type = 't';
    fh(1) = plot_single_ablation(anim, abl_type, resp_type, t_color, val_range);
    abl_print_fig (fh(1), ['3d_plot_' resp_type '_map_' abl_type '_ablation_anim_' anim]);
    resp_type = 'w';
    fh(2) = plot_single_ablation(anim, abl_type, resp_type, w_color, val_range);
    abl_print_fig (fh(2), ['3d_plot_' resp_type '_map_' abl_type '_ablation_anim_' anim]);

function fh = plot_single_ablation(anim, abl_type, resp_type, color, val_range)
    abl_color = [1 0.5 0];
    z0_um = 128; % L1-L2 boundary

    [data data_columns] = get_ablation_data (anim, abl_type);
    [score_pre score_post vali ablati] = get_restricted_neurons (data, data_columns, resp_type, val_range(1), 1);
    xyz = nan*zeros(length(score_pre), 3);
    xyz(:,1) = data(:,find(strcmp(data_columns, 'x_um')));
    xyz(:,2) = data(:,find(strcmp(data_columns, 'y_um')));
    xyz(:,3) = data(:,find(strcmp(data_columns, 'z_um'))) + z0_um;

    fh = figure('Position', [0 0 1500 500]);
    ax(1) = axes('Position', [.1 .1 .1875 .66]);
    ax(2) = axes('Position', [.3 .1 .1875 .66]);
    ax(3) = axes('Position', [.5 .1 .1875 .66]);
    ax(4) = axes('Position', [.7 .1 .1875 .66]);
    plot_single_panel(xyz, score_pre, vali, color, ax(1), val_range(2));
    plot_single_panel(xyz, score_pre, ablati, color, ax(2), val_range(2),0 );
    plot_single_panel(xyz, score_post, vali, color, ax(3), val_range(2), 0);
    plot_single_panel(xyz, score_pre, ablati, abl_color, ax(4), val_range(2), 0, 1);

function plot_single_panel(xyz, score, vali, color, ax, sat_val, labels_shown, highlight_vali)
    size_range = [1 15]; % 2 15
    min_col_frac = 0.6; % color will be at least this fraction of its max
    if (nargin < 7) ; labels_shown = 1 ; end
    if (nargin < 8) ; highlight_vali = 0 ; end

    axes(ax);
    cla;
    hold on;

    score = score/sat_val;
    score(find(score > 1)) = 1;
    score(find(score < 0)) = 0;

    if (highlight_vali) % special case used to outline ablated cells and show range of scores
        for v=1:length(vali)
            i = vali(v);
            plot3(xyz(i,1),xyz(i,2),xyz(i,3),'o','Color',color,'MarkerFaceColor',[1 1 1],'MarkerSize',20);
        end
        plot3(100,100,200,'o','Color',[0 0 0],'MarkerFaceColor',[0 0 0], 'MarkerSize',size_range(2));
        plot3(100,100,300,'o','Color',[0 0 0],'MarkerFaceColor',[0 0 0], 'MarkerSize',size_range(1));
    else
        % non-categorical
        invali = setdiff(1:length(xyz), vali);
        plot3(xyz(invali,1),xyz(invali,2),xyz(invali,3),'.','Color', [1 1 1]*0.5, 'MarkerSize', 1);

        % categorical
        for v=1:length(vali)
            i = vali(v);
              
            col = color*max(min_col_frac,score(i));
            col = min(col,[1 1 1]);
            msize = round(max(size_range(1), score(i)*size_range(2)));
            msize = round(max(size_range(1), score(i)*score(i)*size_range(2)));

            plot3(xyz(i,1),xyz(i,2),xyz(i,3),'o','Color',col,'MarkerFaceColor',col,'MarkerSize',msize);
        end
    end

    panel_color = [1 1 1]*0.95;
    fill3([0 600 600 0 ], [0 0 600 600], [ 1 1 1 1 ]*400, panel_color, 'EdgeColor', panel_color);

    set(gca,'XTick',[],'YTick',[],'YDir','reverse', 'ZDir','reverse','FontSize',15);
    set(gca,'CameraPosition', [-3000 4000 -750]); % was -3000 4000 -500
    ax = gca;
    axis([0 600 0 600 100 400]);
    set(ax.YAxis,'visible','off');
    set(ax.XAxis,'visible','off');

    if (labels_shown)
        ax.ZTick= [100 200 300 400];
        ax.ZLabel.String = 'Depth \mum';
        ax.Title.String = sprintf('Two whisker 3d plot, n=%d cells', length(score));
    else
        set(ax.ZAxis,'visible','off');
        ax.ZTick= [];
    end

