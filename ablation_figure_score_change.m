function ablation_figure_score_change (force_redo)
%
% Will compute score change for touch and whisking populations for all abl types ; set force_redo = 1 to regather data from raw
%
    if (nargin < 1) ; force_redo = 0 ;end
    ablation_figure_score_change_single ('t', force_redo);
    ablation_figure_score_change_single ('w', force_redo);

function ablation_figure_score_change_single (resp_type, force_redo)
%
% Will make the change figures for a given response type (t = touch ; w = whisking) for all ablations (t,w,si)
%
    t_color = [0 0.5 1];
    w_color = [0 0.8 0];
    si_color = [0 0 0 ];

    fh(1) = figure; % per-animal delta
    hold on;
    fh(2) = figure; % delta as a function of net t killed
    hold on;
    fh(3) = figure; % pre vs post by type / animal
    hold on;

    abl_type = 'touch';
    anims = get_anims(abl_type);
    [delta_t score_abl all_deltas_t] = single_type_basic (abl_type, anims, force_redo, fh, 1, t_color, resp_type);
    score_abl_t_only = score_abl;
    overall_t_abl = score_abl;

    abl_type = 'silent';
    anims = get_anims(abl_type);
    [delta_si score_abl all_deltas_si] = single_type_basic (abl_type, anims, force_redo, fh, 2 , si_color, resp_type);
    overall_t_abl = [overall_t_abl score_abl];

    abl_type = 'whisking';
    anims = get_anims(abl_type);
    [delta_w score_abl all_deltas_w]= single_type_basic (abl_type, anims, force_redo, fh, 3 , w_color, resp_type);
    overall_t_abl = [overall_t_abl score_abl];

    %% prettify 
    if (resp_type == 'w') ; y_label = 'Whisking score' ; else ; y_label = 'Touch score'; end

    figure(fh(1))
    axis([0 4 -0.15 0.1])
    plot([0 4],[0 0],'k:');
    ylabel([y_label ' change']);
    xlabel('Ablation type');

    disp(' ' )
    disp('Statistics: ' )
    tsi = ranksum(delta_t, delta_si);
    tw = ranksum(delta_t, delta_w);
    wsi = ranksum(delta_w, delta_si);
    disp(sprintf('T-Si: %0.3f T-W: %0.3f W-Si: %0.3f', tsi, tw, wsi));
    t0 = signrank(delta_t);
    si0 = signrank(delta_si);
    w0 = signrank(delta_w);
    disp(sprintf('T-0: %0.3f Si-0: %0.3f W-0: %0.3f', t0, si0, w0));

    set(gca,'TickDir','out','FontSize',15, 'XTick',[1 2 3], 'XTickLabel', {'Touch', 'Silent', 'Whisking'});
    text(.1,-.1, sprintf('T-Si: %0.3f T-W: %0.3f W-Si: %0.3f', tsi, tw, wsi));

    abl_print_fig(fh(1), ['ablation_figure_score_change_' resp_type]); 

    if (resp_type == 't')
        figure(fh(2));
        aa = axis;
        axis([aa(1) aa(2) -0.15 0.1]);
        plot([aa(1) aa(2)],[0 0],'k:');
        xlabel('Net touch encoding score removed');
        ylabel('Touch score change');
  
        [r p] = nancorr(overall_t_abl, [delta_t delta_si delta_w]);
        text(6,.05,sprintf('Pearson R: %0.3f p: %0.3f', r, p));
        [r p] = nancorr(score_abl_t_only, delta_t);
        text(6,.075,sprintf('T-only R: %0.3f p: %0.3f', r, p));

        set(gca,'TickDir','out','FontSize',15);

        abl_print_fig(fh(2), ['ablation_figure_change_versus_ablation_size_' resp_type]); 
    else
        close(fh(2));
    end

    figure(fh(3))
    axis([0 4 0 0.2])
    ylabel(y_label);
    xlabel('Ablation type');
    set(gca,'TickDir','out','FontSize',15, 'XTick',[1 2 3], 'XTickLabel', {'Touch', 'Silent', 'Whisking'});
    abl_print_fig(fh(3), ['ablation_figure_comparison_' resp_type]); 


function [deltas score_ablated all_deltas] = single_type_basic (abl_type, anims, force_redo, fh, avg_n, avg_color, resp_type)
%
% makes a single plot for a single ablation specified by abl_type (si, w, t) and animal (anims)
%
    anim_labels = 0;
    thresh = 0.1;
    thresh_meth = 1;
    t_color = [0 0.5 1];
    w_color = [0 0.8 1];

    score_pre_overall = [];
    score_post_overall = [];
    score_ablated = [];

    all_deltas = [];

    disp(' ');
    disp('=================================================================================');
    disp(['  Ablation type: ' abl_type ' Response category: ' resp_type]);

    for a=1:length(anims)
        anim = anims{a};
        
        [data data_columns] = get_ablation_data (anim, abl_type, force_redo, 1);
        [score_pre score_post vali ablati] = get_restricted_neurons (data, data_columns, resp_type, thresh, thresh_meth);

        % pre vs post
        score_pre_overall(a) =  nanmedian(score_pre(vali));
        score_post_overall(a) =  nanmedian(score_post(vali));
        score_ablated(a) = nansum(score_pre(ablati)); 

        all_deltas = [all_deltas score_post(vali)'-score_pre(vali)'];

        disp(sprintf('%s Total cells: %d number meeting criteria: %d (frac: %0.3f) number ablated: %d pre score: %0.3f post score: %0.3f delta: %0.3f', ...
                      anim, size(data,1), length(vali), length(vali)/size(data,1), length(ablati), score_pre_overall(a), score_post_overall(a),  score_post_overall(a)-score_pre_overall(a)));
        if (0) % individual scatter plot 
            figure ; 
            plot(score_pre(vali),score_post(vali),'o', 'Color', avg_color) ; 
            hold on ; 
            plot([-.2 .9],[-.2 .9], 'k-'); 
            title(sprintf('%s pre: %0.3g post: %0.3g ablated: %0.3g n_a_b_l=%d n=%d', ...
                  anim, nanmedian(score_pre(vali)), nanmedian(score_post(vali)), score_ablated(a), length(ablati), length(vali) ));
            axis square;
        end

    end

    deltas = score_post_overall-score_pre_overall;

    % change per-animal
    figure(fh(1));
    x_offs = rand(1,a)/8-0.0625;
    plot(avg_n, nanmean(score_post_overall-score_pre_overall), 'o', 'Color', avg_color,'MarkerSize',12,'MarkerFaceColor',avg_color) ;
    single_color = avg_color;
    single_color(find(avg_color== 0)) = 0.5;
    single_color(find(avg_color== 0.5)) = 0.75;
    plot(0*score_pre_overall+avg_n+x_offs, score_post_overall-score_pre_overall, 'o', 'Color', single_color,'MarkerSize', 8,  'LineWidth',1) ;
    if (anim_labels)
        for a=1:length(anims)
            text(avg_n+x_offs(a)+0.25,  score_post_overall(a)-score_pre_overall(a), anims{a});
        end
    end

    % change vs ablated
    figure(fh(2))
    plot(score_ablated, deltas, 'o', 'Color', [1 1 1]*0.8, 'MarkerFaceColor', avg_color, 'MarkerSize',10);

    
    % per-animal pre/post
    figure(fh(3));
    for a=1:length(score_post_overall)
        plot(avg_n+[-0.25 0.25], [score_pre_overall(a) score_post_overall(a)], 'o', 'Color', avg_color,'MarkerSize',10);
        plot(avg_n+[-0.25 0.25], [score_pre_overall(a) score_post_overall(a)], '-', 'Color', avg_color, 'LineWidth', 0.5);
        if (anim_labels) ; text(avg_n-0.5, score_pre_overall(a), anims{a}); end
    end
    text(avg_n-0.25, max([score_pre_overall score_post_overall])+0.02, sprintf('P=%0.3f, signrank', signrank(score_pre_overall,score_post_overall)));

