% assorted stats from the paper
function ablation_summ_stats
    force_redo = 0;

    all_anims = {};
    is_first = [];
    n_cells = [];
    n_t_cells = []; 
    n_w_cells = []; 
    n_ablated = [];
    type_ablated = [];
    score_t_all = {};
    score_w_all = {};
    vali_t_all = {};
    vali_w_all = {};
    ablati_all = {};

    abl_type = 'touch';
    anims = get_anims(abl_type);
    all_anims = [all_anims anims];
    first_anims = {'250220','257220','274424','272761','274577'};
    is_first_this = zeros(1,length(anims)); 
    is_first_this(find(ismember(anims,first_anims))) = 1;
    is_first = [is_first is_first_this];
    [n_cells_this n_t_cells_this n_ablated_this score_t_this ablati_this vali_t_this] = single_type_basic (abl_type, anims, force_redo, 't');
    [irr n_w_cells_this irr2 score_w_this irr3 vali_w_this] = single_type_basic (abl_type, anims, force_redo, 'w');
    n_cells = [n_cells n_cells_this];
    n_t_cells = [n_t_cells n_t_cells_this'];
    n_w_cells = [n_w_cells n_w_cells_this'];
    n_ablated = [n_ablated n_ablated_this];
    score_t_all = [score_t_all score_t_this];
    score_w_all = [score_w_all score_w_this];
    vali_t_all = [vali_t_all vali_t_this];
    vali_w_all = [vali_w_all vali_w_this];
    ablati_all = [ablati_all ablati_this];
    type_ablated = [type_ablated repmat('t',1,length(n_cells_this))];

    abl_type = 'silent';
    anims = get_anims(abl_type);
    all_anims = [all_anims anims];
    first_anims = {'257218','258836','271211','278937','278939','275801','275798','276013'};
    is_first_this = zeros(1,length(anims)); 
    is_first_this(find(ismember(anims,first_anims))) = 1;
    is_first = [is_first is_first_this]; 
    [n_cells_this n_t_cells_this n_ablated_this score_t_this ablati_this vali_t_this] = single_type_basic (abl_type, anims, force_redo, 't');
    [irr n_w_cells_this irr2 score_w_this irr3 vali_w_this] = single_type_basic (abl_type, anims, force_redo, 'w');
    n_cells = [n_cells n_cells_this];
    n_t_cells = [n_t_cells n_t_cells_this'];
    n_w_cells = [n_w_cells n_w_cells_this'];
    n_ablated = [n_ablated n_ablated_this];
    score_t_all = [score_t_all score_t_this];
    score_w_all = [score_w_all score_w_this];
    vali_t_all = [vali_t_all vali_t_this];
    vali_w_all = [vali_w_all vali_w_this];
    ablati_all = [ablati_all ablati_this];
    type_ablated = [type_ablated repmat('s',1,length(n_cells_this))];

    abl_type = 'whisking';
    anims = get_anims(abl_type);
    all_anims = [all_anims anims];
    first_anims = {'281915','278288','278759'};
    is_first_this = zeros(1,length(anims)); 
    is_first_this(find(ismember(anims,first_anims))) = 1;
    is_first = [is_first is_first_this]; 

    [n_cells_this n_t_cells_this n_ablated_this score_t_this ablati_this vali_t_this] = single_type_basic (abl_type, anims, force_redo, 't');
    [irr n_w_cells_this irr2 score_w_this irr3 vali_w_this] = single_type_basic (abl_type, anims, force_redo, 'w');
    n_cells = [n_cells n_cells_this];
    n_t_cells = [n_t_cells n_t_cells_this'];
    n_w_cells = [n_w_cells n_w_cells_this'];
    n_ablated = [n_ablated n_ablated_this];
    score_t_all = [score_t_all score_t_this];
    score_w_all = [score_w_all score_w_this];
    vali_t_all = [vali_t_all vali_t_this];
    vali_w_all = [vali_w_all vali_w_this];
    ablati_all = [ablati_all ablati_this];
    type_ablated = [type_ablated repmat('w',1,length(n_cells_this))];

    disp('=================================================================================');
    disp('NOTE: DOES NOT EXCLUDE CELLS THAT ARE OF BOTH TYPE, WHICH WE DO IN FINAL ANALYSIS');
    disp('=================================================================================');
    disp('Animal  Cell count   Ablation type    Ablation count   Touch count   Wh count');
    for a=1:length(all_anims)
        disp(sprintf('%s  %08d     %c                %08d         %08d      %08d', all_anims{a}, n_cells(a), type_ablated(a), n_ablated(a), n_t_cells(1,a), n_w_cells(1,a)));
    end
    
    disp('=================================================================================');
    disp(sprintf('SUM     %08d     %c                %08d         %08d      %08d', sum(n_cells), 'x', sum(n_ablated), sum(n_t_cells(1,:)), sum(n_w_cells(1,:))));
    disp('=================================================================================');
    firsti = find(is_first);

    disp(sprintf('Cell count: %d +/- %d ; n=%d', round(nanmean(n_cells(firsti))), round(nanstd(n_cells(firsti))), length(firsti)));
    disp(sprintf('T cell count: %d +/- %d range: %d to %d frac: %0.3f +/- %0.3f', round(nanmean(n_t_cells(1,firsti))), round(nanstd(n_t_cells(1,firsti))), ...
                 min(n_t_cells(1,firsti)), max(n_t_cells(1,firsti)), nanmean(n_t_cells(1,firsti)./n_cells(firsti)), nanstd(n_t_cells(1,firsti)./n_cells(firsti))));
    disp(sprintf('W cell count: %d +/- %d range: %d to %d frac: %0.3f +/- %0.3f', round(nanmean(n_w_cells(1,firsti))), round(nanstd(n_w_cells(1,firsti))), ...
                 min(n_w_cells(1,firsti)), max(n_w_cells(1,firsti)), nanmean(n_w_cells(1,firsti)./n_cells(firsti)), nanstd(n_w_cells(1,firsti)./n_cells(firsti))));

    disp('=================================================================================');
    touchi = find(type_ablated == 't');
    disp(sprintf('Ablated T cell count: %0.3f +/- %0.3f %d mice ; frac: %0.3f ; frac barrel: %0.3f', nanmean(n_ablated(touchi)), nanstd(n_ablated(touchi)), length(touchi), nanmean(n_ablated(touchi)./n_t_cells(1,touchi)),  nanmean(n_ablated(touchi)./287)));
    q_cell = [];
    for aa=1:length(touchi)
        a = touchi(aa);
        score_t = score_t_all{a};
        score_w = score_w_all{a};
        med_touch(a,:) = nanmedian(score_t(vali_t_all{a},:));
        med_whisk(a,:) = nanmedian(score_w(vali_w_all{a},:));

        score_sorted = sort(score_t(vali_t_all{a},1), 'ascend');
        for ci=1:length(ablati_all{a})
            ii = min(find(score_sorted > score_t(ablati_all{a}(ci))));
            if (length(ii) == 0) ; ii = 1; end
            q_cell(end+1) = ii/length(score_sorted);
        end
    end
    disp(sprintf('  ablated percentinle: %d +/- %d', round(100*nanmean(q_cell)), round(100*std(q_cell))));
    disp(sprintf('  T score pre grand median +/- corrected MAD: %0.3f +/- %0.3f post: %0.3f +/- %0.3f', nanmedian(med_touch(:,1)), 1.4826*mad(med_touch(:,1)), nanmedian(med_touch(:,2)), 1.4826*mad(med_touch(:,2))));
    disp(sprintf('  n T cells pre: +/- SD: %0.3f +/- %0.3f post: %0.3f +/- %0.3f', nanmean(n_t_cells(1,touchi)), nanstd(n_t_cells(1,touchi)),  nanmean(n_t_cells(2,touchi)), nanstd(n_t_cells(2,touchi)))); 
    disp(sprintf('  W score pre grand median +/- corrected MAD: %0.3f +/- %0.3f post: %0.3f +/- %0.3f', nanmedian(med_whisk(:,1)), 1.4826*mad(med_whisk(:,1)), nanmedian(med_whisk(:,2)), 1.4826*mad(med_whisk(:,2))));
    disp(sprintf('  n W cells pre: +/- SD: %0.3f +/- %0.3f post: %0.3f +/- %0.3f', nanmean(n_w_cells(1,touchi)), nanstd(n_w_cells(1,touchi)),  nanmean(n_w_cells(2,touchi)), nanstd(n_w_cells(2,touchi)))); 
    disp(sprintf('  # touch cells pre, total: %d whisking: %d', sum(n_t_cells(1,touchi)), sum(n_w_cells(1,touchi))));

    disp('=================================================================================');
    whiski = find(type_ablated == 'w');
    disp(sprintf('Ablated W cell count: %0.3f +/- %0.3f %d mice ; frac: %0.3f ; frac barrel: %0.3f', nanmean(n_ablated(whiski)), nanstd(n_ablated(whiski)), length(whiski), nanmean(n_ablated(whiski)./n_t_cells(1,whiski)),  nanmean(n_ablated(whiski)./287)));
    
    q_cell = [];
    for aa=1:length(whiski)
        a = whiski(aa);
        score_t = score_t_all{a};
        score_w = score_w_all{a};
        med_touch(a,:) = nanmedian(score_t(vali_t_all{a},:));
        med_whisk(a,:) = nanmedian(score_w(vali_w_all{a},:));

        score_sorted = sort(score_t(vali_t_all{a},1), 'ascend');
        for ci=1:length(ablati_all{a})
            ii = min(find(score_sorted > score_w(ablati_all{a}(ci))));
            if (length(ii) == 0) ; ii = 1; end
            q_cell(end+1) = ii/length(score_sorted);
        end
    end
    disp(sprintf('  ablated percentinle: %d +/- %d', round(100*nanmean(q_cell)), round(100*std(q_cell))));
    disp(sprintf('  T score pre grand median +/- MAD: %0.3f +/- %0.3f post: %0.3f +/- %0.3f', nanmedian(med_touch(:,1)), 1.4826*mad(med_touch(:,1)), nanmedian(med_touch(:,2)), 1.4826*mad(med_touch(:,2))));
    disp(sprintf('  n T cells pre: +/- SD: %0.3f +/- %0.3f post: %0.3f +/- %0.3f', nanmean(n_t_cells(1,whiski)), nanstd(n_t_cells(1,whiski)),  nanmean(n_t_cells(2,whiski)), nanstd(n_t_cells(2,whiski)))); 
    disp(sprintf('  W score pre grand median +/- corrected MAD: %0.3f +/- %0.3f post: %0.3f +/- %0.3f', nanmedian(med_whisk(:,1)), 1.4826*mad(med_whisk(:,1)), nanmedian(med_whisk(:,2)), 1.4826*mad(med_whisk(:,2))));
    disp(sprintf('  n W cells pre: +/- SD: %0.3f +/- %0.3f post: %0.3f +/- %0.3f', nanmean(n_w_cells(1,whiski)), nanstd(n_w_cells(1,whiski)),  nanmean(n_w_cells(2,whiski)), nanstd(n_w_cells(2,whiski)))); 
    disp(sprintf('  # touch cells pre, total: %d whisking: %d', sum(n_t_cells(1,whiski)), sum(n_w_cells(1,whiski))));


    disp('=================================================================================');
    silenti = find(type_ablated == 's');
    disp(sprintf('Ablated Si cell count: %0.3f +/- %0.3f %d mice ; frac: %0.3f', nanmean(n_ablated(silenti)), nanstd(n_ablated(silenti)), length(silenti), nanmean(n_ablated(silenti)./n_t_cells(1,silenti))));
    q_cell = [];
    for aa=1:length(silenti)
        a = silenti(aa);
        score_t = score_t_all{a};
        score_w = score_w_all{a};
        med_touch(a,:) = nanmedian(score_t(vali_t_all{a},:));
        med_whisk(a,:) = nanmedian(score_w(vali_w_all{a},:));
    end
    disp(sprintf('  T score pre grand median +/- MAD: %0.3f +/- %0.3f post: %0.3f +/- %0.3f', nanmedian(med_touch(:,1)), 1.4826*mad(med_touch(:,1)), nanmedian(med_touch(:,2)), 1.4826*mad(med_touch(:,2))));
    disp(sprintf('  n T cells pre: +/- SD: %0.3f +/- %0.3f post: %0.3f +/- %0.3f', nanmean(n_t_cells(1,silenti)), nanstd(n_t_cells(1,silenti)),  nanmean(n_t_cells(2,silenti)), nanstd(n_t_cells(2,silenti)))); 
    disp(sprintf('  W score pre grand median +/- corrected MAD: %0.3f +/- %0.3f post: %0.3f +/- %0.3f', nanmedian(med_whisk(:,1)), 1.4826*mad(med_whisk(:,1)), nanmedian(med_whisk(:,2)), 1.4826*mad(med_whisk(:,2))));
    disp(sprintf('  n W cells pre: +/- SD: %0.3f +/- %0.3f post: %0.3f +/- %0.3f', nanmean(n_w_cells(1,silenti)), nanstd(n_w_cells(1,silenti)),  nanmean(n_w_cells(2,silenti)), nanstd(n_w_cells(2,silenti)))); 
    disp(sprintf('  # touch cells pre, total: %d silenting: %d', sum(n_t_cells(1,silenti)), sum(n_w_cells(1,silenti))));


function [n_cells n_class_cells n_ablated score_all ablati_all vali_all] = single_type_basic (abl_type, anims, force_redo, resp_type)
    anim_labels = 0;
    thresh = 0.1;
    thresh_meth = 1;

    score_pre_overall = [];
    score_post_overall = [];
    score_ablated = [];

    disp(' ');
    disp('=================================================================================');
    disp(['  Ablation type: ' abl_type ' Response category: ' resp_type]);

    for a=1:length(anims)
        anim = anims{a};
        
        [data data_columns] = get_ablation_data (anim, abl_type, force_redo, 1);
        [score_pre score_post vali ablati] = get_restricted_neurons (data, data_columns, resp_type, thresh, thresh_meth);

        score_all{a}= [score_pre score_post];
        ablati_all{a} = ablati;
        vali_all{a} = vali;

        % pre vs post
        score_pre_overall(a) =  nanmedian(score_pre(vali));
        score_post_overall(a) =  nanmedian(score_post(vali));
        score_ablated(a) = nansum(score_pre(ablati)); 

        disp(sprintf('%s Total cells: %d number meeting criteria: %d (frac: %0.3f) number ablated: %d pre score: %0.3f post score: %0.3f delta: %0.3f', ...
                      anim, size(data,1), length(vali), length(vali)/size(data,1), length(ablati), score_pre_overall(a), score_post_overall(a),  score_post_overall(a)-score_pre_overall(a)));

        % very careful -- this, unlike everything else, does NOT exclude intersections; use this to report touch cell count, but as noted in manuscript, analyses exclude intersection
        n_cells(a) = size(data,1);
        n_ablated(a) = length(ablati);
        n_class_cells(a,1) = length(find(score_pre > thresh));
        n_class_cells(a,2) = length(find(score_post > thresh));
    end

