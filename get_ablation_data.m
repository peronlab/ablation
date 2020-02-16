function [data data_columns] = get_ablation_data (anim, abl, force_redo, suppress_output)
% 
% Returns a matrix with the requested data:
%
%   anim: string representation of animal -- EXCLUDES an
%   abl: which ablation ('touch', 'silent', 'whisking')
%   force_redo: regenerate data file (usually it just loads it and returns it)
%   suppress_output: silence text
%
% Columns in matrix specified in data_columns; data is matrix itself
%
%   abl_type: 0 = silent ; 1 = touch ; 2 = whisking
%   
    settings = get_settings;
    
    if (nargin < 3) ; force_redo = 0 ;end
    if (nargin < 4) ; suppress_output = 1 ;end

    % output directory for storing reduced data for fast access
    if (~exist(settings.processed_data_path,'dir'))
        disp([settings.processed_data_path ' does not exist; will create, or, abort then edit get_settings.m to put another directory in there instead.']);
        pause;
        mkdir(settings.processed_data_path);
    end
    
    % settings
    score_thresh = 0.99; % must be this percentile of scores or higher to be valid member of a representation
    xy_min_dist_rad = 10; % excluded, min
    z_min_dist_vert = 30; % excluded, min

    if (strcmp(abl, 'silent'))  
        abl_type = 0;
    elseif (strcmp(abl,'touch'))
        abl_type = 1;
    elseif (strcmp(abl,'whisking'))
        abl_type = 2;
    end

    data_file = sprintf('%s%c%s_%s_ablation_data.mat', settings.processed_data_path, filesep, anim, abl);

    data_columns = {'anim_id', 'abl_type', 'cell_id', 'subvol', 't_pre', 't_post', 't_di_pre', 't_di_post', 'w_pre', 'w_post', 't_pre_bad', 't_post_bad', 'w_pre_bad', 'w_post_bad', ...
                    'er_pre', 'er_post', 'is_ablated', 'x_um', 'y_um', 'z_um', 'proximity_exclude', 'd_ablated_um'};

    % creation block
    if (~exist(data_file) | force_redo)
        disp(sprintf('%s not found (or redo forced); creating.', data_file));

        [pre_path post_path abl_curation_file_path subvol_restrict] = get_data_path(anim, abl);

        % pull the raw data
        [t_pre t_dff_pre t_di_pre w_pre w_dff_pre t_pre_sh w_pre_sh er_pre x y z ids subvol] = pull_struct(pre_path, anim);
        if (~issorted(ids))
            disp('SORT FAIL ; fixing');
            [ids sorti] = sort(ids,'ascend');
            vars = {'t_pre','t_di_pre','w_pre','t_pre_sh','w_pre_sh','er_pre','x','y','z','subvol'};
            for v=1:length(vars)
                eval([vars{v} '=' vars{v} '(sorti);']);
            end
        end

        % apply subvolume exclusion here
        goodi = find(ismember(subvol, subvol_restrict));
        if (length(goodi) ~= length(subvol) && ~isnan(subvol_restrict(1)))
            disp('Applying subvolume restriction');
            vars = {'t_pre','t_di_pre','w_pre','t_pre_sh','w_pre_sh','er_pre','x','y','z','ids','subvol'};
            for v=1:length(vars)
                eval([vars{v} '=' vars{v} '(goodi);']);
            end
        end
            
        [t_post t_dff_post t_di_post w_post w_dff_post t_post_sh w_post_sh er_post irr_x irr_y irr_z ids_post] = pull_struct(post_path, anim);

        if (length(ids_post) ~= length(ids))
            disp('Pre/post mismatch ; fixing');
            [irr i_pre i_post] = intersect(ids, ids_post);
            vars = {'ids','t_pre','t_di_pre','w_pre','t_pre_sh','w_pre_sh','er_pre','x','y','z','subvol'};
            for v=1:length(vars)
                eval([vars{v} '=' vars{v} '(i_pre);']);
            end
            vars = {'t_post','t_di_post','w_post','t_post_sh','w_post_sh','er_post', 'ids_post'};
            for v=1:length(vars)
                eval([vars{v} '=' vars{v} '(i_post);']);
            end
        end

        if (sum(ids == ids_post) ~= length(ids)) 
            disp('IDS MISMATCH - fixing'); 
            if (length(ids) ~= sum(sort(ids_post,'ascend') == ids)) ; disp('this is well and truly Fd') ; pause; end
            vars = {'t_post','t_di_post','w_post','t_post_sh','w_post_sh','er_post'};
            [irr sorti] = sort(ids_post,'ascend');
            for v=1:length(vars)
                eval([vars{v} '=' vars{v} '(sorti);']);
            end
        end

        % compute shuffle-matched quantile for each cell -- nan out the relevant guys this way
        usv = unique(subvol);
        t_pre_badi = [];
        t_post_badi = [];
        w_pre_badi = [];
        w_post_badi = [];
        for u=1:length(usv)
            subi = find(subvol == usv(u));
            if (length(subi) > 0)
                t_pre_badi = [t_pre_badi get_rejects_based_on_shuffle(t_pre, t_pre_sh, er_pre, score_thresh, subi)];
                t_post_badi = [t_post_badi get_rejects_based_on_shuffle(t_post, t_post_sh, er_post, score_thresh, subi)];
                w_pre_badi = [w_pre_badi get_rejects_based_on_shuffle(w_pre, w_pre_sh, er_pre, score_thresh, subi)];
                w_post_badi = [w_post_badi get_rejects_based_on_shuffle(w_post, w_post_sh, er_post, score_thresh, subi)];
            end
        end
        t_pre_bad = zeros(length(ids),1);
        t_post_bad = zeros(length(ids),1);
        w_pre_bad = zeros(length(ids),1);
        w_post_bad = zeros(length(ids),1);
        t_pre_bad(t_pre_badi) = 1;
        t_post_bad(t_post_badi) = 1;
        w_pre_bad(w_pre_badi) = 1;
        w_post_bad(w_post_badi) = 1;

        % ablated ids get
        abl_ids = get_ablated_ids(anim, abl);
        abli = find(ismember(ids, abl_ids));
        abli_vec = zeros(length(ids),1);
        abli_vec(abli) = 1;

        % compute positional stuff
        prox_excludi = [];
        d_ablated_um = inf*ones(length(x),1);
        proximity_exclude = zeros(length(x),1);
        for a=1:length(abli)
            d_ablated_um = min(d_ablated_um, sqrt((x-x(abli(a))).^2 + (y-y(abli(a))).^2 + (z-z(abli(a))).^2)');
                
            d_rad = sqrt((x-x(abli(a))).^2 + (y-y(abli(a))).^2);
            d_z = (z-z(abli(a)));
            prox_excludi = union(prox_excludi, find(d_rad < xy_min_dist_rad & d_z < z_min_dist_vert));
        end
        proximity_exclude(prox_excludi) = 1;

        % populate matrix
        data = nan*zeros(length(ids), length(data_columns));
        data(:,find(strcmp(data_columns, 'anim_id'))) = ones(length(t_pre),1)*str2num(anim);
        data(:,find(strcmp(data_columns, 'abl_type'))) = ones(length(t_pre),1)*abl_type;
        data(:,find(strcmp(data_columns, 'cell_id'))) = ids;
        data(:,find(strcmp(data_columns, 'subvol'))) = subvol;
        data(:,find(strcmp(data_columns, 't_pre'))) = t_pre;
        data(:,find(strcmp(data_columns, 't_post'))) = t_post;
        data(:,find(strcmp(data_columns, 't_di_pre'))) = t_di_pre;
        data(:,find(strcmp(data_columns, 't_di_post'))) = t_di_post;
        data(:,find(strcmp(data_columns, 'w_pre'))) = w_pre;
        data(:,find(strcmp(data_columns, 'w_post'))) = w_post; 
        data(:,find(strcmp(data_columns, 't_pre_bad'))) = t_pre_bad;
        data(:,find(strcmp(data_columns, 't_post_bad'))) = t_post_bad;
        data(:,find(strcmp(data_columns, 'w_pre_bad'))) = w_pre_bad;
        data(:,find(strcmp(data_columns, 'w_post_bad'))) = w_post_bad;
        data(:,find(strcmp(data_columns, 'er_pre'))) = er_pre;
        data(:,find(strcmp(data_columns, 'er_post'))) = er_post;

        data(:,find(strcmp(data_columns, 'x_um'))) = x;
        data(:,find(strcmp(data_columns, 'y_um'))) = y;
        data(:,find(strcmp(data_columns, 'z_um'))) = z;

        data(:,find(strcmp(data_columns, 'is_ablated'))) = abli_vec;
        data(:,find(strcmp(data_columns, 'proximity_exclude'))) = proximity_exclude;
        data(:,find(strcmp(data_columns, 'd_ablated_um'))) = d_ablated_um;

        save(data_file, 'data', 'data_columns', 'abl', 'anim');

        disp(sprintf('Nan entries for tpre: %d tpost: %d wpre: %d wpost: %d', length(find(isnan(t_pre))), length(find(isnan(t_post))), length(find(isnan(w_pre))), length(find(isnan(w_post)))));
    else
        if(~suppress_output) ; disp(['Loading file ' data_file]); end
        load(data_file);
    end

% 
% pulls key data from session object(s) ; x y z should be in MICRONS ; this method is in use
%
function [t_score t_dff t_di w_score w_dff t_score_sh w_score_sh e_rate x y z ids subvol] = pull_struct(dat_path, anim)
    pix2um = 600/512;
    
    % GLM feature names
    whisk_feat = 'eventBasedDff_k16t14_MeanWhiskerAmplitudeGLMScore';
    touch_feat_base = 'eventBasedDff_k16t14_%sAbsMaxKappaZeroNotouchGLMScore';

    t_score = [];
    t_dff = [];
    w_score = [];
    w_dff = [];
    t_score_sh = [];
    w_score_sh = [];
    e_rate = [];
    x = [];
    y = [];
    z = [];
    ids = [];
    subvol = [];
    t_di = [];

    cd (dat_path);
    fl = dir('an*sess_struct.mat'); 
    for f=1:length(fl) ; 
        load(fl(f).name) ; 
        s = sess;

        if (f == 1) % get some stuff
            wtag = s.whiskerTag{1};

            touch_feat = sprintf(touch_feat_base, lower(wtag));
        end

        whisk_feat_sh = strrep(whisk_feat, 'GLMScore','GLMShuffleScore');
        touch_feat_sh = strrep(touch_feat, 'GLMScore','GLMShuffleScore');
        touch_feat_kern = strrep(touch_feat, 'GLMScore','GLMStimKernelWeight');

        t_score = [t_score s.caTSA.cellFeatures.values{find(strcmp(s.caTSA.cellFeatures.keys, touch_feat))}]; 
        t_score_sh = [t_score_sh s.caTSA.cellFeatures.values{find(strcmp(s.caTSA.cellFeatures.keys, touch_feat_sh))}];
        w_score = [w_score s.caTSA.cellFeatures.values{find(strcmp(s.caTSA.cellFeatures.keys, whisk_feat))}] ; 
        w_score_sh = [w_score_sh s.caTSA.cellFeatures.values{find(strcmp(s.caTSA.cellFeatures.keys, whisk_feat_sh))}];
        
        % directionality idx
        t_kern = s.caTSA.cellFeatures.values{find(strcmp(s.caTSA.cellFeatures.keys, touch_feat_kern))};
        if (0) % use means
            nkp = floor(size(t_kern,1)/2);
            t_pro_resp = nanmean(t_kern(1:nkp,:));
            t_ret_resp = nanmean(t_kern((end-nkp+1):end,:));
            t_di_mean = (t_pro_resp-t_ret_resp)./(t_pro_resp+t_ret_resp);
            t_di = [t_di t_di_mean];
        else % use ends
            t_di_ends = (t_kern(1,:)-t_kern(end,:))./(t_kern(1,:)+t_kern(end,:));
            t_di = [t_di t_di_ends];
        end

        e_rate = [e_rate s.caTSA.cellFeatures.values{find(strcmp(s.caTSA.cellFeatures.keys, 'eventRateHz'))} ];

        image_bounds = [512 512];
        for v=1:length(s.caTSA.roiArray)
            coms = nan*zeros(length(s.caTSA.roiArray{v}.rois), 2);
            for r=1:length(s.caTSA.roiArray{v}.rois)
            	x_ = ceil(s.caTSA.roiArray{v}.rois{r}.indices/image_bounds(1));
                y_ = s.caTSA.roiArray{v}.rois{r}.indices - (floor(s.caTSA.roiArray{v}.rois{r}.indices/image_bounds(1))*image_bounds(1));

                if (length(x_) > 0)
                    coms(r,:) = [mean(x_) mean(y_)]*pix2um;
                end
            end

            z = [z ones(1,size(coms,1))*((f-1)*45 + (v-1)*15);];
            x = [x coms(:,1)'];
            y = [y coms(:,2)'];
            ids = [ids s.caTSA.roiArray{v}.roiIds];
            subvol = [subvol f+0*s.caTSA.roiArray{v}.roiIds];
        end
    end

% 
% pulls key data from session object(s) ; x y z should be in MICRONS ; this is an OLD method using a different data format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOT USED - LEFT IN FOR POSTERITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t_score t_dff t_di w_score w_dff t_score_sh w_score_sh e_rate x y z ids subvol] = pull_session_objects(dat_path, anim)
    pix2um = 600/512;
    
    % GLM feature names
    whisk_feat = 'eventBasedDff_k16t14_MeanWhiskerAmplitudeGLMScore';
    touch_feat_base = 'eventBasedDff_k16t14_%sAbsMaxKappaZeroNotouchGLMScore';

    t_score = [];
    t_dff = [];
    w_score = [];
    w_dff = [];
    t_score_sh = [];
    w_score_sh = [];
    e_rate = [];
    x = [];
    y = [];
    z = [];
    ids = [];
    subvol = [];
    t_di = [];

    cd (dat_path);
    fl = dir('an*sess.mat'); 
    for f=1:length(fl) ; 
        load(fl(f).name) ; 

        if (f == 1) % get some stuff
            wtag = s.whiskerTag{1};

            touch_feat = sprintf(touch_feat_base, lower(wtag));
        end

        whisk_feat_sh = strrep(whisk_feat, 'GLMScore','GLMShuffleScore');
        touch_feat_sh = strrep(touch_feat, 'GLMScore','GLMShuffleScore');
        touch_feat_kern = strrep(touch_feat, 'GLMScore','GLMStimKernelWeight');

        t_score = [t_score s.caTSA.cellFeatures.get(touch_feat)]; 
        t_score_sh = [t_score_sh s.caTSA.cellFeatures.get(touch_feat_sh)];
        w_score = [w_score s.caTSA.cellFeatures.get(whisk_feat)] ; 
        w_score_sh = [w_score_sh s.caTSA.cellFeatures.get(whisk_feat_sh)];
        
        % directionality idx
        t_kern = s.caTSA.cellFeatures.get(touch_feat_kern);
        if (0) % use means
            nkp = floor(size(t_kern,1)/2);
            t_pro_resp = nanmean(t_kern(1:nkp,:));
            t_ret_resp = nanmean(t_kern((end-nkp+1):end,:));
            t_di_mean = (t_pro_resp-t_ret_resp)./(t_pro_resp+t_ret_resp);
            t_di = [t_di t_di_mean];
        else % use ends
            t_di_ends = (t_kern(1,:)-t_kern(end,:))./(t_kern(1,:)+t_kern(end,:));
            t_di = [t_di t_di_ends];
        end

        e_rate = [e_rate s.caTSA.cellFeatures.get('eventRateHz')];

        for v=1:length(s.caTSA.roiArray)
            coms = s.caTSA.roiArray{v}.getCOMs*pix2um;
            z = [z ones(1,size(coms,1))*((f-1)*45 + (v-1)*15);];
            x = [x coms(:,1)'];
            y = [y coms(:,2)'];
            ids = [ids s.caTSA.roiArray{v}.roiIds];
            subvol = [subvol f+0*s.caTSA.roiArray{v}.roiIds];
        end
    end

function mat = nanjoin(A, B, fixed_dim)
    sA = size(A);
    sB = size(B);
    
    if (fixed_dim == 2)
        if (sA(1) < sB(1))
            pad = (sA(1)+1):sB(1);
            A(pad,:) = nan;
        elseif (sB(1) < sA(1))
            pad = (sB(1)+1):sA(1);
            B(pad,:) = nan;
        end
    end

    mat = [A B];

%
% gets ablated ids
%
function abl_ids = get_ablated_ids(anim, abl)
    [pre_path post_path abl_curation_file_path] = get_data_path(anim, abl);

    if (~iscell(abl_curation_file_path))
        abl_curation_file_path = {abl_curation_file_path};
    end
    
    abl_ids = [];
    for a=1:length(abl_curation_file_path)
        abdat = load(abl_curation_file_path{a});
        if (isfield(abdat, 'ablationSuccessful'))
            abl_ids = [abl_ids unique(abdat.ablationTargetIds(find(abdat.ablationSuccessful)))];
        elseif (isfield(abdat, 'ablated_ids'))
            abl_ids = [abl_ids unique(abdat.ablated_ids)];
        else
            disp([abl_curation_file_path{a} ' not a valid ablation data file']);
        end
   end
     
