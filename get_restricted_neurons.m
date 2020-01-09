function [score_pre score_post vali ablati] = get_restricted_neurons (data, data_columns, resp_type, thresh, thresh_meth)
%
% tells you who to use for a given 
%
    [t_pre t_post vali_t ablati] = get_restricted_neurons_one_type (data, data_columns, 't', thresh, thresh_meth);
    [w_pre w_post vali_w ablati] = get_restricted_neurons_one_type (data, data_columns, 'w', thresh, thresh_meth);

    % remove intersection
    vali_both = intersect(vali_t, vali_w);
    vali_t = setdiff(vali_t, vali_both);
    vali_w = setdiff(vali_w, vali_both);

    switch resp_type
        case 't'
            score_pre = t_pre;
            score_post = t_post;
            vali = vali_t;

        case 'w'
            score_pre = w_pre;
            score_post = w_post;
            vali = vali_w;
    end

function  [score_pre score_post vali ablati] = get_restricted_neurons_one_type (data, data_columns, resp_type, thresh, thresh_meth)

    % apply criteria: 1) thresh pre or post, and not pre/post _bad 2) no proximity exclude 3) not ablated
    score_pre = data(:,find(strcmp(data_columns, [resp_type '_pre'])));
    score_post = data(:,find(strcmp(data_columns, [resp_type '_post'])));

    % test dataset -- if > 25% points are bad, report
    badi = union(find(isnan(data(:,find(strcmp(data_columns, [resp_type '_pre']))))),find(isnan(data(:,find(strcmp(data_columns, [resp_type '_post']))))));
    if (length(badi) > 0.25*size(data,1))
        disp(sprintf('  WARNING: LOTS OF BAD NEURONS: %d of %d ; %0.3f', length(badi), size(data,1), length(badi)/size(data,1)));
    end

    % simple threshold - bulky but easy to explain / understand
    if (thresh_meth == 1)
        vali_pre = find(data(:,find(strcmp(data_columns, [resp_type '_pre']))) > thresh & ~(data(:,find(strcmp(data_columns, [resp_type '_pre_bad'])))));
        vali_post = find(data(:,find(strcmp(data_columns, [resp_type '_post']))) > thresh & ~(data(:,find(strcmp(data_columns, [resp_type '_post_bad'])))));
        vali = union(vali_pre,vali_post);

        % remove missing values (nan)
        vali = setdiff(vali,badi);
    else % sum is more elegant, but this is only used in an analysis not in the manuscript (shift in score distribution) as it does not have gap in distribution
        vali = find(data(:,find(strcmp(data_columns, [resp_type '_pre']))) + data(:,find(strcmp(data_columns, [resp_type '_post'])))  > thresh );
    end
    vali = setdiff(vali,find(data(:,find(strcmp(data_columns, 'proximity_exclude')))));
    ablati = find(data(:,find(strcmp(data_columns, 'is_ablated'))));
    vali = setdiff(vali,ablati);


