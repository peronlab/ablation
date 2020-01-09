function badi = get_rejects_based_on_shuffle (real_vec, sh_vec, er, score_thresh, subi)
% gets rejects baed on event-rate matched bins of the shuffled data
%   real_vec: actual score vectore
%   sh_vec: shuffled vector for all cells
%   er: event rate
%   score_thresh: threshold to use
%   subi: idx of cells to consider (typically subvolume)
    er_cutoffs = quantile(er(subi), 0:.1:1);
    er_cutoffs(end) = er_cutoffs(end)+.0001; % so below logic works
    badi = [];
    for e=1:length(er_cutoffs)-1
        er_bini = subi(find(er(subi) >= er_cutoffs(e) & er(subi) < er_cutoffs(e+1)));

        thresh_sh = quantile(sh_vec(er_bini),score_thresh); 
        badi = [badi er_bini(find(real_vec(er_bini) < thresh_sh))];
    end
    
