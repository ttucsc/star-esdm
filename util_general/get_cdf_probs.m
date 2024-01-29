function [out_cdf, out_okflags, anchor_probs, av2_anchors, av2_anchor_ix] = get_cdf_probs(anchor_pts, do_exponential, ref_cdf, ref_bins, ref_okflags, ref_mean, av2_cdf, av2_bins, av2_okflags, av2_mean)
% Appends an empirical distribution onto the calculated probability surface beyond the point where the local
% calculations are likely correct.  This as about +,-3 sigma.  Empirical distributions for each model,
% (passed in as ref_cdf), come from aggregating data for ~1000 GHCN sites across N america (sites w/ > 55 yrs data)
%   Empirical distributions are generated from the MODEL data, not the GHCN observations;  but the list of locations was
%   generated to match a list of 1000 stations with lots of data.
%   Emprical distributions are made for each model, but combine data across the entire year, using 55 years' data from
%   1950-2005, merging data from all 1000 sites.  Merging is done by normalizing to the 2.5-sigma-equivalent point
%   (i.e., the point where the CDF's probability is ~2.3e-4, or 2.5 sigma on a normal distribution.  This removes
%   dependence on time of year.

%   Code interpolates the reference distribution to the ARRM_V2 bins, based on the ratio to the inner anchor points 
%   (usually +,- ~3 sigmas).  
%   Beyond the second (outer)anchor points (usually 4.4 or 4.5 sigmas, about 5 or3e-6), it assumes a straight 
%   exponential decay, down to 1e-17.  Flags everything beyond the outer anchor points as mapping to NAs, but other
%   parts of the code need a reasonably smooth probability out to 1e-15 or so to avoid discontinuities.
%
%   Inputs:
%       ref_cdf         model's overall reference cdf
%       ref_bins        bins where ref_cdf is evaluated
%       ref_okflags     flags of where to use the cdf (increasing points...i.e., all places where pdf is non-zero))
%       anchor_pts      places to resample relative to.  Usually -3-sigma-equivalent and +3-sigma-equivalent points
%       av2_cdf         ARRM V2's cdf for site
%       av2_bins        ARRM V2's bins for av2_cdf
%       av2_okflags     flags of where to use the cdf

    ref_bins_0 = ref_bins - ref_mean;        % adjust bins relative to mean.
    av2_bins_0 = av2_bins - av2_mean;
    av2_ctr_ix = find(av2_cdf >= 0.5,1);
    ref_ctr_ix = find(ref_cdf >= 0.5,1);
    anchor_probs = normcdf(anchor_pts);
    av2_nbins = length(av2_bins);
    ref_anchors = interp1(ref_cdf(ref_okflags), ref_bins_0(ref_okflags), anchor_probs(1,:), 'makima');
    av2_anchors = interp1(av2_cdf(av2_okflags), av2_bins_0(av2_okflags), anchor_probs(1,:), 'makima');
%   ref_anchor_ix = [find(ref_bins <= ref_anchors(1),1,'last'), find(ref_bins >= ref_anchors(2),1)]; 
    av2_anchor_ix = [find(av2_bins <= av2_anchors(1),1,'last'), find(av2_bins >= av2_anchors(2),1)]; 
        
    ref_ratios_neg = ref_bins_0(1:(ref_ctr_ix-1)) / ref_anchors(1);         % left-tail  ratios of each reference bin to the anchor pt (e.g., 3-sigma-equivalent point)
    ref_ratios_pos = ref_bins_0(ref_ctr_ix:end)   / ref_anchors(2);         % right-tail ratios of each reference bin to the anchor pt (e.g., 3-sigma-equivalent point)
    av2_ratios_neg = av2_bins_0(1:(av2_ctr_ix-1)) / av2_anchors(1);         % left-tail  ratios of each ARRM_V2 bin to the anchor pt (e.g., 3-sigma-equivalent point)
    av2_ratios_pos = av2_bins_0(av2_ctr_ix:end)   / av2_anchors(2);         % right-tail ratios of each ARRM_V2 bin to the anchor pt (e.g., 3-sigma-equivalent point)
%     ref_ratios_neg = ref_bins(1:(ref_ctr_ix-1)) / ref_bins(ref_anchor_ix(1));         % left-tail  ratios of each reference bin to the anchor pt (e.g., 3-sigma-equivalent point)
%     ref_ratios_pos = ref_bins(ref_ctr_ix:end)   / ref_bins(ref_anchor_ix(2));         % right-tail ratios of each reference bin to the anchor pt (e.g., 3-sigma-equivalent point)
%     av2_ratios_neg = av2_bins(1:(av2_ctr_ix-1)) / av2_bins(av2_anchor_ix(1));         % left-tail  ratios of each ARRM_V2 bin to the anchor pt (e.g., 3-sigma-equivalent point)
%     av2_ratios_pos = av2_bins(av2_ctr_ix:end)   / av2_bins(av2_anchor_ix(2));         % right-tail ratios of each ARRM_V2 bin to the anchor pt (e.g., 3-sigma-equivalent point)
    out_cdf = nan(size(av2_cdf));
    out_cdf(1:(av2_ctr_ix-1)) = interp1(ref_ratios_neg, ref_cdf(1:(ref_ctr_ix-1)), av2_ratios_neg, 'makima',0);
    out_cdf(av2_ctr_ix:end)   = interp1(ref_ratios_pos, ref_cdf(ref_ctr_ix:end),   av2_ratios_pos, 'makima',1);
    
    if (~do_exponential)
        out_okflags = ~isnan(out_cdf);
        return; 
    end
    
        % find where to start exponential decay (anchor_probs(2,:)
    min_ix = max(10,ceil(.25*av2_anchor_ix(1)));    % Kludgey safety check. start exponential decay no more than 1/4  of the way in from the edge to the first anchor point.
    ix0_neg = max(min_ix, find(out_cdf >= anchor_probs(2,1),1));

    prob0_neg = out_cdf(ix0_neg);
    x0_neg = av2_bins(ix0_neg);
    
    max_ix = min(av2_nbins-9,floor(av2_anchor_ix(2) + .75*(av2_nbins-av2_anchor_ix(2))));
    ix0_pos = min(max_ix,find(out_cdf <= anchor_probs(2,2),1,'last')+1);

    prob0_pos = 1-out_cdf(ix0_pos);     % note:  right tail, but we want exponential decay down to zero,, then subtract that from 1, so we subtract 1 here.
    x0_pos = av2_bins(ix0_pos);

%     out_cdf1 = out_cdf;
            % now calculate exponential decay from end point of empirical out to end of av2 bins.
            
%     k_neg = (log(prob0_neg) - log(3e-16))/(av2_bins(1)  -x0_neg);       % below 3e-16, 1-prob rounds to 1...
%     k_pos = (log(prob0_pos) - log(3e-16))/(x0_pos-av2_bins(end));
% 
    k_neg = (log(prob0_neg) - log(1e-17))/(av2_bins(1)  -x0_neg);       % below 3e-16, 1-prob rounds to 1...
    k_pos = (log(prob0_pos) - log(1e-17))/(x0_pos-av2_bins(end));

    out_cdf(1:(ix0_neg-1))   =     prob0_neg * exp(k_neg * (x0_neg - av2_bins(1:(ix0_neg-1))));
    out_cdf((ix0_pos+1):end) = 1 - prob0_pos * exp(k_pos * (av2_bins((ix0_pos+1):end) - x0_pos));
    
%     if (do_plot)
%         figure(98);  
%         clf; 
%         semilogy(av2_bins,av2_cdf, av2_bins, out_cdf1, av2_bins, 1-av2_cdf, av2_bins, 1-out_cdf1)
%         hold on;
%         semilogy(av2_bins, out_cdf, 'g', av2_bins,1-out_cdf,'g','linewidth',2)
%         hold off;
%     end
%     
    out_okflags = true(size(out_cdf));
    out_okflags(1:(ix0_neg-1))   = false;
    out_okflags((ix0_pos+1):end) = false;
end

