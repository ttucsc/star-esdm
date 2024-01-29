function  [excess_thresh, excess_frac, wet_frac_obs, obs_clim_orig] = calc_excess_prcp(prcp_obs, prcp_mdl, obs_prcp_min, mdl_prcp_min, yrlen, ignore_nans, fignum, lbl2) % , pdf_yrlen
%
%   Calculates threshold values to Removes excess wet days from model we have either approximately the same 
%   number of wet days in the model (historical) period.
%
%   This is intended as input to trim_excess_prcp(...), which can iterate to get the exact number removed.
%
%   Inputs:  
%       prcp_obs        obs precipitation, size ( yrlen*nyrs x 1 )
%       prcp_mdl        model precip, historical period
%       prcp_min        min precip value;  values below are treated as dry days
%       yrlen           year length.  Usually 365...occasionally 360...
%       ignore_nans     do not correct counts for nans
%       fignum          (optional)  If present, draw figures as we iterate, and output steps to console.
%
%   Outputs
%       excess_thresh   threshold precip values (yrlen x 1) at which to threshold model precip to match obs precip.
%       excess_frac     Fraction of precip to remove to match obs precip.
%                           These are daily values, smoothly varying, and can't be applied directly to any day.
%                           They are based on the wet-day climatologies.  They are needed to avoid reducing wet-day
%                           counts during times when the model produces too few wet days.  If adjusting to match counts
%                           exactly, you don't want to change the thresholds for periods when the model doesn't have
%                           enough wet days to start with.
%       wet_frac        fraction of days which should be wet days.  This is a single value applied over all.
%                           It is adjusted to account for the number of nan's in the model input (usually none!).
%
%       (could pass back mdl histogram & prob. surface to avoid recalculating...)
% NEEDS:  adjust for NaNs in obs data, and possibly in model data as well.
%

    if (~exist("fignum", "var") || isempty(fignum)), fignum = 0; end
      
    if (isa(prcp_obs, "single")), fudge_obs = 1e-7;  else, fudge_obs = 1e-16; end  % fudge used to avoid rounding errors.
    if (isa(prcp_mdl, "single")), fudge_mdl = 1e-7;  else, fudge_mdl = 1e-16; end  % fudge used to avoid rounding errors.
    prcp_obs( prcp_obs < obs_prcp_min-fudge_obs) = 0;      % trim off trace obs amounts.  Use fudge to avoid computer arithmetic errors
    prcp_mdl( prcp_mdl < mdl_prcp_min-fudge_mdl) = 0;
            
    obs_nnans = sum(isnan(prcp_obs));       % # of missing points in obs, mdl
    mdl_nnans = sum(isnan(prcp_mdl));
    
    npts_obs = length(prcp_obs);
    npts_mdl = length(prcp_mdl);
    
    if (ignore_nans)
        nan_fact_obs = 1;
        nan_fact_mdl = 1;
    else
        nan_fact_obs = 1-obs_nnans/npts_obs;    % ratio of available obs & mdl data to length of data.  We'll scale stuff by this later
        nan_fact_mdl = 1-mdl_nnans/npts_mdl;
    end

    wet_frac_obs = sum(prcp_obs>0)/npts_obs/nan_fact_obs;       % overall fraction of wet-days in obs data.  This will set the target for thresholding.
    wet_frac_mdl = sum(prcp_mdl>0)/npts_mdl/nan_fact_mdl;
    mdl_excess_per_yr = max(0, wet_frac_mdl-wet_frac_obs) * yrlen;     % # of mdl precips/year which need to be removed.
    
    
        % reshape to 2-D matrices so we can get climatologies and histograms.
        
    prcp_obs = reshape(prcp_obs, yrlen, []);
    prcp_mdl = reshape(prcp_mdl, yrlen, []);
    
    nyrs_obs = size(prcp_obs, 2);
    nyrs_mdl = size(prcp_mdl, 2);
    
    mdl_prcps_per_year = sum(prcp_mdl(:)>0)/nyrs_mdl/nan_fact_mdl;
    obs_prcps_per_year = sum(prcp_obs(:)>0)/nyrs_obs/nan_fact_obs;
    diff_count_per_year = mdl_prcps_per_year - obs_prcps_per_year;  % if >= 0, then if we have more model wet days per year than obs wet days

            % some later code needs the obs count climatology, so we go ahead and calculate it here if we had to bail out because of insufficient model precip.
    obs_counts = sum(prcp_obs>=(obs_prcp_min-fudge_obs),2);
    obs_clim_orig = climatology(obs_counts/nyrs_obs, 8,3, yrlen)/nan_fact_obs;        
    
    if (diff_count_per_year <=  0)      % fewer model precips than obs precips, can't manufacture wet days, so just bail out.        
        excess_thresh = zeros(yrlen,1);
        excess_frac = zeros(yrlen,1);
        return;
    end
   
    [edges_obs, bins_obs] = calc_prcp_raw_edges(prcp_obs(:), obs_prcp_min);
    
    [edges_mdl, bins_mdl] = calc_prcp_raw_edges(prcp_mdl(:), mdl_prcp_min);
    
    
    nbins_obs = length(bins_obs);
    nbins_mdl = length(bins_mdl);
        
        % get histograms, and histogram the precip values. 
    obs_hist = zeros(yrlen, nbins_obs);
    mdl_hist = zeros(yrlen, nbins_mdl);
    
    for i=1:yrlen
        obs_hist(i,:) = histcounts(prcp_obs(i,:), edges_obs);
        mdl_hist(i,:) = histcounts(prcp_mdl(i,:), edges_mdl);
    end
    
        % zero out the counts for dry-days.
    obs_hist(:,1)=0;
    mdl_hist(:,1)=0;
    
        % total precips/year on each day
    obs_counts = sum(obs_hist,2)/nyrs_obs;      % precip events per year.
    mdl_counts = sum(mdl_hist, 2)/nyrs_mdl;
    
%   obs_clim_orig = climatology(obs_counts, 8,3, yrlen)/nan_fact_obs;  % this is done earlier now.    
%     [obs_eq_short, dayPos]  = remap_time_dimension(obs_hist', pdf_yrlen, obs_clim_orig);  % Note:  remap_time_dimension expects histogram as transpose of how we have it here.
%     [mdl_eq_short, ~     ]  = remap_time_dimension(mdl_hist', pdf_yrlen, obs_clim_orig);
%     
%     obs_eq = remap_time_dimension(obs_eq_short, yrlen, [], dayPos);
%     mdl_eq = remap_time_dimension(mdl_eq_short, yrlen, [], dayPos);
%     
%         % Scale the equalized counts back by the width of the day in the
%         % histogram-equalizing.
%     dayWidth = diff([0,dayPos]);
%     for i=1:yrlen
%         obs_eq(:,i) = obs_eq(:,i) .* dayWidth(i)/nyrs_obs;      % daily obs counts/year after smoothing
%         mdl_eq(:,i) = mdl_eq(:,i) .* dayWidth(i)/nyrs_mdl;      % daily mdl counts/year after smoothing.
%     end
% 
%         % total precips/year on each day
%     obs_eq_counts = sum(obs_eq, 1);      % precip events per year.
%     mdl_eq_counts = sum(mdl_eq, 1);
%     
    
            %   these should be scaled by daily nan counts, rather than by single number for year.
            %   probably best to fold nans into compressed year first.  Dealing with days w/ all or large # of NAs could
            %   get rather tricky, so ignoring the problem for now.  Ignoring == assuming nans are randomly located
            %   throughout the year.
    obs_clim   = climatology(obs_counts, 8,3, yrlen)/nan_fact_obs;  % /nan_fact is scaling to account for missing data.
    mdl_clim   = climatology(mdl_counts, 8,3, yrlen)/nan_fact_mdl;
%   clim_mdl_excess = sum(mdl_clim) - sum(obs_clim);  % how many to remove per year
    
       % Just in case the climatology goes negative...bump it up just above 0.
       %    shouldn't be needed now, Ian, with histogram equalization.
       %
       %        This probably should be higher, maybe 1%?  
    obs_clim = max(obs_clim, .0001);  
    mdl_clim = max(mdl_clim, .0001);
    
    if (fignum>0)
        figure(fignum);
        subplot(2,1,1);
        plot(1:yrlen, obs_clim,   'b-', 1:yrlen, mdl_clim, 'r-', 1:yrlen,  mdl_clim - obs_clim, "g-", "linewidth",2);
        hold on;
        plot(1:yrlen, obs_counts, 'b-', 1:yrlen, mdl_counts, 'r-');
        hold off;
        grid on;
        title("wet day climatology");
    %   legend("obs","mdl","obs-mdl", "mdl_excess");
       legend("obs","mdl","obs-mdl");
       xlabel("day of yr");
       ylabel("prob of wet day");
%         subplot(3,1,2);
%         plot(1:yrlen, obs_clim_eq, 'b', 1:yrlen, mdl_clim_eq, 'r', 1:yrlen, mdl_clim_eq -  obs_clim_eq, "g-", "linewidth",2);
%         grid on;
%         title("wet day climatology, equalized");
%     %   legend("obs","mdl","obs-mdl", "mdl_excess");
%        legend("obs","mdl","obs-mdl");
%        xlabel("day of yr");
%        ylabel("prob of wet day");
    end
    
%             % this shouldn't happen, because it should have been trapped by the if(diff_count_per_year<=0)...end above.
%     if (round(mdl_excess) <= 0)      % fewer model precips (climatologically) than obs precips, so just bail out.   
%         excess_thresh = zeros(yrlen,1);
%         excess_frac = zeros(yrlen,1);       % this is wrong, Ian
%         return;
%     end

    daily_keep_frac = find_daily_excess(mdl_clim, obs_clim, obs_counts, mdl_counts);    % fraction of precips to keep on each individual day.
%   mask = daily_keep_frac < 1.0;       % binary mask for days we can adjust the fractions on.  We don't want to 
                                        % fiddle with days where we need to keep 100% of the points.
    excess_frac = 1-daily_keep_frac;

%     mdl_precip_prob_surf = find_precip_probs(mdl_hist, mdl_clim);       % this will need to adjust for nan's as well.
%     excess_thresh     = calc_excess_precip_threshold(daily_keep_frac, mdl_precip_prob_surf, bins_mdl, prcp_min);
%     
    mdl_precip_prob_surf = find_precip_probs2(mdl_hist);                 % this will need to adjust for nan's as well.
    excess_thresh = calc_excess_precip_threshold(daily_keep_frac, mdl_precip_prob_surf, bins_mdl, mdl_prcp_min);
%   ethresh = excess_thresh;
%   excess_thresh = climatology(excess_thresh,6,2);     % this is now done in calc_excess_precip_threshold...
    
    
%     prcp_mdl = reshape(prcp_mdl,[],1);
%     ethresh = repmat(excess_thresh, nyrs_mdl,1) - fudge_mdl;
%     prcp_mdl(prcp_mdl < ethresh) = 0;

    if (fignum > 0)
        surf_threshy = zeros(yrlen,1);
        surf_threshz = zeros(yrlen,1);
        for i=1:yrlen
            ix = find(bins_mdl>excess_thresh(i),1);
            if (isempty(ix)), ix=1; end
            surf_threshy(i) = bins_mdl(ix);
            surf_threshz(i) = mdl_precip_prob_surf(i, ix);
        end
        subplot(2,1,2);
        surf(bins_mdl, 1:yrlen, mdl_precip_prob_surf, "edgecolor","none");
        hold on;
        plot3(surf_threshy, 1:yrlen, surf_threshz,'r-',"linewidth",2);
%       view([143, 30]);
        lights_on;
        hold off;
        if (strlength(lbl2)>0 && ~strcmp(lbl2,"test"))
            figname=sprintf("/Volumes/lacie_1/data/prcp_figs/calc_prcp_thresh_%s.png", fix_filename(lbl2));
            h=figure(fignum);
            saveas(h, figname);    
        end
    end
    
%     no_trim_mask = daily_keep_frac == 1;
%     if (sum(no_trim_mask) > 0)
%         notrim_obs = prcp_obs(no_trim_mask,:);
%         notrim_mdl = prcp_mdl(no_trim_mask,:);
%         missing = sum(notrim_obs(:)>0) - sum(notrim_mdl(:)>0); 
% %       wet_frac_obs = (sum(prcp_obs(:)>0)-missing)/npts_obs/nan_fact_obs;
%     end
%     
%     
%     excess_thresh = remap_time_dimension(excess_thresh, yrlen, [], dayPos);
%     excess_frac   = remap_time_dimension(excess_frac,   yrlen, [], dayPos);
end

function   hsurf = find_precip_probs2(mdl_hist)
%   returns a reverse CDF:  probability of 1 is at 0, probability of 0 is at max precip.

    hsurf = rect_filt(mdl_hist, 42, 1, true);
    
    hsurf = cumsum(hsurf, 2, "reverse");
    yrlen = size(hsurf,1);
    for i=1:yrlen
        hsurf(i,:) = hsurf(i,:)/hsurf(i,1);
    end
end

function [keep_frac, pos_ix] = find_daily_excess(mdl_clim, obs_clim,  obs_counts, mdl_counts)
    %
    %   mdl_excess is # of excess wet days/year that the model produces.
    %
    % returns daily fraction of model's precip events that should be kept.
    % If all the model's daily counts are greater than the observation's daily counts, then we simply keep obs_counts/mdl_counts
    % If there are times of the year where the model produces fewer precip events than the observations (obs_counts/mdl_counts > 1) 
    %   ** change, 4/13/22
    % We need to only trim from the areas where there are excess events.
    %
    %   
    % obsolete:  then we need to include some additional precips in the positive area.
    %
    %   fraction to keep is based on the difference between the wet-day climatologies.  But this can produce a wet-day
    %   total that is slightly different (generally within 1 percent) from the actual counts.  
    %   So we adjust that with a scale factor.
    %   
    yrlen = size(mdl_clim,1);
    if (sum(mdl_clim) < sum(obs_clim))      % If more obs than mdl wet days, keep everybody.
        keep_frac = ones(yrlen,1);
        pos_ix = false(yrlen,1);
        return;
    end

    keep_frac = min(1.0,obs_clim./mdl_clim);
    excess_frac = max(0,1-keep_frac);
    pos_ix = excess_frac > 0;
    
    if (~all(pos_ix))
        fprintf("here!\n");
    end

    est_excess = sum(excess_frac .* mdl_clim);
    mdl_excess = sum(mdl_counts(pos_ix)) - sum(obs_counts(pos_ix));
    
    scale_factor = mdl_excess/est_excess;       % adjustment to match the exact amount needed.
    keep_frac = 1 - excess_frac * scale_factor;
    
%                 %ss=(sum(y1(pos))-sum(diff(~pos)))/sum(y1(pos))
%     if (all(pos_ix))
%         keep_frac = 1 - excess_frac * scale_factor;
%     else
%             % a bit more work if some days don't have enough precips.  Then we need
%             % to distribute the excess over only those days with too many precips.
%         difs = mdl_clim - obs_clim;
%         posix = difs > 0;                         % region where we can add in extra precips.
% 
%         negix = ~posix;     % the days where mdl counts are greater than 
%         neg_area = sum(abs(difs(negix)));
%         pos_area = sum(    difs(posix));
%         scale_factor = scale_factor * (pos_area+neg_area)/pos_area;        
% 
%         excess_frac = excess_frac * scale_factor;
%         keep_frac = 1-excess_frac;
%         keep_frac(posix) = min(1.0, keep_frac(posix) * scale_factor);
%   end    
end
