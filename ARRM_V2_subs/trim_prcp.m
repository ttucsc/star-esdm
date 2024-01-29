function [DA_obs, DA_mdl, DA_hist, obs_fut, all_anoms,excess_thresh] = trim_prcp(DA_obs, DA_mdl, DA_hist, obs_fut, all_raw_data, siteLbl, show_figs, gridptLbl, lbl2)
%   Trims precip so DA_obs & DA_hist have same number of precip events.
%   DA_mdl is trimmed at the same threshold as DA_hist.,
%   Also sets all zeros to NAs so they won't be included in histogramming.
%   1.  obs precip is truncated at obs_prcp_min. i.e. obs precip < obs_prcp_min is set to 0.
%   2.  mdl precip is thresholded so it matches in count to obs, and then
%       mdl precip is shifted so min value is obs_prcp_min.
%       That way, identical scaling can be used on both obs & model data. 
%       (Except we don't use identical scaling, now. But since mdl & hist is shifted to obs prcp_min, we update DA_obs
%       and DA_hist's prcp_min.
%
%   3.  DA_obs, DA_mdl & DA_hist's zero_maps are updated to reflect all the zeros.

% need:
%   1.  trim model precip @ 1/2 or 1/4 of obs' min precip.
%^  2.  histogram-equalize obs precips to get appropriate counts in each "day".  
%           use same compression to equalize model precips
%           probably use same rules for setting equalizing yrlen.

    yrlen = DA_obs.RP.yrlen;
    obs_prcp_min = DA_obs.DP.prcp_min;
    mdl_prcp_min = DA_mdl.DP.prcp_min*.25;    % we'll allow the trimming process to reach down to 1/4 the mdl prcp_min if necessary to arrive at the right number of precip events.
%   pdf_yrlen = DA_obs.RP.pdf_yrlen;
    
%   fudge = 1e-7;   % to avoid any rounding error problems.
%   tot_prcp_events = sum(DA_obs.anoms >= obs_prcp_min - fudge); 
%   pdf_yrlen = min(DA_obs.RP.pdf_yrlen, DA_mdl.RP.pdf_yrlen);
%   pdf_yrlen = 128;
%   pdf_yrlen = floor(tot_prcp_events/30);      % We use a less conservative pdf_yrlen here than we do for the main downscaling.
                                                % Aim for at least 30 precip events in each compressed day.    
    
%    [excess_thresh, excess_frac, excess_tot_frac, obs_clim] = calc_excess_prcp(DA_obs.anoms, DA_hist.anoms, obs_prcp_min, yrlen, true);

    if (strlength(lbl2)>0)
        [excess_thresh, excess_frac, tot_wet_frac, obs_clim] = calc_excess_prcp(DA_obs.anoms, DA_hist.anoms, obs_prcp_min, mdl_prcp_min, yrlen, true, 101, lbl2); %, pdf_yrlen
    else
        [excess_thresh, excess_frac, tot_wet_frac, obs_clim] = calc_excess_prcp(DA_obs.anoms, DA_hist.anoms, obs_prcp_min, mdl_prcp_min, yrlen, true); %, pdf_yrlen
    end
   % do we need to trim DA_obs.anoms still, ian?
    
    if (show_figs)
        fignum=[61,62,63,64];
    else
        fignum=[0,0,0,0];
    end
    fidlog = DA_obs.DP.fidlog;
%   prcp,        = trim_excess_prcp(prcp,        excess_thresh, excess_frac, wet_frac,        exact, prcp_min,     yrlen, do_shift, ignore_nans, fignum,    obs_clim, maxtries, fidlog, lbl,             lbl2)
    DA_obs.anoms = trim_excess_prcp(DA_obs.anoms, obs_prcp_min,   [],          [],            false, obs_prcp_min, yrlen, true,     true,        fignum(1), [],       [],       fidlog, "obs "+gridptLbl, "obs"+lbl2); % clears out below prcp_min and shifts to where we need it.   
    if (~isempty(obs_fut))
        obs_fut     = trim_excess_prcp(obs_fut,   obs_prcp_min,   [],          [],            false, obs_prcp_min, yrlen, true); % clears out below prcp_min and shifts to where we need it.   
    end
%   fo_tot_frac = sum(obs_fut>0)/length(obs_fut);
zobs=reshape(DA_obs.anoms,365,[]);
[zhist1, ~, ~] = trim_excess_prcp(DA_hist.anoms, excess_thresh, excess_frac, tot_wet_frac,    false,  [mdl_prcp_min, obs_prcp_min], yrlen,     true,        true,      fignum(2), obs_clim,  [], fidlog, "      hist "+gridptLbl, "hist_"+lbl2);

%   [prcp,          excess_thresh, excess_frac, ntries] = trim_excess_prcp(prcp,          excess_thresh, excess_frac, excess_tot_frac, exact, prcp_min,                     yrlen,     do_shift, ignore_nans, fignum, maxtries)
    [DA_hist.anoms, ethresh,       etotfrac           ] = trim_excess_prcp(DA_hist.anoms, excess_thresh, excess_frac, tot_wet_frac,    true,  [mdl_prcp_min, obs_prcp_min], yrlen,     true,        true,      fignum(2), obs_clim,  [], fidlog, "      hist "+gridptLbl, "hist_"+lbl2);
%   [DA_mdl.anoms, ethresh2,       etotfrac2          ] = trim_excess_prcp(DA_mdl.anoms,  excess_thresh, excess_frac, excess_tot_frac, true,  [mdl_prcp_min, obs_prcp_min], yrlen,     true,        true,      63, obs_clim); %#ok<ASGLU>
%   [all_anoms,           ~,               ~          ] = trim_excess_prcp(all_anoms,     excess_thresh, excess_frac, excess_tot_frac, true,  [mdl_prcp_min, obs_prcp_min], yrlen,     true,        true,      63, obs_clim);
    [DA_mdl.anoms, ethresh2,       etotfrac2          ] = trim_excess_prcp(DA_mdl.anoms,        ethresh, etotfrac,    tot_wet_frac,    false, [mdl_prcp_min, obs_prcp_min], yrlen,     true,        true,      fignum(3), obs_clim,  [], fidlog, "       mdl "+gridptLbl, " mdl_"+lbl2); %#ok<ASGLU>
    [all_anoms,           ~,               ~          ] = trim_excess_prcp(all_raw_data,        ethresh, etotfrac,    tot_wet_frac,    false, [mdl_prcp_min, obs_prcp_min], yrlen,     true,        true,      fignum(4), obs_clim,  [], fidlog, "all precip "+gridptLbl, " all_"+lbl2);
zhist1=reshape(zhist1,365,[]);
zhist2=reshape(DA_hist.anoms,365,[]);
zmdl=reshape(DA_mdl.anoms,365,[]);
ssns = {[335:365,1:59]; 60:151; 152:243; 244:335; 1:365};
DA_hist.DP.print_log("prcp_counts (obs,hist1,hist_match,mdl):  %s :", gridptLbl)
for i=1:5  
    DA_hist.DP.print_log(" %2d %6d %6d %6d %6d:", i, sum(zobs(ssns{i},:)>0,"all"), sum(zhist1(ssns{i},:)>0,"all"), sum(zhist2(ssns{i},:)>0,"all"), sum(zmdl(ssns{i},:)>0,"all")); 
end
DA_hist.DP.print_log("\n")

        % we've shifted hist & mdl's precip so their min values match the model's prcp_min, so update their prcp_min
    DA_hist.DP.prcp_min = obs_prcp_min;
    DA_mdl.DP.prcp_min  = obs_prcp_min;
    
    nyrs_obs  = length( DA_obs.anoms)/yrlen;
    nyrs_hist = length(DA_hist.anoms)/yrlen;
    nyrs_mdl  = length( DA_mdl.anoms)/yrlen;
    nyrs_fobs = length(obs_fut)/yrlen;
    wd_obs  = sum( DA_obs.anoms > 0)/nyrs_obs;
    wd_hist = sum(DA_hist.anoms > 0)/nyrs_hist;
    wd_mdl  = sum( DA_mdl.anoms > 0)/nyrs_mdl;
    wd_fobs = max(0, sum(obs_fut > 0)/nyrs_fobs);
    fprintf("total wet days/year:  obs  %5.1f   hist %5.1f   mdl %5.1f   fobs %5.1f   %s\n", wd_obs, wd_hist, wd_mdl, wd_fobs, siteLbl)
    DA_obs.zeros_map  = DA_obs.anoms == 0; 
    DA_hist.zeros_map = DA_hist.anoms == 0;
    DA_mdl.zeros_map  = DA_mdl.anoms == 0;
    
        % move anoms back into raw_fix so we can use it to calculate theclimatology  (for precip amount, not wet-day climatology) 
    DA_obs.raw_fix  = DA_obs.anoms;
    DA_hist.raw_fix = DA_hist.anoms;
    DA_mdl.raw_fix  = DA_mdl.anoms;
       
        % now set all the zeros to NAs.
    DA_obs.set_zeros_to_NAs("anoms");
    DA_mdl.set_zeros_to_NAs("anoms");
    DA_hist.set_zeros_to_NAs("anoms");
    all_anoms(all_anoms<=0) = nan;
end

