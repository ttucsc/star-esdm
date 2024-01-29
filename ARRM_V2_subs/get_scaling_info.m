function  scale_info = get_scaling_info(DA, distrib, prcp_min,  DA_lbl, stnID, fignum, mnp, offset, fidlog, lbl2)
%function scale_info = get_scaling_info(DA_obs, DA_mdl, DA_hist, prcp_min, stnID, fignum) %#ok<INUSL>
%   Finds best precip scaling power to use to get as close to a gaussian as possible.
%
%   Joins all the input raw data together and finds the best scaling overall.  This is because we need the scaling to be
%   identical across all the inputs.  This isn't needed for mapping between distributions. But for outliers, where we
%   will use an alternate mapping based on model precip and ratio of obs/model 2-sigma-equivalent points, we need the
%   scaling to be the same.
%
%   assumes 

    if (~exist("lbl2","var")), lbl2 = ""; end
    
    nbins = 1500;
    nsigmas = 10;
%   offset = 0;
    if (strcmp(distrib,"log"))
        min_pwr = distrib;
    else
        min_pwr = 0.2;
    end
    figlbl = sprintf("%s, %s", DA_lbl, stnID);
    if (~exist("mnp","var")), mnp = []; end
    
    if (~exist("fidlog","var")), fidlog = []; end

    if (isa(DA, "ARRM_V2_DisaggregateSignal"))
        prcp  = DA.anoms(~isnan(DA.anoms));
    else
        prcp  = DA;
    end
    
    min_prcp = min(prcp(prcp>.01));
    try
        prcp_min = min(prcp_min, min_prcp);
    catch
        fprintf("oops get_scaling_info(%s)\n", DA.DA_type);
    end
%     obs_prcp  = DA_obs.anoms(~isnan(DA_obs.anoms));
%     hist_prcp = DA_hist.anoms(~isnan(DA_hist.anoms));     % PROBLEM:  should avoid using duplicate data here, Ian!
%     mdl_prcp = DA_mdl.anoms(~isnan(DA_mdl.anoms));
%     hist_prcp = DA_hist.anoms(~isnan(DA_hist.anoms));     % PROBLEM:  should avoid using duplicate data here, Ian!
%     all_prcp = [obs_prcp; mdl_prcp; hist_prcp];
%   [    pwr_best, found, skew_best, kurt_best, sig_best,   mu_best,      offset, bins, hcnts] = prcp_kurtosis_search_2(prcp, prcp_min, min_pwr, lbl, offset, fignum, mnp)
%   [prcp_scaling, found, prcp_kurt, prcp_sigma] = prcp_kurtosis_search(all_prcp, prcp_min, nsigmas);
%   [prcp_scaling, ~,     ~,         prcp_sigma, ~,      prcp_offset      ] = prcp_kurtosis_search_1(all_prcp, prcp_min, min_pwr, stnID, offset, fignum);       % min_pwr set to 1/5 (fifth root).  (using default offset of .2.)
    try
%       [prcp_scaling,    found,        ~,          ~, prcp_sigma, prcp_mu, prcp_offset      ] = prcp_scaling_search(    prcp, prcp_min, min_pwr, "skew", 1, offset, figlbl, fignum, mnp);       % min_pwr set to 1/5 (fifth root).  (using default offset of .2.)
        [prcp_scaling,        ~,        ~,          ~, prcp_sigma, prcp_mu, prcp_offset      ] = prcp_scaling_search(    prcp, prcp_min, min_pwr, "skew", 1, offset, figlbl, fignum, mnp, fidlog, lbl2);       % min_pwr set to 1/5 (fifth root).  (using default offset of .2.)
    catch me
        oops("failed prcp_scaling_search", me);
    end
%   prcp_scaling = 1/prcp_scaling;  % we'll hang on to scaling as its reciprocal, so 2 means we're using square root, etc.

%     if (~found)
%         [prcp_scaling,    found,        ~,          ~, prcp_sigma, prcp_mu, prcp_offset      ] = prcp_scaling_search(    prcp, prcp_min, "log", "",[.1,.9], 1e-8, figlbl, fignum, mnp);       % min_pwr set to 1/5 (fifth root).  (using default offset of .2.)
%     end   

    edge1 = prcp_mu - nsigmas * prcp_sigma;
    edge2 = prcp_mu + nsigmas * prcp_sigma;
    edges = linspace(edge1, edge2, nbins+1);    % bins going from scaled mean to +/-10 sigmas.
    dx = edges(2)-edges(1);
    zix = find(edges>=0,1);
%   edges = edges + dx/2;                    % shift by 1/2 delta so edges will straddle 0.
    lin_edges = rescale_data_1(edges(zix:end), prcp_scaling, prcp_min, "reverse", prcp_offset);   % edges scaled back in mm    
    linear_edges = [repmat(lin_edges(1),1,zix-1), lin_edges];
    linear_bins  = (linear_edges(1:end-1)+linear_edges(2:end))/2;
    bins = rescale_data_1(linear_bins, prcp_scaling, prcp_min, "forward", prcp_offset);    
%     edges = [fliplr(-edges), edges]; % reflect edges so we can make a vaguely gaussian distribution later by reflecting the precipitation.
%     bins =  [fliplr(-bins),0,bins];
            % edges now straddle 0, going from -binmax to +binmax
       
%     linear_edges = [fliplr(-linear_edges), linear_edges];        % do I need or want linear bins as well or instead?  And should the linear bins be centered on the linear edges?
%     linear_bins  = [fliplr(-linear_bins), 0, linear_bins];
        % to scale precip data:
        %   1.  subtract prcp_min, & threshold at 0  (i.e., any thing less tha 0 --> nan)
        %   2.  add offset (default 0.2 is good) to shift slightly away from 0 before scaling.
        %   2.  raise to scaling_pwr
        %   3.  subtract offset^scaling_pwr (currently using default, of .2). so prcp_min is scaled to 0.
    scale_info = struct("prcp_scaling",prcp_scaling, "sigma",prcp_sigma,"edges",edges, "bins", bins, "dx", dx, "prcp_offset", prcp_offset, "prcp_min", prcp_min, "nbins", nbins, "nsigmas", nsigmas, "linear_edges", linear_edges, "linear_bins", linear_bins);
        % to scale back to regular precip values:
        %   1.  add offset^scaling_pwr
        %   2.  raise to 1/scaling_pwr
        %   3.  subtract offset, and add prcp_min back so original prcp_min's are scaled back to prcp_min.
        
%       disp(scale_info);
     
end

