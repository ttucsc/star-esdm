function [DA, RP] = ARRM_V2_merge_DAs(DAs, RP, DP, do_pdfs, ttl, merge_seed)
% moerges several DAs together to  create a single data stream.  Usually the DAs are either 4 gridcells being
% interpolated to a point inside the rectangle, or several weather stations in close vicinity.
%
%   Creates a weighted average of the base and rolling histograms and climatology, and does a weighted random selection
%   of the anomalies.  Uses the same random seed (1) every time, so the results will be identical if the run is repeated
%   with the same data.
%
%   Inputs:
%       DAs     two or more DisAggregations
%       RP      RunParams object
%       DP      DataParams object
%       do_pdfs boolean.  If true, calculates base and rolling pdfs and CDFs on merged.
%       ttl     label or identifier for the output DA.  (Used by plotting software for titles)
%       merge_seeed
%               integer between 0 & 2^32-1.  should be deterministic if you want repeatable results.
%                   recommendation:  floor(lat * 100 * lon * 100)
%
%   Outputs
%       DA      merged DisAggregation
%       RP      updated RunParams.
%
%   NOTE:  merged histograms, climatology, pdfs, CDFs, etc., use combined inputs, while anomalies are a random sampling of the
%          anomalies from all the inputs.  Therefore if you reassemble the output signal and DisAggregate it, you will
%          get a slightly different set of values.  
%          Reason:  The merged histograms, climatology, etc. are based on a larger set of points.

    if (~exist('ttl','var')), ttl='ARRM_V2_merge_DAs'; end
    
    nsets = numel(DAs);      % if data only 1-D or 2-D, size(data,3) is 1.
    
    if (nsets == 1)
        if (iscell(DAs))
            DA = DAs{1};
        else
            DA = DAs;
        end
        return
    end
    
    DAs = reshape(DAs, 1,nsets);   % take it from a grid of size 2 x 2 to a matrix of data, size npts x 4.
    if (numel(RP.weights)==1)
        RP.weights=ones(1,nsets);
    elseif (numel(RP.weights) ~= nsets)
        error('error:  weights dimensions don''t match size of DAs');
    end
    RP.weights = reshape(RP.weights, 1, nsets);
                                                        % scale weights so histograms represent true number of data points.
    histogram_scaling = 1/max(RP.weights);              % when combining histograms, we want the counts to reflect all
                                                        % the data points so we can calculate the Silverman-Sigmas
                                                        % properly.
                                                        
    
    DA= ARRM_V2_DisaggregateSignal([],RP, DP, ttl); % create an empty DisAggregation object.
    
        % weighted merge of the various components.
    [DA.anoms, DA.na_map_in, DA.grid_pt_used] = weighted_select(DAs, RP.weights, true, 'anoms', merge_seed);      % you should seed with a repeatable number so we always get the same weighted selection...
    DA.na_map_unknown = false;
    [bstart, bend] = DP.using_range('base_yrs');
    DA.lohi_base = [sum(DA.anoms(bstart:bend)<RP.edges(1)),sum(DA.anoms(bstart:bend)>RP.edges(end))];
    DA.avg_base           = weighted_sum(DAs, RP.weights, 'avg_base');
    DA.base_clim          = weighted_sum(DAs, RP.weights, 'base_clim');
    DA.base_hist          = weighted_sum(DAs, RP.weights, 'base_hist') * histogram_scaling;  % scale up so we can calculate the Siverman-Sigmas properly.
    DA.moving_clim        = weighted_sum(DAs, RP.weights, 'moving_clim');
    DA.moving_clim_1D     = weighted_sum(DAs, RP.weights, 'moving_clim_1D');
    
    %   recalculate the phase shifts for base and moving climatologies.
    %   recalc rather than weighted sum, because phase may not match w/ weighted sum.
    
    DA.phase_clim = calc_phase(DA.base_clim, RP.clim_nterms, 10, RP.yrlen);
    
    mycos = cos(2*pi*(0:DA.RP.yrlen-1)/RP.yrlen);
    moving_clim = reshape(DA.moving_clim, RP.yrlen, []);
    nyrs = size(moving_clim,2);
    DA.phase_moving = nan(nyrs,1);
    for i=1:nyrs
        y = mycos * range(moving_clim(:,i)/2) +  mean(moving_clim(:,i));
        DA.phase_moving(i) = calc_phase(moving_clim(:,i), RP.clim_nterms, 10, RP.yrlen, y);
    end


    
    if (~isempty(DP.rolling_yrs))
        [rstart, rend, nyrs] = DP.using_range('rolling_yrs');
        DA.rolling_clim = mean(reshape(DA.moving_clim(rstart:rend),RP.yrlen,nyrs),2);     
        DA.rolling_hist = weighted_hist(DAs, RP.weights,'rolling_hist') * histogram_scaling;
        DA.rolling_clim = weighted_sum(DAs, RP.weights, 'rolling_clim');
    end

        % recreate trend_params by calculating trend from weighted trend.
        % calc_trend_params doesn't remove any mean, so we can use the weighted mean calculated above.
    DA.trend = weighted_sum(DAs, RP.weights, 'trend');
%    trend_order = 3;    % kludge for now...
%    for i=1:nsets; trend_order = max(trend_order, DAs{i}.trend_order); end
%    DA.trend_params = calc_trend_params(DA.trend, trend_order, RP.yrlen);
%   DA.trend_params.avg_trend = DA.avg_trend;
    DA.lohi_all = [sum(DA.anoms(:)<RP.edges(1)),sum(DA.anoms(:)>RP.edges(end))];
    DA.raw_data=DA.reassemble();  
%     DA.merge_data = cell(nsets,1);
%     for i=1:nsets
%         DA.merge_data = DAs{i}.raw_data;
%     end        

    [tr_start, tr_end] = DP.using_range('trend_yrs');
    trdata = DA.raw_data(tr_start:tr_end) - DA.avg_base;
    DA.trend_params = DA.calc_trend_params(trdata);
    if (DA.RP.do_daily_trends)
        DA.trend_params_daily = DA.calc_daily_trend_params();
    end
        % calculate pdfs on merged dataset.
    DA.calc_pdfs(do_pdfs);   % s/b 3 long, for do_base, do_rolling and do_problines.
    
    if (RP.do_base_probs)
        DA.calc_base_probs();
    end
    if (RP.do_rolling_probs)
        DA.calc_rolling_probs();
    end

end

