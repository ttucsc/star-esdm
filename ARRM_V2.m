 function [mapped_output, DA_out, DA_obs, DA_hist, DA_mdl, DA_redisaggregated] = ARRM_V2( DA_obs, DA_hist, DA_mdl, varargin)
%function [mapped_output, DA_out, DA_obs, DA_hist, DA_mdl, DA_redisaggregated] = ARRM_V2( DAs_obs, DAs_hist, DAs_mdl, varargin)
%
%       ARRM_V2
%   
%       This version:  receives 3 DisAggregateSignal objects:  Obs, Hist and Model.  Hist & Model may contain some of
%       the same data (if downscaling in the historical period as well as future), but for simplicity are separate
%       objects.  (Originally code used a single DA for both...just made things too complicated.)
%       Hist and Model will generally be resampled from larger gridded data. 
%       Calling code (ARRM_V2_run_gridded or ARRM_V2_run_stations) will read gridded data, find four closest gridcells, and 
%       resample/interpolate into a single data stream and set of histograms.  
%       Hist data is used for the base period's histograms/pdfs/CDFs;  model data for the rolling histograms/pdfs/CDFs.  
%       Obs data will either be station data or gridded observation data, presumably covering the same base period as is
%       used for the Hist data.
%
%       Inputs:
%           input data streams:  can be cell arrays if more than 1 data set for obs, hist or mdl.
%                                   (for mdl, will generally be 4 nearest model gridpoints)
%
%               DAs_obs      Observation data.  will use base data
%               DAs_hist     Historical model data.  (may be empty, if hist is in model DA).  will use base data
%               DAs_mdl      future model data. will use rolling data.  Will downscale using 
%                                                       DA.RP   run params
%                                                       DA.DP   data params
%                               each DA contains its own RP and RP, but they should have the same parameters.
%                               Only differences will be in the data years.
%
%           optional kwd/value pairs (in varargin):
%
%               "title", ttl                title for run
%               "obs_weights",  obs_wts     weights for merging obs, hist or model DAs.
%               "hist_weights", hist_wts        provide these when there are more than 1 DA for any of                                       
%               "mdl_weights",  mdl_wts         the input datasets
%               "do_histograms", do_histograms  boolean.  If true, recalculate histograms even if present.
%                                                   can be 3 elements long (one each for obs, hist & mdl)
%
%                       See  ARRM_V2_RunParams.m and ARRM_V2_DataParams.m for descriptions of the RP and RP parameters
%
%       Outputs:
%         Production run:
%           mapped_output       bias-corrected data values.  (Also in DA_out struct, as DA_out.mapped).
%           DA_out              DisAggregateSignal object with mapped_output:
%               DA_out.mapped_output       = moving_clim_adjusted(:) + trend_model + avg_model + anom_mapped;
%               DA_out.lohi                = counts of # of anomalies outside mapping range (usually -50 to +50 C)
%               DA_out.na_counts           = counts of # of NAs in each dataset & # of model NA's mapped to NA in mapped results via mapping process
%                                       Note:  if mapped na_count is more than a few, then an analysis should be done to see
%                                       where the NAs are occurring on the globe, and determine the cause, as it
%                                       probably indicates a weakness or problem in the mapping algorithm.
%                                       Note 2:  if obs, hist, model or futobs NAs exist, analysis should be done to
%                                       make sure these are not biasing the results.
%                                       Note 3:  NAs in the input data are initially replaced with an estimate of their
%                                       likely value, in order to have a complete set for working in fourier domain.
%                                       However, the NAs are excluded from the histogramming/pdf creation so as not to
%                                       bias the pdfs.
%           DA_redisaggregated  re-analyzed disaggregation set created from the mapped output.
%                                       will be very close, but not identical, to DA_out
%
%           DA_obs, DA_hist, DA_mdl     Original inputs, but with PDFs and CDFs calculated, if they weren't already done
%                                       on input.

%______________________________________________________________________________________________________________________
%                    Math Notes:
%
%   If range(y_src) > possible range of new values to be mapped, then no checking for out-of-range indexes is needed.
%
%       Note:  thresh is a cumulative probability:
%               Sigma       Thresh      single-tailed likelihood
%               6           9.86e-10    (~1 in 1 billion)
%               5           2.87e-7     (~1 in 3.5 million)
%               4           3.17e-5     (~1 in 30,000)
%               3.09        1.0e-3      ( 1 in 1000)
%               3           1.35e-3     (~1 in 750)
%               2.5         .00621      (~1 in 150)
%               2.3         .01         ( 1 in 100)
%               2           .0228       (~1 in 40)
%               1           .159        (~1 in 6)
%
%       When checking for unrealistic output results, in theory, a reasonably safe threshold to use is 1e-4 (~4.3 sigma).  
%           1 in 58000 will fall in the tails; on average, 1 point every 160 years will fall outside this range for 
%           a given day.
%           However, that would apply to a normal (gaussian) distribution;  certainly will not be true for precip, and
%           temperature data, which is fairly normal, is probably not normal out in the wings.
%

%______________________________________________________________________________________________________________________

    ARRM_V2_setpath(mfilename('fullpath'));
        
                                                                        % parse and verify that data and run parameters are all consistent. 
    [DA_obs, DA_hist, DA_mdl, RP, DP, ttl, do_histograms] = init_params(DA_obs, DA_hist, DA_mdl, varargin{:});       
                                                                        % throws an exception if anything looks wrong.     
                                                                        
                % make sure everybody has been histogrammed.
                % will also calculate binning.
                % histogramming will also calculate the climatologies, fix NAs, etc.
    DA_histogram(DA_obs,  do_histograms(1), "DA_obs");
    DA_histogram(DA_hist, do_histograms(2), "DA_hist");
    DA_histogram(DA_mdl,  do_histograms(3), "DA_mdl", DA_mdl.RP.pdf_yrs);
    
    delta_hist_mdl = DA_mdl.avg_base - DA_hist.avg_base;    % offset between hist and model.  Non-zero only if model data doesn't overlap hist data.
    
    
                % make sure we have PDFs
    if (isempty(DA_obs.CDF_base)),     DA_obs.calc_pdfs(  true, false); end   % calculate pdf & CDF for base period only
    if (isempty(DA_hist.CDF_base)),    DA_hist.calc_pdfs( true, false); end   % base period only for hist
    
        % s/b ~isPrecipRun?
%   if (any(strcmp(DA_mdl.RP.prcp_distrib,["pwr","log"])))        % if not doing a precip run, then calculating binning from range of data.  Use at least 500 bin.
    if (~DA_mdl.RP.isPrecipRun)
        %   get edges from hist disaggregation and use for model
        rr = .25*range(DA_hist.RP.out_edges);
        ee = [DA_hist.RP.out_edges(1)-rr, DA_hist.RP.out_edges(end)+rr];
    %   ee = mean(DA_hist.RP.out_edges) + [-rr/2,rr/2];
        dx = min(.1, range(ee)/500);
        nedges = ceil(range(ee)/dx-1e-12)+1;
    %   if (mod(nedges,2)==0), nedges = nedges+1; end
        DA_mdl.RP.out_edges = linspace(ee(1), ee(2), nedges);
        
    end        
    
    if (isempty(DA_mdl.CDF_rolling)),  DA_mdl.calc_pdfs( false, true);  end   % rolling pdfs only for mdl.
    
            % this should probably be done in ARRM_V2_Disaggregate at end of calc_pdfs?
            % and do we need to do this for DA_obs or DA_hist?
    if (DA_obs.sigma_normalized)
        DA_obs.anoms = DA_obs.normalize_anoms_by_sigmas(DA_obs.anoms, DA_obs.DP.yrlen, false, DA_obs.norm_sigmas);
        DA_obs.sigma_normalized = false;
%       [obj.anoms, obj.norm_sigmas] = obj.normalize_anoms_by_sigmas(obj.anoms, obj.DP.yrlen, true, []);
    end
    if (DA_hist.sigma_normalized)
        DA_hist.anoms = DA_hist.normalize_anoms_by_sigmas(DA_hist.anoms, DA_hist.DP.yrlen, false, DA_hist.norm_sigmas);
        DA_hist.sigma_normalized = false;
%       [obj.anoms, obj.norm_sigmas] = obj.normalize_anoms_by_sigmas(obj.anoms, obj.DP.yrlen, true, []);
    end
    if (DA_mdl.sigma_normalized)
        DA_mdl.anoms = DA_mdl.normalize_anoms_by_sigmas(DA_mdl.anoms, DA_mdl.DP.yrlen, false, DA_mdl.norm_sigmas);
        DA_mdl.sigma_normalized = false;
%       [obj.anoms, obj.norm_sigmas] = obj.normalize_anoms_by_sigmas(obj.anoms, obj.DP.yrlen, true, []);
    end
    
    
                % make sure we have the anomaly probabilities
                % we don't need obs & hist anomaly probs for downscaling, but we want to log the number of NAs and
                % outliers, so we calculate them here.
    if (isempty(DA_obs.anom_probs_base)),    DA_obs.calc_base_probs(); end
    if (isempty(DA_hist.anom_probs_base)),   DA_hist.calc_base_probs(); end
    if (isempty(DA_mdl.pthresh_base))
        DA_mdl.pthresh_base = DA_hist.pthresh_base;
    end
    if (isempty(DA_mdl.anom_probs_rolling)), DA_mdl.calc_rolling_probs(); end
    
                % for precips, put zeros back in the anomalies and anomaly probabilities.
    if (DP.isPrecipRun)
        DA_mdl.anoms(DA_mdl.zeros_map) = nan;               % only want to work on wet days, so set all zero-precip days to NAs.
        DA_mdl.anom_probs_rolling(DA_mdl.zeros_map) = nan;  % 
        DA_mdl.anoms(DA_mdl.na_map_in) = nan;
        DA_mdl.anom_probs_rolling(DA_mdl.na_map_in) = nan;
    end        
        
    DA_out = create_output_DA(DA_obs, DA_mdl, RP, DP, ttl);             % start the DA object for the output by selecting the 
                                                                        % appropriate fields from inputs.
     
    if (~DP.isPrecipRun)                                                                
        delta_clim = (DA_obs.base_clim -DA_obs.avg_base) - (DA_hist.base_clim-DA_hist.avg_base);  
                                                            % difference between observations and historical model.
                                                            % NOTE:  the means have NOT been removed from the signal, so
                                                            % these are not zero-mean.  Thus we subtract the means from both
                                                            % to get the residual change.  

                    % adjust model climatology so it fits the observations during the historical period
                    % this subtracts the mean model climatology from moving_climatology and adds in the mean obs climatology
                    % Note:  this still needs to have the trend added back in;  we'll do that later.

        % DA_out.basic_clim = DA_obs.basic_clim;
        DA_out.moving_clim = adjust_climatology(delta_clim, DA_mdl.moving_clim);       % This is now reshaped to look like
                                                                        % the observations during the base period.
                                                                        % Moving_clim is detrended and zero-mean.
                                                                        % 
        DA_out.moving_clim_1D = adjust_climatology(delta_clim, DA_mdl.moving_clim_1D);    % do the same for the non-year_smoothed moving climatology.
    end
                % map the model anomalies to the observation probabilities, optionally adjusting for change between hist
                % CDF and model's rolling CDF.
                
    
    
    [mapped_anoms_im, mapped_nas, ~ ] = map_anomalies(DA_mdl, DA_obs, DA_hist); % get intermediate mapped anomalies.
    
    if (DP.isPrecipRun)
        mapped_anoms_im(DA_out.zeros_map) = 0;     % put the dry days back...
    end
    
    DA_out.insert_data(mapped_anoms_im,'mapped_anoms_im','rolling_yrs');
    if (DP.isPrecipRun)
        new_nas = mapped_nas & ~DA_mdl.zeros_map;
    else
        new_nas = mapped_nas;
    end
    nnew_nas = sum(new_nas);
    if (nnew_nas > 0)
%       DA_out.insert_data(new_nas,  'na_map_mapped',   'rolling_yrs');  adjust_outliers now will adjust both the far outliers and any mapped NAS.
        DA_out.DP.warn_log("warning:  %d mapped NAs\n", nnew_nas);
    end
%     else        
%         if (sum(mapped_nas) > 0)
%             DA_out.insert_data(mapped_nas,  'na_map_mapped',   'rolling_yrs');
%             DA_out.DP.warn_log("warning:  %d mapped NAs\n", sum(mapped_nas));
%         end
%     end
    
    [mapped_anoms, outlier_mapped_anoms] = adjust_outliers(DA_obs, DA_hist, DA_mdl, DA_mdl.anoms, mapped_anoms_im, new_nas);     % and fix the outliers.
    DA_out.insert_data(mapped_anoms,'mapped_anoms','rolling_yrs');
    DA_out.anoms = DA_out.mapped_anoms;     % replace anoms with the final anoms.
    DA_out.grid_pt_used = DA_mdl.grid_pt_used;
    DA_out.insert_data(outlier_mapped_anoms,'outlier_mapped_anoms','rolling_yrs');
            
    DA_out.mapped_output = DA_out.reassemble(DA_out.mapped_anoms, delta_hist_mdl);        % reassemble mapped anomalies, trend, avg_base & moving climatology into an adjusted time series.
    
        % for precip, mapping can map zeros to minisicule amounts (< 1e-12)
        % So we mask things back to zero using the zeros_map.
    if (DP.isPrecipRun)
        zeros_map = DA_out.zeros_map;
        DA_out.anoms(zeros_map) = 0;
        DA_out.mapped_anoms_im(zeros_map) = 0;
        DA_out.mapped_anoms(zeros_map) = 0;
        DA_out.mapped_output(zeros_map) = 0;
                % And make sure anything that mapped to below 0 gets set back to zero as well.        
        DA_out.anoms            = max(0, DA_out.anoms, "includenan");       % includenan so nan's stay nans. Default for max:  max(5, nan) is 5...
        DA_out.mapped_anoms_im  = max(0, DA_out.mapped_anoms_im, "includenan");
        DA_out.mapped_anoms     = max(0, DA_out.mapped_anoms, "includenan");
        DA_out.mapped_output    = max(0, DA_out.mapped_output, "includenan");
    end
            % mask off the original NAs back to NAs
    DA_out.anoms(DA_out.na_map_in) = nan;
    DA_out.mapped_anoms_im(DA_out.na_map_in) = nan;
    DA_out.mapped_anoms(DA_out.na_map_in) = nan;
    DA_out.mapped_output(DA_out.na_map_in) = nan;
    
            % finally, trim the output down to just the rolling years.
    DA_out.trim_data_yrs(DA_out.DP.rolling_yrs); 
            % and grab the mapped_output for returning to calling code.
    mapped_output = DA_out.mapped_output;

            % update the DP with all the proper filenames, etc.
    DA_out.DP = revise_output_DP(DA_out.DP, DA_obs.DP, DA_hist.DP, DA_mdl.DP);
    
            % This is  minor fudge.  We didn't actually stretch the obs. pdf's, so we don't have the true pdfs and problines
            % With a bit of work, we could actually stretch the obs probabilities to create a set of future-obs-hat 
            % distributions, rather than recalculating them.  
            % this would take more time, but would be more correct.            
            % Instead we recalculate the pdfs, which incorporates any noise in the resulting mapped values, rather than
            % having the actual PDFs that we used to generate the mapped values.
    if (DP.recalc_pdfs)
        DA_out.calc_pdfs(true, true);
        if (RP.do_problines)
            DA_out.calc_problines();
        end
    end
    if (nargout > 5)
        DA_mdl.DP.print_log("\n---Disaggregating downscaled output\n\n");
        
        DA_redisaggregated  = ARRM_V2_DisaggregateSignal(mapped_output, RP, DA_out.DP, ttl, false, false, false, false, false, false);
        DA_redisaggregated.prcp_scaling = DA_obs.prcp_scaling;
        DA_redisaggregated.calc_anomalies();
        DA_redisaggregated.set_zeros_to_NAs("anoms");
        DA_redisaggregated.run_disaggregation(false, true, true, false, true, true);
    end
 end
  
 function DA = create_output_DA(DA_obs, DA_mdl, RP, DP, runLbl)
    npts = (range(DA_mdl.DP.data_yrs)+1)*DA_mdl.DP.yrlen;           %same size as DA_mdl.
    DP_out = DP.update('dspType','downscaling');
    DA = ARRM_V2_DisaggregateSignal(npts,RP, DP_out, runLbl, false, false, false, false, false, false);       % create DA with blank raw_data array and appropriate flag arrays.
    DA.anoms            = DA_mdl.anoms;
    DA.orig_anoms       = DA_mdl.orig_anoms;
    DA.na_map_in        = DA_mdl.na_map_in;
    DA.na_map_unknown   = false;
    DA.zeros_map        = DA_mdl.zeros_map;
    
            % might want to copy some more pieces over, like moving_clim, moving_clim_1D, etc.
            % OR, always disaggregate results?
    DA.base_clim = DA_obs.base_clim;
    DA.base_hist = [];
    DA.avg_base  = DA_obs.avg_base;
    DA.trend = DA_mdl.trend;
    DA.trend_params = DA_mdl.trend_params;
    DA.rolling_hist = [];
    
    DA.raw_data  = DA_mdl.raw_data;
        % output DA's mapped data points will still have the same probabilities as the model data,
        % so these numbers are valid based on the weighted, merged PDFs.  Note, though, that if you redisaggregate the
        % output, the PDFs will change slightly, since DA_out will be based on a subset of the input to the downscaling
        % process, so these values will differ slightly after re-disaggregating.
    DA.anom_probs_rolling = DA_mdl.anom_probs_rolling;
    DA.anom_sdevs_rolling = DA_mdl.anom_sdevs_rolling;
    DA.outlier_map = DA_mdl.outlier_map;
    DA.far_outlier_map = DA_mdl.far_outlier_map;
    if (DP.isPrecipRun)
        DA.prcp_scaling = DA_obs.prcp_scaling;
    end
    
    DA.norm_sigmas = DA_obs.norm_sigmas;
    DA.sigma_normalized = DA_obs.sigma_normalized;
    
 end

 function DP = revise_output_DP(DP, DPobs, DPhist, DPmodel)   % updates data years and filenames for the output DP
 
    DP.obsnames = basename(DPobs.obsnames);
    DP.mdlnames = basename(DPmodel.mdlnames);
    DP.histnames = basename(DPhist.histnames);
    DP.obsdir = DPobs.obsdir;
    DP.mdldir = DPmodel.mdldir;
    DP.histdir = DPhist.histdir;
    DP.stninfo = DPobs.stninfo;
    DP.stnID = DPobs.stnID;
    DP.stnName = DPobs.stnName;
    
%    DP.trend_yrs = DP.rolling_yrs;
 end
 
 function     [mapped_anoms, outlier_mapped_anoms] = adjust_outliers(DA_obs, DA_hist, DA_mdl, anoms, mapped_anoms_im, new_nas)
 %  Returns:  mapped anomalies, with outliers adjusted via outlier-scaling.
 %  Also returns the outlier-scaled values for all anomalies, so we can run some statistics on them later.
 %
 %  Outlier scaling creates an alternate mapping based simply on the ratio of the outlier to the pthresh
 %  (acceptable-probability threshold).  This is essentially assuming that the observation and historical-model
 %  distributions are identically-shaped in the tails, except for a scaling factor.
 %
 %  The pthresh is simply the probability quantile line for the acceptable-probability threshold
 %  We scale each outlier based on it's ratio to the pthresh.  If a far-outlier's original anomaly (anoms()  is 
 %  50% beyond DA_mdl's pthresh line, then we scale it (relative to it's median) so it is 50% beyond the output's
 %  pthresh line.
 %
 %  For temperature data, which is zero-mean, this is a simple scaling.
 %  For precip data:  the anomalies are scaled differently between obs & model, and also not zero-mean.  So we have to
 %  rescale back to regular precip value, then scale into the obs' precip scaling space, then we can scale based on the
 %  pthresh ratios.  Then we subtract the median, scale based on pthresh's, and add back the median.
 %
 %  IAN:  might make sense to store the median line separately and have zero-median-based anomalies...That would make
 %  this routine much simplier.

    [ixstart,ixend] = DA_mdl.using_range('rolling_yrs');
    anoms = anoms(ixstart:ixend);   % we only want to work on the anomalies in the rolling years.
    npts = length(anoms);
    yrlen = DA_mdl.RP.yrlen;
%   pdf_yrstep = DA_mdl.RP.pdf_yrstep;
    nyrs = npts / yrlen;
    mapped_anoms = mapped_anoms_im;         % copy all intermediate anomalies into final anomalies
    outlier_mapped_anoms = zeros(npts,1);
    far_outliers = DA_mdl.far_outlier_map | new_nas;    % far outliers and anything mapped to NA
    
    nfars = sum(far_outliers);
    
    if (nfars > 0)
 
              % outlier-scaling for each day of the year
        if (~DA_hist.DP.isPrecipRun || any(strcmp(DA_hist.RP.prcp_distrib,["generalizedpareto", "inversegaussian", "loglogistic", "lognormal"])))
                                        % these precip "scalings" scale the bins rather than the precip values, so
                                        % no need to rescale the numbers, as is needed for pwr and log scaling.
                                        % NOTE:  there's no reason the pwr and log scaling couldn't be using the 
                                        % bins to make the numbers look gaussian.  Just needs to change the code
                                        % earlier.
            
            if (~DA_hist.DP.isPrecipRun)
                negs = anoms < 0;
            else  
                medians = repmat(DA_hist.pthresh_base(:,2),nyrs,1);
                negs = anoms < medians;
                far_outliers = far_outliers(~negs);
            end
                    % for temp runs, we assume the distributions have a median of 0.  May be slightly off...this should
                    % be updated?
                    % Also:  minor bug here.  Using the hist's pthresh line.  should really figure out which 
                    % yrstep we're in and use the pthresh of the yrstep...
            doy_scaling_neg = DA_obs.pthresh_base(:,1) ./ DA_hist.pthresh_base(:,1);    % for left tail corrections
            doy_scaling_pos = DA_obs.pthresh_base(:,3) ./ DA_hist.pthresh_base(:,3);    % for right tail corrections.

%                     % and replicate the scaling to be as long as the actual data
            doy_scaling_neg = repmat(doy_scaling_neg,nyrs,1);
            doy_scaling_pos = repmat(doy_scaling_pos,nyrs,1);
            
            if (~DA_hist.DP.isPrecipRun)
                outlier_mapped_anoms( negs) = anoms( negs).*doy_scaling_neg( negs);
            else
                outlier_mapped_anoms( negs) = anoms( negs);
            end
            outlier_mapped_anoms(~negs) = anoms(~negs).*doy_scaling_pos(~negs);

%               Same as above, but we're removing the medians to get anomaly from median, not anomaly from 0   

%             doy_medians_obs  = repmat(DA_obs.pthresh_base(:,2),nyrs,1);
%             doy_medians_hist = repmat(DA_hist.pthresh_base(:,2),nyrs,1);
%             doy_medians_mdl  = repmat(DA_mdl.pthresh_base(:,2),nyrs,1);
% 
%             doy_scaling_neg = (doy_medians_obs - DA_obs.pthresh_base(:,1)) ./ (doy_medians_hist - DA_hist.pthresh_base(:,1));    % for left tail corrections
%             doy_scaling_pos = (DA_obs.pthresh_base(:,3) - doy_medians_obs) ./ (DA_hist.pthresh_base(:,3) - doy_medians_hist);    % for right tail corrections.
%             
%                     % and replicate the scaling to be as long as the actual data
%             doy_scaling_neg = repmat(doy_scaling_neg,nyrs,1);
%             doy_scaling_pos = repmat(doy_scaling_pos,nyrs,1);
% 
%             outlier_mapped_anoms( negs) = (anoms( negs) - doy_medians( negs)).*doy_scaling_neg( negs) + doy_medians( negs);
%             outlier_mapped_anoms(~negs) = (anoms(~negs) - doy_medians(~negs)).*doy_scaling_pos(~negs) + doy_medians(~negs);

            % replace any far-outliers
            mapped_anoms(far_outliers) = outlier_mapped_anoms(far_outliers);
           
%       this version didn't work very well.  replaced with above code,
%       which is simpler AND more conservative.  icsf 7/5/23
%
%                     % for temp runs, we assume the distributions have a median of 0.  May be slightly off...this should
%                     % be updated?
%                     % Also:  minor bug here.  Using the hist's pthresh line.  should really figure out which 
%                     % yrstep we're in and use the pthresh of the yrstep...
% %             doy_scaling_neg = DA_obs.pthresh_base(:,1) ./ DA_hist.pthresh_base(:,1);    % for left tail corrections
% %             doy_scaling_pos = DA_obs.pthresh_base(:,3) ./ DA_hist.pthresh_base(:,3);    % for right tail corrections.
% 
%             doy_scaling = zeros(npts,2);
%                     % we can make this faster by working in sets instead of on each point individually.
%             for ix=1:npts
%                 yr = floor((ix-1)/yrlen) + 1;
%                 doy = mod(ix-1, yrlen) + 1;
%                 setnum = floor((yr-1)/pdf_yrstep) + 1;
%                 
%                 pthresh_obs = squeeze(DA_obs.pthresh_base(doy,:));
%                 try
%                     pthresh_mdl = squeeze(DA_mdl.pthresh_rolling(doy,:,setnum));
%                 catch
%                     fprintf("oops\n");
%                 end
%                 pthresh_hist= squeeze(DA_hist.pthresh_base(doy,:));
%                 obs_mean    = DA_obs.means_base(doy);
%                 hist_mean   = DA_hist.means_base(doy);
%                 mdl_mean    = DA_mdl.means_rolling(setnum, doy);
% %               doy_median_obs  = pthresh_obs(ix,2);
% %               doy_median_mdl  = pthresh_mdl(2);
%                 
%                             % goal:  scale far-outliers to equivalent distance from the obs mean.
%                             % 1st ratio:  scaling based on historical pdfs and means.
%                             % 2nd ratio:  scaling for change between historical and future model pdfs.
%                             %                     scale to historical obs/hist                       scale for model change from historical period. 
%                 doy_scaling_neg = (obs_mean - pthresh_obs(1)) ./(hist_mean - pthresh_hist(1)) .* ((pthresh_mdl(1)-mdl_mean)/(pthresh_hist(1)-hist_mean));    % for left tail corrections
%                 doy_scaling_pos = (obs_mean - pthresh_obs(3)) ./(hist_mean - pthresh_hist(3)) .* ((pthresh_mdl(3)-mdl_mean)/(pthresh_hist(3)-hist_mean));    % for right tail corrections.
%                 
%                 doy_scaling(ix, :) = [doy_scaling_neg, doy_scaling_pos];
% 
% 
%                 anom = anoms(ix);
% 
%                 if (anom < 0)
%                     outlier_mapped_anoms(ix) = (anom - mdl_mean).*doy_scaling_neg + obs_mean;
%                 else
%                     outlier_mapped_anoms(ix) = (anom - mdl_mean).*doy_scaling_pos + obs_mean;
%                 end
%                 % replace any far-outliers
%             end   
%             
%             mapped_anoms(DA_mdl.far_outlier_map) = outlier_mapped_anoms(DA_mdl.far_outlier_map);            

        else
                % doing a precip run.  
                % If prcp_distrib was not log or pwr, we don't need to do
                % anything here.
            % IAN:  same minor bug as for temperature data:  should be using rolling years' pthresh's, not hist's.
            if (DA_obs.RP.isPrecipRun && any(strcmp(DA_obs.RP.prcp_distrib,["pwr","log"])))  % if user didn't specify a distribution as one of ["generalizedpareto", "inversegaussian", "loglogistic", "lognormal"], to fit to, then
                                                                    % we need to rescale DA_hist's pthresh_base values into the scaled dataspace of the obs data first
                obs_scaling   = DA_obs.prcp_scaling.prcp_scaling;
                obs_offset    = DA_obs.prcp_scaling.prcp_offset;
                obs_prcp_min  = DA_obs.prcp_scaling.prcp_min;
                hist_scaling  = DA_hist.prcp_scaling.prcp_scaling;
                hist_offset   = DA_hist.prcp_scaling.prcp_offset;
                hist_prcp_min = DA_hist.prcp_scaling.prcp_min;

                    % for debugging
    %             aa=anoms;                           % original anoms
    %             mim=mapped_anoms_im;                % original mapped_anoms.
    %           far_ix=find(DA_mdl.far_outlier_map);% list of outliers' indexes.

                    % pthresh lines:    pthresh(1)  is the prob line for pts to left of median
                    %                   pthresh(2)  is the median line (50th percentile line)
                    %                   pthresh(3)  is the prob line for pts to right of median
                    % map the original anomalies and thresholds into obs' data space, and subtract the medians.
                        % first take model stuff out of scaled space back into linear space
                zeros_map = anoms==0;
                anoms(zeros_map)= nan;  % mask these off 
                anoms           = rescale_data_1(anoms,                     hist_scaling, hist_prcp_min, "reverse", hist_offset);
    %             mapped_anoms_im = rescale_data_1(mapped_anoms_im,            obs_scaling,  obs_prcp_min, "reverse",  obs_offset);

    %           mdl_pthresh_low = rescale_data_1(DA_hist.pthresh_base(:,1), hist_scaling, hist_prcp_min, "reverse", hist_offset);
                mdl_medians     = rescale_data_1(DA_hist.pthresh_base(:,2), hist_scaling, hist_prcp_min, "reverse", hist_offset);
                mdl_pthresh_hi  = rescale_data_1(DA_hist.pthresh_base(:,3), hist_scaling, hist_prcp_min, "reverse", hist_offset);

    %           obs_pthresh_low = rescale_data_1(DA_obs.pthresh_base(:,1), obs_scaling, obs_prcp_min, "reverse", obs_offset);
    %           obs_medians     = rescale_data_1(DA_obs.pthresh_base(:,2), obs_scaling, obs_prcp_min, "reverse", obs_offset);
                obs_pthresh_hi  = rescale_data_1(DA_obs.pthresh_base(:,3), obs_scaling, obs_prcp_min, "reverse", obs_offset);

    %           doy_scaling_neg = obs_pthresh_low ./mdl_pthresh_low;    
                doy_scaling_pos = obs_pthresh_hi  ./ mdl_pthresh_hi;

    %           doy_scaling_neg = repmat(doy_scaling_neg, nyrs,1);
                doy_scaling_pos = repmat(doy_scaling_pos, nyrs,1);

                    % don't worry about outliers below the median.  These are all within a fraction of their actual value,
                    % because they are all very small precip values.
                negs = anoms < repmat(mdl_medians, nyrs, 1);

    %           anoms( negs) = anoms( negs).* doy_scaling_neg( negs);
                pos = ~negs;
                pos_outliers_ix = find(pos & DA_mdl.outlier_map);
                pos_far_outliers_ix = find(pos & DA_mdl.far_outlier_map);
                outlier_mapped_anoms(pos_outliers_ix) = anoms(pos_outliers_ix).* doy_scaling_pos(pos_outliers_ix);

    %                 % for debugging:  show the corrections we've made.
    %             [~,jx] = sort(mapped_anoms_im(pos_outliers_ix));
    %             % orig_anoms, orig_mapped,adj_outliers, final values.  These last two columns should be identical. 
    %             z=[anoms(pos_outliers_ix),mapped_anoms_im(pos_outliers_ix),outlier_mapped_anoms(pos_outliers_ix)];
    %             pix = pos_outliers_ix(jx);
    %             [~,ia,ib] = intersect(pos_far_outliers_ix, pix);
    %             z=z(jx,:);
    %                 % re-order, low to high
    %             n=size(z,1);
    %             figure; %(98);
    %             subplot(2,1,1);
    %             plot((1:n)',z,"linewidth",2);
    %             hold on;
    %             scatter(ib, mapped_anoms_im(pos_far_outliers_ix(ia)),50, "b","filled");
    %             hold off;
    %             legend(["anom orig","mapped orig","adjusted","far outliers"],"location","northwest");    
    %             s=split(DA_obs.DA_title,":");
    %             siteLbl = s{3};
    %             title(sprintf("%s, Outliers (sorted), mm", siteLbl));
    %             ylabel("precip, mm");
    %             grid on;


                outlier_mapped_anoms = rescale_data_1(outlier_mapped_anoms, obs_scaling, obs_prcp_min, "forward", obs_offset);

    %             oma = outlier_mapped_anoms;

                    % replace any far-outliers
                mapped_anoms(pos_far_outliers_ix) = outlier_mapped_anoms(pos_far_outliers_ix);

    %                 % for debugging:  show the corrections we've made.
    % %           [~,jx] = sort(mim(pos_outliers));
    %             % orig_anoms, orig_mapped,adj_outliers, final values.  These last two columns should be identical. 
    %             z=[aa(pos_outliers_ix),mim(pos_outliers_ix),oma(pos_outliers_ix)];
    %             z=z(jx,:);                              % re-order, low to high
    %             n=size(z,1);
    %             subplot(2,1,2);
    %             plot((1:n)',z(:,1:3),"linewidth",2);
    %             hold on;
    %             scatter(ib, mim(pos_far_outliers_ix(ia)),50, "b","filled");
    %             hold off;
    %             legend(["anom orig","mapped orig","adjusted","far outliers"],"location","northwest");
    %             title("precip, scaled space");
    %             grid on;
            else
                error("error:  we shouldn't get here unless we've added a precip scaling that isn't being handled here");
            end
        end
    end
 end



function adjusted_climatology = adjust_climatology(delta_clim, moving_clim_model)
% adjusted_climatology = adjust_climatology(delta_clim_365, clim_model_moving)
%   adjusts modeled climatology so modeled hist. period matches climatology of observations
%   Adds delta_clim to moving climatology.  Delta_clim should be obs_clim - hist_clim.
%
%       Inputs:  
%           delta_clim              difference between observation climatology and modeled climatology for same period
%                                       delta = obs_clim - model_hist_clim  (365x1 vector)
%           clim_model_moving       model climatology from moving LPF  (2-D surface)
%
%       Outputs
%           adjusted_climatology    just that--adjusted climatology. (2D surface)
%   

        % get size for reshaping, and make sure we've got column matrices for the climatology.
    yrlen = length(delta_clim);
    
    nyrs = length(moving_clim_model)/yrlen;
    % 
    if (mod(nyrs, 1) ~= 0)
        throw(MException('ARRM_V2:BAD_INPUT','adjust_climatology:  length of model data not even multiple of hist. climatology''s year length'));
    end
    if (~iscolumn(delta_clim)), delta_clim = delta_clim'; end
    
    delta_clim = repmat(delta_clim, nyrs, 1);
    
    adjusted_climatology = moving_clim_model + delta_clim; % moving_clim_mdl - hist_base_clim + obs clim  ( - obs_avg_base + hist_avg_base)

end

function [mapped_anoms, mapped_nas, adjustments_gt_2] = map_anomalies(DA_mdl, DA_obs, DA_hist)     
% Maps DA_mdl's rolling_years' anomalies from DA_mdl's probability distribution to to DA_obs's probability distribution.
%   normally, DA_mdl is the disaggregated model, and DA_obs is the disaggregated observation data.
%
%   NOTE:  returned data is only for the rolling years, so these should be inserted into the correct place in DA_out,
%   not simply assigned into DA_out.

% for non-rolling analysis, simply use the entire rolling years as 1 rolling step.
%
% Returns the mapped anomalies.
%   also returns an a map  of any new na's created in the mapping process.
%   adjustment_fails is a count of the number of distribution points which moved by a factor of 2 or more from the
%   historical model period to the rolling future model period.
%   These should all be near zero. If too many of those are showing up, adjust the parameter for the 
%   sigmoid used to limit the change for points near the zero point.  
%   If it is occurring away from zero, then there's a problem somewhere... this implies that the distributions are
%   drastically different between the historical and future periods;  likely to happen only if there is a discontinuity
%   in one of the distribution surfaces.

    yrlen      = DA_mdl.DP.yrlen;
    if (DA_obs.DP.isPrecipRun && any(strcmp(DA_obs.RP.prcp_distrib,["pwr","log"])))
        obs_bins   = DA_obs.prcp_scaling.bins;
        hist_bins  = DA_mdl.prcp_scaling.bins;
    else
        obs_bins   = DA_obs.RP.out_bins;
        hist_bins  = DA_hist.RP.out_bins;
    end
    pdf_yrstep = DA_mdl.RP.pdf_yrstep;
%   cdf_thresh = DA_mdl.RP.cdf_thresh;    
%   mdl_bins   = DA_mdl.RP.bins;
%   bin_delta  = 0; % bins(2)-bins(1);
    
    
        % extract the data for the rolling years only.
        % we can only map the rolling_yrs anomalies because the rolling year PDFs only use the rolling years' data...
   [rstart,rend, nyrs] = DA_mdl.using_range('rolling_yrs');
    mdl_anoms      = reshape(DA_mdl.anoms(rstart:rend),              yrlen, nyrs);
    mdl_anom_probs = reshape(DA_mdl.anom_probs_rolling(rstart:rend), yrlen, nyrs);
    na_map_in = reshape(DA_mdl.na_map_in(rstart:rend), yrlen, nyrs);
    if (DA_mdl.DP.isPrecipRun)
        zeros_map = reshape(DA_mdl.zeros_map(rstart:rend), yrlen, nyrs);
    end
    CDF_model = DA_mdl.CDF_rolling;
%   pdf_model = DA_mdl.pdf_rolling;
    nsets = size(CDF_model,3); 
    
    mapped_anoms   = zeros(yrlen, nyrs);
    hist_anoms     = zeros(yrlen, nyrs);
    
        % get flags of which points in the CDF to use.
        % Matlab's interp1 gets upset if the x-values aren't strictly monotonically increasing, so we need to use a
        % subset of the CDF.
%     obsflags  = 	~isnan(DA_obs.CDF_base)  & (DA_obs.CDF_base  >= DA_obs.RP.cdf_thresh)  & (DA_obs.CDF_base  <= 1-DA_obs.RP.cdf_thresh)  & DA_obs.pdf_base  > 0;
%     histflags = 	~isnan(DA_hist.CDF_base) & (DA_hist.CDF_base >= DA_hist.RP.cdf_thresh) & (DA_hist.CDF_base <= 1-DA_hist.RP.cdf_thresh) & DA_hist.pdf_base > 0;
% 
    obsflags  = 	~isnan(DA_obs.CDF_base)   & DA_obs.pdf_base  > DA_obs.RP.na_thresh;
    histflags = 	~isnan(DA_hist.CDF_base) &  DA_hist.pdf_base > DA_hist.RP.na_thresh;


%       Step 1:  get the probability value for each day of model data in rolling period.
%           we actually did this during the disaggregation process, so no longer done here.
         
            % now map the model anomalies via CDF_model probabilities to DA_obs.CDF_base probabilities.
%       Step 2:  get obs. value for each pr_mdl's probability.
       
    for i=1:yrlen   % for each day of year...
        oflags = obsflags(:,i);
        hflags = histflags(:,i);
        for j = 1:nsets     % for each rolling probability set
            try
%               mflags = squeeze(mdlflags(i,:,j)); %  ???
                ix1 = (j-1)*pdf_yrstep + 1;  % year start index
                if (j==nsets)
                    ix2 = nyrs;               % year end index
                else
                    ix2 = j*pdf_yrstep;
                end
                            % get observation value for each model anomaly's probability
                mapped_anoms(i,ix1:ix2) = interp1(DA_obs.CDF_base(oflags, i),       obs_bins(oflags),  mdl_anom_probs(i,ix1:ix2), 'makima',nan);
                            % get historical model CDF's value for each model anomaly's probability
                hist_anoms(i,ix1:ix2)   = interp1(DA_hist.CDF_base(hflags, i),      hist_bins(hflags), mdl_anom_probs(i,ix1:ix2), 'makima',nan);
                            % (no need to get model's value for it's probability...we have it already!
%               mdl_vals  = interp1(DA_mdl.CDF_rolling(mflags,i,j), bins(mflags), mdl_anom_probs(i,ix1:ix2), 'pchip',nan);         
            catch me
                DA_mdl.DP.warn_log( '--------%s: oops: caught exception----------\n', mfilename);     % shouldn't happen...for debugging!
                msgtext=getReport(me);
                DA_mdl.DP.warn_log( '%s\n', msgtext);
                DA_mdl.DP.warn_log( '------\n');
            end 

            
                        % for analysis...show where we mapped to NAs.
                        % we'll report this properly in check_nas.
            if (~DA_mdl.DP.isPrecipRun)
                newfails = isnan(mapped_anoms(i,ix1:ix2)) & ~(na_map_in(i, ix1:ix2)); 
            else
                newfails = isnan(mapped_anoms(i,ix1:ix2)) & ~(na_map_in(i, ix1:ix2)) & ~zeros_map(i, ix1:ix2); 
            end

            if (any(newfails(:)))               
                yix = find(newfails);
                doy = i;
                for k=1:length(yix)
                    yy = ix1+yix(k)-1;
                    actual_yr = DA_mdl.DP.rolling_yrs(1) + yy - 1;
                    prb = min( mdl_anom_probs(doy,yy), 1- mdl_anom_probs(doy,yy));
                    nanix = (yy-1)*yrlen+doy;
                    DA_mdl.DP.print_log( 'mapped nan:          %20s year: %4d doy: %3d (ix %5d yr_off %3d day %4d)  anom %8.5f mapped to %8.5f prob %7.4e\n', DA_mdl.DP.runLbl, actual_yr, doy, nanix, yy, doy, mdl_anoms(doy,yy), mapped_anoms(doy,yy), prb);
                end
            end
        end
    end
    
                % now, if we're using rolling probabilities, see how much the value changed from historical to rolling
                % future probability, and adjust the mapped anomalies by the same factor.
    if (~DA_mdl.DP.isPrecipRun || any(strcmp(DA_mdl.RP.prcp_distrib,["generalizedpareto", "inversegaussian", "loglogistic", "lognormal"])))
% %       adjustment = mdl_anoms ./ hist_anoms;   % scaling correction for change in probability distributions over time
                                                % this unfortunately blows up as hist_anoms goes to zero, so we scale
                                                % it with a sigmoid so scaling goes to 1 as hist_anoms goes to 0.
        adjustment = 1 + (mdl_anoms - hist_anoms)./hist_anoms.*jc_sigmoid(hist_anoms, 2.5,'abs');
        adjustments_gt_2 = sum(abs(adjustment(:)) > 2 );
        mapped_anoms = mapped_anoms .* adjustment;    
    else
%         adjustments_gt_2 = 0;
%         prcp_min    = prcp_scale_info.prcp_min;
%         scaling     = prcp_scale_info.prcp_scaling;
%         prcp_offset   = prcp_scale_info.prcp_offset;
%         scaled_prcp = rescale_data_1(prcp, scaling, prcp_min, "forward", prcp_offset);
        
        prcp_min = DA_obs.prcp_scaling.prcp_min;
        
        sigmoids = prcp_sigmoid(hist_anoms, DA_hist.pthresh_base);
        
            % For precip, we want to caculate the shift back in the regular precip space, rather than the scaled space.
            % mdl & hist *should* be using the same scaling, but just to be cautious in case someone changes that down
            % the road, am using each DA's specific scaling.
        mdl_anoms  = rescale_data_1( mdl_anoms, DA_mdl.prcp_scaling.prcp_scaling,  DA_mdl.prcp_scaling.prcp_min, "reverse",  DA_mdl.prcp_scaling.prcp_offset);  
        hist_anoms = rescale_data_1(hist_anoms,DA_hist.prcp_scaling.prcp_scaling, DA_hist.prcp_scaling.prcp_min, "reverse", DA_hist.prcp_scaling.prcp_offset);
        
        
        adjustment = 1 + (mdl_anoms - hist_anoms)./hist_anoms.*sigmoids;   % shouldn't need "abs", since all values should be positive.
        adjustments_gt_2 = sum(abs(adjustment(:)) > 2 );
%       mdl_anoms  = rescale_data_1( mdl_anoms, DA_mdl.prcp_scaling.prcp_scaling,  DA_mdl.prcp_scaling.prcp_min, "forward",  DA_mdl.prcp_scaling.prcp_offset);  
        
        mapped_anoms = rescale_data_1( mapped_anoms, DA_obs.prcp_scaling.prcp_scaling,  DA_obs.prcp_scaling.prcp_min, "reverse",  DA_obs.prcp_scaling.prcp_offset);
        mapped_anoms = max(prcp_min, mapped_anoms .* adjustment);    
        mapped_anoms = rescale_data_1( mapped_anoms, DA_obs.prcp_scaling.prcp_scaling,  DA_obs.prcp_scaling.prcp_min, "forward",  DA_obs.prcp_scaling.prcp_offset);

% figure;
% xx=repmat(1:365,1,nyrs)';
% yy=reshape(repmat(1:nyrs,365,1),[],1);
% scatter3(xx,yy,adjustment(:),10,"filled")  
% grid on;
% title(extractAfter(DA_mdl.DA_title,7),"interpreter","none");
% view([0,00]);
%         adjustment = 1;
%         adjustments_gt_2 = 0;
    end

%                     % for precip, if mapped any small precips below
%                     % prcp_min, or to NA, adjust it to prcp_min. 
%             if (isPrecipRun)
%                 negs = mdl_anom_probs(i, ix1:ix2) < .5;
%                 

    
        
% %         manoms = mapped_anoms;
    mapped_anoms = reshape(mapped_anoms,nyrs*yrlen,1);

    mapped_nas = isnan(mapped_anoms) & ~na_map_in(:);
        
% %         deltas = mapped_anoms - manoms;       % code to plot the adjustments over entire year and over each season.
% %         del1 = deltas(60:150,:);              % this shows whether the shift-over-time varies by season, and also
% %         del2 = deltas(151:242,:);             % shows if something looks wierd about the adjustments.
% %         del3 = deltas(253:334,:);
% %         del4 = deltas([1:59,335:365],:);        % Add color markings for any adjustment_fails, Ian!
% %         man1 = manoms(60:150,:);
% %         man2 = manoms(151:242,:);
% %         man3 = manoms(253:334,:);
% %         man4 = manoms([1:59,335:365],:);
% %         figure(21); scatter(manoms(:),deltas(:),20,'filled'); grid on        
% %         figure(22); subplot(2,2,1); scatter(man1(:),del1(:),10,'filled'); grid on;  xlabel('mapped anomaly'); ylabel('change over time'); title('spring'); 
% %         figure(22); subplot(2,2,2); scatter(man2(:),del2(:),10,'filled'); grid on;  xlabel('mapped anomaly'); ylabel('change over time'); title('summer');
% %         figure(22); subplot(2,2,3); scatter(man3(:),del3(:),10,'filled'); grid on;  xlabel('mapped anomaly'); ylabel('change over time'); title('fall');      
% %         figure(22); subplot(2,2,4); scatter(man4(:),del4(:),10,'filled'); grid on;  xlabel('mapped anomaly'); ylabel('change over time'); title('winter');
% %         grid on;
% % %       figure(3); scatter(manoms(:),deltas1(:),20,'filled'); grid on        
% % %       figure(4); scatter(manoms(:),deltas2(:),20,'filled'); grid on  
% % %       fprintf('adjustment_fails:  0: %3d    -: %3d    1: %3d    2: %3d\n', adjustment_fails0, adjustment_fails, adjustment_fails1, adjustment_fails2); 

end

function sigmoids = prcp_sigmoid(anoms, pthresh)

    % med = 6;
    % figure(99);  plot(x,normcdf(12/med*(x-med/2))); xlim([0,med+1]); grid on;
    %   This generates a scaling value for the correction for each anom point.
    %   the formula above (used below), produces a sigmoid from 0-1, centered at med which follows a normal CDF starting
    %   at about 1/4 med and going to 3/4 med.  
    %   We'll use this to apply the full correction for points above 3/4(median - min);
    %   For points below 1/4(median - min), we'll apply no correction to the mapped point.
    %   for values inbetween, we'll apply a fraction of the correction, as given by the formula above.
    %   Median is given by pthresh(:,2);
    %   min    is given by pthresh(:,1);  and is the point where the probability goes to about 1e-8.
    %
    %   Returns the correction specific to each anom point.

    yrlen = size(anoms,1);
    medians = pthresh(:,2);
    mins    = pthresh(:,1);
    deltas  = medians - mins;
    mids = deltas/2;
    
    aa = anoms - mins;
    
    sigmoids = zeros(size(anoms));
    for i=1:yrlen
        sigmoids(i,:) = normcdf(12/mids(i)*(aa(i,:)-mids(i)/2));
    end
end
% function flags = mk_flags(CDF, pdf, thresh)
% % Flags the data points in the CDF to use for interp1 so x-values are monotonically increasing
% %
% %     yrlen = size(CDF,2);
% %     fstart = zeros(1,yrlen);
% %     fend   = zeros(1,yrlen);
%     
%             % mask off everything outside thresh 
% 	CDFflags = ~isnan(CDF) & (CDF >= thresh) & (CDF <= 1-thresh);
%     pdfflags = ~isnan(pdf) &  pdf >= thresh; 
%     flags = CDFflags & pdfflags;
%     
% %             % find index of 1st and last points inside the NA threshold
% %     for i=1:yrlen
% %         fstart(i) = find(flags(:,i),1);
% %         fend(i)   = find(flags(:,i),1,'last');
% %     end
% end
% 
function [DA_obs, DA_hist, DA_mdl, RP, DP, ttl, do_histograms] = init_params(DA_obs, DA_hist, DA_mdl, varargin)
% parse inputs, and verify that data and run parameters are all consistent.
% Returns RP & DP of DA_mdl, so that DA_out will be same size, shape, etc. as DA_mdl.  
% We'll limit DA_out to rolling_yrs only right at the end.
% throws an exception if anything looks wrong.
%
%   DA_checks  = [];
    DP_checks  = ["runType","yrlen"];
 
%   RP_checks  = ["runType","yrlen","bins","cdf_thresh","clim_nterms","clim_sig_terms","anom_nterms","anom_sig_terms", "sigrange"];
%   RP_checks  = ["runType","yrlen",       "cdf_thresh","clim_nterms","clim_sig_terms","anom_nterms","anom_sig_terms", "sigrange"];
    RP_checks  = ["runType","yrlen",       "cdf_thresh",                                                                         ]; %"sigrange"];
    
    p = inputParser;
    p.StructExpand = true;
    p.KeepUnmatched = true;    
    addParameter(p,'title',[]);
    addParameter(p,'RP',[]);
    addParameter(p,'DP',[]);
    addParameter(p,'do_histograms',true);
    parse(p, varargin{:});
    
    if (~isempty(fieldnames(p.Unmatched))), error("ARRM_V2_error:  unmatched input parameters!"); end
    
    ttl = p.Results.title;
    RP = p.Results.RP;
    DP = p.Results.DP;                    
    
    do_histograms = p.Results.do_histograms;
    if (length(do_histograms) == 1), do_histograms = repmat(do_histograms,1,3); end
    
        % check the DownScalingParams and RunParams objects to make sure they match between datasets.
    for i=1:length(DP_checks)
        check_param(DP_checks(i), DA_obs.DP, DA_hist.DP, DA_mdl.DP)
    end
    
    for i=1:length(RP_checks)
        check_param(RP_checks(i), DA_obs.RP, DA_hist.RP, DA_mdl.RP)
    end    
        % we'll assume that the RPs and DPs are identical across all inputs.  
        % Maybe need to code a check for that here, Ian.
    if (isempty(RP)), RP = DA_mdl.RP; end
    if (isempty(DP)), DP = DA_mdl.DP; end
    
    if (RP.yrlen ~= DP.yrlen), error("ARRM_V2_error:  yrlen's do not match across RP & DP"); end  
end

function check_param(param, s1, s2, s3)     % makes sure that the run params match between obs, hist & model

    for i=1:numel(s1.(param))
        if (s1.(param)(i) ~= s2.(param)(i))
            error("ARRM_V2_error:  %s: input parameters don't match", param);
        end
        if (s1.(param)(i) ~= s3.(param)(i))
            error("ARRM_V2_error:  %s: input parameters don't match", param);
        end
    end
end

function DA_histogram(DA, do_histograms, lbl, pdf_yrs) %#ok<INUSD>      % Ian:  check that we don't need pdf_yrs here.
% does histogramming on DA (obs, hist or model).

% 1.  check that DA is disaggregated.

    if (isempty(DA.anoms))
        DA.calc_anomalies();            
        if (DA.insufficient_data), error("insufficient data: %s", lbl); end            
    end
            
% 2.  check that all sets have the same binning

    do_binning = isempty(DA.RP.edges);
    
    if (DA.DP.isPrecipRun)
            % for precip, get the scaling from the DA's scaling info.
        DA.RP.do_calc_binning = false;
        if (any(strcmp(DA.RP.prcp_distrib,["pwr","log"])))
            DA.RP.edges = DA.prcp_scaling.edges;
        end
        
    else
        if (do_binning)
            do_histograms = true;
%           [edges,         ~,      ~, dx] = DA.calc_binning();  %       [edges, edgerange, nedges, dx]
%           [edges, edgerange, nedges, dx] = DA.calc_binning();
            [edges,         ~,      ~, dx] = DA.calc_binning();
            minedge = edges(1);
            maxedge = edges(end);

            nbins = ceil((maxedge - minedge)/dx);
            nedges = nbins+1;
            edges = linspace(minedge,maxedge, nedges);
                    % now find max bin range and number of bins, so all sets use the same binning.

            DA.RP.do_calc_binning = false;     % make sure we don't redo the binning
            DA.RP.edges = edges;
        end
    end                            
% 3.  calc histograms if needed, and clear out the pdfs so we calculate them later.

    if (do_histograms)
        DA.calc_hists();
        DA.pdf_base                = [];
        DA.CDF_base                = [];
        DA.CDF_okflags_base        = false(0,0);
        DA.pdf_rolling             = [];
        DA.CDF_rolling             = [];
        DA.CDF_okflags_rolling     = false(0,0);
    end
end

