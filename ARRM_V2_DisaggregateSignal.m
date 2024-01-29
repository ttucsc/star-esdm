classdef ARRM_V2_DisaggregateSignal < handle
    %Object to hold disaggregated climate data
    %   Detailed explanation goes here
        
    properties
        DA_title        = "empty";
        insufficient_data = false;      % flag set to true if insufficient data available to process.
        anoms           = [];
        na_map_in       = false(0,0);       % flags NAs on input
        zeros_map       = [];               % for precip.  keep track of zeros on input.
        outlier_map     = false(0,0);       % flags data with probability beyond RP.outlier_thresh.
        far_outlier_map = false(0,0);       % flags data with probability beyond RP.far_outlier_thresh
        na_map_mapped   = false(0,0);       % flags data mapped to NAs by ARRM_V2 downscaling process.
        avg_base        = [];       % average of base years only.
        grid_pt_used    = [];
        trend           = [];
        base_clim       = [];
        base_hist       = [];
        moving_clim     = [];       % 2-D climatology, unwrapped to single time series
        moving_clim_1D  = [];       % like moving_clim, but not smoothed along day-of-year.
        phase_clim      = [];       % phase term (in days) relative to pure cosine for base climatology. 
        phase_moving    = [];       % phase term (in days) relative to pure cosine for moving climatology.
        trend_params    = [];       % struct with trend parameters for entire set.
        trend_params_daily = [];    % array of polyval parameters for trend calculated for every day of the year
        rolling_clim    = [];       % climatology of each rolling_set (usually 31 years of data).
        rolling_hist    = [];
        lohi_base       = [];
        lohi_all        = [];

        RP              = [];       % ARRM_V2_RunParams with run parameters used to disaggregate the data.
        DP              = [];       % ARRM_V2_DataParams with data parameters used to disaggregate the data.
        
        sigma_normalized= false;
        norm_sigmas     = [];           % sigmas used for normalizing.

                                    % some precip-specific stuff for debugging:
        wet_day_clim    = [];           % wet-day  climatology.  Will be adjusted after trimming off drizzle
        wet_day_clim_raw= [];           % unfiltered climatology.
        wet_day_clim_0  = [];           % for saving the pre-trimming climatology               drizzle-trimming is done in  pre-processing for ARRM_V2
        wet_day_rolling_clim     = [];
        wet_day_rolling_clim_raw = [];
        wet_day_phase   = [];
        wet_day_avg     = [];
        excess_thresh   = [];           % thresholds used to trim precip.                       this is also done in  pre-processing for ARRM_V2
        
                                    % struct for scaling info when scaling/unscaling precip data.
        prcp_scaling    = struct("scaling",1, "dx",0,"sigma",0,"edges",[], "pr_min", 0);
        
        raw_data        = [];       % will contain raw data on output.
        raw_fix         = [];       % raw data, with NAs replaced by low-pass-filtered value.
        
                                        % the following 4 are used for downscaling, otherwise empty
        orig_anoms              = [];   % original anomalies;  only set in output DAs.
        mapped_anoms_im         = [];   % intermediate mapped (downscaled) anomalies, with no corrections to anomalies                                
        outlier_mapped_anoms    = [];   % outlier-corrected values (for all points, even if not outliers)
        mapped_anoms            = [];   % final mapped anomalies, with far-outliers replaced by outlier-mapped values.
        mapped_output           = [];   % final mapped (downscaled) output.  Mapped_anoms + 
                                    
        
        % the following are only calculated if RP.do_pdfs is  true.
            % probability surfaces.  
        pdf_base                = [];
        CDF_base                = [];
        CDF_okflags_base        = false(0,0);
        pdf_rolling             = [];
        CDF_rolling             = [];
        CDF_okflags_rolling     = false(0,0);
        
        % the following are only calculated if DP.do_base_probs or DP.do_rolling_probs is true
        
        anom_probs_base         = [];
        anom_probs_rolling      = [];
        anom_sdevs_base         = [];       % std deviation equivalents of each anomaly : norminv(anom_probs_base | anom_probs_rolling)
        anom_sdevs_rolling      = [];
        
            % useful stats      
            % NOTE:  These stats are calculated from the smoothed output pdf's, not on the raw histograms.
        means_base              = [];
        sigmas_base             = [];
        skewness_base           = [];
        kurtosis_base           = [];
        pthresh_base            = [];       % outlier threshold line (365 x 2, for low & high side)
        
        means_rolling           = [];
        sigmas_rolling          = [];
        skewness_rolling        = [];
        kurtosis_rolling        = [];
        pthresh_rolling         = [];
        

        % the following are only calculated if RP.do_pdfs is true and RP.do_problines is true.
            % probability Isolines.
        problines_base      = [];
        pdf_xlines_base     = [];
        CDF_xlines_base     = [];
        probcounts_base     = [];
        totcount_base       = [];
        
        problines_rolling   = [];
        pdf_xlines_rolling  = [];
        CDF_xlines_rolling  = [];
        probcounts_rolling  = [];
        totcounts_rolling   = [];

            % the following are used to hold the result from the entire data read for plotting and visualizing.
            % they are filled only when a call to DA.trim_data_yrs( ) is made.  Otherwise are left blank.
        anoms_all           = [];       % over entire data read
        trend_all           = [];       % over entire data read
        moving_clim_all     = [];       % over entire data read
        moving_clim_1D_all  = [];       % over entire data read
        
        na_map_unknown      = true;
        binning_checked     = false;
        Userdata            = struct;   % so users can tack on additional info as needed
        
    end
        
    properties (Dependent)
        DAType;                 % from DA.DP.dspType:                obs,       hist,       model,         downscaling or   global
        DA_yrs;                 % data years for the DA's DAType:    obs_yrs,   base_yrs,   rolling_yrs,   rolling_yrs,     data_yrs
        DA_yrlbl                % returns field name for DA_yrs:    'obs_yrs', 'base_yrs', 'rolling_yrs', 'rolling_yrs',   'data_yrs'
    end
    
    properties(Constant)
        DA_types = {'obs',     'hist',     'model',     'global',  'downscaling', '',         'unknown'};
        yr_types = {'base_yrs','base_yrs','rolling_yrs','data_yrs','rolling_yrs', 'data_yrs', 'data_yrs'};
    end
    
    methods
        
        function DA = clone(obj, RP, DP)                            % clones a DA 
            if (~exist('RP','var')), RP=obj.RP; end
            if (~exist('DP','var')), DP=obj.DP; end
            DA = ARRM_V2_DisaggregateSignal([],RP,DP);
            mc=?ARRM_V2_DisaggregateSignal;  % get metadata about DA class.
            props = properties(obj);
            for i=1:length(props)
                if (~mc.PropertyList(i).NonCopyable)
                    DA.(props{i}) = obj.(props{i});
                end
            end
            if (nargin >= 2)
                DA.RP = RP;
            end
            if (nargin >= 3)
                if (~isempty(DA.RP))
                    DP.rolling_steps = DP.calc_rolling_steps(DA.RP.pdf_yrstep);
                end
                DA.DP = DP;
            end
            return;
        end
        
        function obj = ARRM_V2_DisaggregateSignal(raw_data, RP, DP, ttl, do_anomalies, do_histograms, do_pdfs, do_base_probs, do_rolling_probs, do_problines)

            if ( exist('RP', 'var') && ~isempty(RP)),  obj.RP       = RP;  end
            if ( exist('DP', 'var') && ~isempty(DP)),  obj.DP       = DP;  end
            if ( exist('ttl','var') && ~isempty(ttl)), obj.DA_title = ttl; end
            
                % put RunParams & DataParams into output
            obj.RP = RP;
            obj.DP = DP;
            
                    % if just passed a size, create working arrays and return
            if (numel(raw_data) == 1)
                npts = raw_data;
                raw_data            = [];
                obj.raw_data        = [];
                obj.anoms           = nan(npts,1);
                obj.mapped_output   = nan(npts,1);
                obj.na_map_in       = false(npts,1);
                obj.zeros_map       = [];
                obj.outlier_map     = false(npts,1);
                obj.far_outlier_map = false(npts,1);
                obj.na_map_mapped   = false(npts,1);
                obj.grid_pt_used    = zeros(npts,1,"uint8");
            end
            
            if (isempty(raw_data))      % creating an empty one?
                return;                 % if so, just return now.
            end
            
            if (~exist('do_anomalies',    'var') || isempty(do_anomalies    )), do_anomalies     = RP.do_anomalies;     end
            if (~exist('do_histograms',   'var') || isempty(do_histograms   )), do_histograms    = RP.do_histograms;    end
            if (~exist('do_pdfs',         'var') || isempty(do_pdfs         )), do_pdfs          = RP.do_pdfs;          end
            if (~exist('do_base_probs',   'var') || isempty(do_base_probs   )), do_base_probs    = RP.do_base_probs;    end
            if (~exist('do_rolling_probs','var') || isempty(do_rolling_probs)), do_rolling_probs = RP.do_rolling_probs; end
            if (~exist('do_problines',    'var') || isempty(do_problines    )), do_problines     = obj.RP.do_problines; end
            
            obj.raw_data = double(raw_data);    % make sure it's double.
                                        % set up flag arrays for NAs and outliers.
            if (obj.DP.isPrecipRun)                 
                obj.zeros_map   = raw_data <= 0;    % dry days.  Mapping process may map some zeros to tiny values (<1e12), so we'll mask them back out at the end.
            else
                obj.zeros_map   = false(0,0);
            end
            obj.na_map_in       = isnan(raw_data);
            obj.na_map_unknown  = false;
            obj.outlier_map     = false(size(raw_data));
            obj.far_outlier_map = false(size(raw_data));
            obj.na_map_mapped   = false(size(raw_data));
            
            raw_yrs = length(raw_data)/obj.DP.yrlen;
            if (raw_yrs ~= range(obj.DP.data_yrs)+1)
                obj.DP.error_log("error:  input raw data doesn't cover entire data period:  raw_yrs:  %.2f  Data_yrs:  %d - %d", raw_yrs, obj.DP.data_yrs(1), obj.DP.data_yrs(2));
            end
            
            if (~isempty(obj.DP.base_yrs))
                ok = obj.check_nans('base_yrs');
                if (~ok)
                    obj.insufficient_data = true;
                    return;
                end
            end
            
            if (~isempty(obj.DP.rolling_yrs))
                ok = obj.check_nans('rolling_yrs');
                if (~ok)
                    obj.insufficient_data = true;
                    return;
                end
            end
            
            if (~any([do_anomalies, do_histograms, do_pdfs, do_base_probs, do_rolling_probs, do_problines]))
                return
            end
            obj.run_disaggregation(do_anomalies, do_histograms, do_pdfs, do_base_probs, do_rolling_probs, do_problines)
            
        end
        
        function run_disaggregation(obj, do_anomalies, do_histograms, do_pdfs, do_base_probs, do_rolling_probs, do_problines)
            
            if (~exist('do_anomalies',    'var') || isempty(do_anomalies    )), do_anomalies     = obj.RP.do_anomalies;     end
            if (~exist('do_histograms',   'var') || isempty(do_histograms   )), do_histograms    = obj.RP.do_histograms;    end
            if (~exist('do_pdfs',         'var') || isempty(do_pdfs         )), do_pdfs          = obj.RP.do_pdfs;          end
            if (~exist('do_base_probs',   'var') || isempty(do_base_probs   )), do_base_probs    = obj.RP.do_base_probs;    end
            if (~exist('do_rolling_probs','var') || isempty(do_rolling_probs)), do_rolling_probs = obj.RP.do_rolling_probs; end
            if (~exist('do_problines',    'var') || isempty(do_problines    )), do_problines     = obj.RP.do_problines;     end
            
            if (do_anomalies)
                obj.calc_anomalies();       % anoms zero-mean (base period).  Also calculated climatology, moving climatology, avg_base, & trend,
                                            % all of which are stored in the object's data members.
                if (obj.insufficient_data), return; end
            end
            if (obj.DP.isPrecipRun && (~do_histograms && ~do_pdfs && ~do_base_probs && ~do_rolling_probs && ~do_problines)), return; end
            if (isempty(obj.orig_anoms))
                obj.orig_anoms = obj.anoms;
            end
            
            obj.calc_binning();
            if (obj.insufficient_data), return; end
            
            if (do_histograms)
                obj = obj.calc_hists(); 
                if (obj.insufficient_data), return; end
            end
            
            if (do_pdfs)                
                obj.calc_pdfs(obj.RP.do_pdfs);
                if (obj.insufficient_data), return; end
            end
            
            if (do_problines)
                obj.calc_problines();
            end
            
            if (do_base_probs)
                obj.calc_base_probs();
                if (obj.insufficient_data), return; end
            end                

            if (do_rolling_probs)
                obj.calc_rolling_probs();
                if (obj.insufficient_data), return; end
            end
     
        end
        
        function obj = calc_anomalies(obj)
        % function to calculate climatologies, trend and anomalies from the raw data.
        
            if (obj.na_map_unknown)
                obj.na_map_in = isnan(obj.raw_data);
                obj.na_map_unknown = false;
            end
            
            npts = numel(obj.raw_data);
            if (npts==0), obj.insufficient_data = true; return; end
            
                % abort run if user hasn't set sigma_normalize and median_normalize.            
                %   Both can be false, but must be explicitly set to false if neither is set to true.
                %   Setting either one will turn off the other, as they are mutually exclusive.  
            if (~isempty(obj.RP.sigma_normalize))
               if (obj.RP.sigma_normalize), obj.RP.median_normalize = false; end
%                if (obj.RP.sigma_normalize), obj.RP.median_normalize = true; end  % This is probably more correct...but needs some creative thinking.
            else
                if (obj.DP.isPrecipRun)             % added icsf 7/24/23, to default to true if temp run, false if precip run.  BUT:  there's probably a better place to set this default, Ian!
                    obj.RP.sigma_normalize = false;
                else
                    obj.RP.sigma_normalize = true;
                end
                if (~isempty(obj.RP.median_normalize))
                    if (obj.RP.median_normalize), obj.RP.sigma_normalize = false; end
                end
            end

                    % commented out icsf 7/24/23.  see above, where defaults are now set.
%             if (isempty(obj.RP.sigma_normalize))
%                 obj.DP.error_log("error:  RunParams' sigma_normalize is empty.  It must be explicitly set");
%             end
%             if (isempty(obj.RP.median_normalize))
%                 obj.DP.error_log("error:  RunParams' median_normalize is empty.  It must be explicitly set");
%             end                                
            
            nyrs = npts/obj.RP.yrlen;
            if (mod(nyrs,1) ~= 0), obj.DP.error_log("ARRM_V2_DisaggregateSignal.calc_clims_anoms:  length(raw_data) is not multiple of yrlen"); end
%             if (nyrs ~= obj.RP.nyrs), error("ARRM_V2_DisaggregateSignal:  data length (%d yrs) doesn't match RunParams nyrs (%d)",nyrs,obj.RP.nyrs); end

                % fix NAs and get NA maps.
            if (obj.DP.isPrecipRun)
                obj.raw_fix = obj.raw_data;
                obj.raw_fix(obj.raw_fix < obj.DP.prcp_min) = 0;     % "Raw Fix" for precip simply assumes NAs are dry days.
                obj.raw_fix(isnan(obj.raw_fix)) = 0;        % set NAs to zero for climatology calculations.  We'll reset them back to NAs down below.
                obj.anoms   = obj.raw_data;                 % all dry days will later be set to NAs, so we don't set them to zeros here.
                                                            % we might want to make use of the trace days later, so we
                                                            % leave them in the anoms for now.
                obj.orig_anoms = obj.anoms;
                obj.calc_precip_climatology("raw_data");    % gets climatology of the precip counts, for histogram equalization later.
                                                            % we're using raw_data to calculate the climatology so it
                                                            % can adjust for NAs to avoid sampling bias (such as 
                                                            % skipping most Sundays and holidays) in observation data.  
                    % skip getting climatologies for precip (other than climatology of #wet-days)
  %             return;     % skip
            elseif (obj.RP.do_fixNAs)
                obj.raw_fix   = obj.fix_nans(obj.raw_data, obj.RP.yrlen, nyrs);
            else
                obj.raw_fix = obj.raw_data;
            end

                % For precip, the following code calculates the average
                % precip on each day of the year, same as for temp data.
                % But we use the wet-day climatology (calc'd above) for
                % histogram-equalization.
            
                        % get basic climatology (straight average, smoothed) for base years
                        % straight-average climatology, smoothed by low-pass filter.
                        % NOTE:  using data with NA's replaced with reasonable values, so there is no bias if there are
                        % more missing days at certain times of the year.

            if (~isempty(obj.DP.base_yrs))
                [bstart, bend] = obj.using_range('base_yrs');
                [obj.base_clim, obj.avg_base, obj.phase_clim]  = obj.DA_climatology(obj.raw_fix(bstart:bend), obj.RP.clim_nterms, obj.RP.clim_sig_terms, obj.RP.yrlen);         
            else
                [rstart, rend] = obj.using_range('rolling_yrs');
                [obj.base_clim, obj.avg_base, obj.phase_clim]  = obj.DA_climatology(obj.raw_fix(rstart:rend), obj.RP.clim_nterms, obj.RP.clim_sig_terms, obj.RP.yrlen);       % must be doing model w/ separate hist.
                                                                                                                                % So for model, use rolling years 
                                                                                                                                % to get base climatology and average.
                
            end
                            
                        % climatology for each rolling step
            if (~isempty(obj.DP.rolling_yrs))
                nsets = length(obj.DP.rolling_steps);
                obj.rolling_clim = zeros(nsets,obj.RP.yrlen);
                for i=1:nsets
                    [rstart, rend] = obj.using_range('rolling_yrs', i);
                    try
                        obj.rolling_clim(i,:) = obj.DA_climatology(obj.raw_fix(rstart:rend), obj.RP.clim_nterms, obj.RP.clim_sig_terms, obj.RP.yrlen); 
                    catch
                        fprintf("oops rolling clim\n");
                    end
                end
            end                        
                        
                        % NOTE: anomalies and climatology include ALL data read from file, which may be longer than
                        %       either the base_yrs or the rolling_yrs.  That's so we use as much real data as possible
                        %       when filtering surfaces and signals.

            if (~obj.DP.isPrecipRun && ~isempty(obj.DP.trend_yrs))
                mid_data = obj.detrend_valid_only(); 
                if (obj.insufficient_data), return; end       % insufficient data to process.
                [obj.moving_clim, obj.moving_clim_1D, obj.phase_moving] = obj.climatology_3D(mid_data);   % moving clim is zero-mean.
                obj.anoms   = obj.raw_fix - obj.avg_base - obj.trend - obj.moving_clim;
                if (obj.RP.do_daily_trends)
                    obj.trend_params_daily = obj.calc_daily_trend_params();  % calc trend params for each day of the year.  
                end
            else
                [obj.moving_clim, obj.moving_clim_1D, obj.phase_moving] = obj.climatology_3D(obj.raw_fix);
                if (~obj.DP.isPrecipRun)
                    obj.anoms   = obj.raw_fix;      % don't muck with the anoms if working with precip!
                end
            end            
            
            obj.orig_anoms = obj.anoms;
            if (any(obj.na_map_in))
                obj.orig_anoms(obj.na_map_in) = nan;        % need to make sure orig nans' anoms are set back to nan's.
            end
            

            if (obj.RP.sigma_normalize && ~obj.sigma_normalized)
%                 if (obj.DAType=="model") 
%                     dd=1;  
%                 elseif (obj.DAType=="hist")
%                     dd=2;  
%                 else
%                     dd=3; 
%                 end
%                 ix=strfind(obj.DA_title,"_");
%                 ix0 = ix(end-3) + 1;
%                 ix1 = ix(end-1) + 1;
%                 ix2 = ix(end) + 1;
%                 ttl = sprintf("%s %s", obj.DAType,string(extractBetween(obj.DA_title,ix0,ix0+12)));
%                 figno=5000+100*dd*str2double(extractBetween(obj.DA_title,ix1, ix1+3)) +  str2double(extractBetween(obj.DA_title,ix2, ix2+3));
%                 [obj.anoms, obj.norm_sigmas] = obj.normalize_anoms_by_sigmas(obj.anoms, obj.DP.yrlen, true, [], ttl, figno);
                obj.orig_anoms = obj.anoms;         % save the original anomalies, un-normalized.
                [obj.anoms, obj.norm_sigmas] = obj.normalize_anoms_by_sigmas(obj.anoms, obj.DP.yrlen, true, []);
                obj.sigma_normalized = true;
            end
            

                        % create the mapping from model to obs data for historical period

            if (~obj.DP.isPrecipRun)
                obj.anoms(obj.na_map_in) = nan;        % set original missing data back to NAs
                if (~isempty(obj.orig_anoms))
                    obj.orig_anoms(obj.na_map_in) = nan;        % set original missing data back to NAs
                end
            else
                                                        % for precip, we've already copied the raw data into anoms array
    %           obj.raw_fix(obj.na_map_in) = nan;       % We leave the original NAs as zeros in raw_fix.  We don't
                                                        % actually use raw_fix for anything else in Precip runs.  
                                                        % (It is used in detrending for temperature data runs.)
            end


        end        
        
        function obj = calc_hists(obj)
        % function to calculate base and rolling histograms.  
            
                    % if anomalies not calculated yet, go calculate them first.
            if (isempty(obj.anoms))
                obj.calc_anomalies();
                if (obj.insufficient_data), return; end  % abort if insufficient data to process.
            end
            
            if (~obj.binning_checked)
                obj.calc_binning();
                if (obj.insufficient_data), return; end  % abort if insufficient data to process.
            end
            
            npts = numel(obj.raw_data);
            nyrs = npts/obj.RP.yrlen;
                                        % NOTE:  trying to disaggregate data with leap years won't work because we need to 
                                        % use fourier transforms, which needs all years to be the same length.  So we'll
                                        % bail out under all circumstances if data is isn't a multiple of the yearlen.
            if (mod(nyrs,1) ~= 0), obj.DP.error_log("ARRM_V2_DisaggregateSignal:  length(obj.raw_data) is not multiple of yrlen"); end      
%             if (nyrs ~= obj.RP.nyrs), error("ARRM_V2_DisaggregateSignal:  data length (%d yrs) doesn't match RunParams nyrs (%d)",nyrs,obj.RP.nyrs); end


            if (~isempty(obj.DP.base_yrs))
                [bstart, bend, nbase_yrs] = obj.using_range('base_yrs');               
                [obj.base_hist, obj.lohi_base] = histogram_anomalies(obj.anoms(bstart:bend), obj.na_map_in(bstart:bend), [1,nbase_yrs], [], nbase_yrs, obj.RP.edges, false, false, obj.RP.yrlen);
%               [obj.base_hist, obj.lohi_base] = histogram_anomalies(obj.anoms(bstart:bend), obj.na_map_in(bstart:bend), [1,nbase_yrs], [], nbase_yrs, obj.RP.edges, true,  true,  obj.RP.yrlen);
                if (~obj.check_valid_counts(obj.base_hist)), obj.insufficient_data = true; return; end  % abort if insufficient data to process.
            end
            
            if (~isempty(obj.DP.rolling_yrs))
                yr_start = obj.DP.rolling_yrs(1)-obj.DP.data_yrs(1)+1;      % year # to start at.  1-based (1 is 1st year.)
                yr_end   = obj.DP.rolling_yrs(2)-obj.DP.data_yrs(1)+1;    % year # to end at.
                yr_range = [yr_start, yr_end];

                obj.DP.rolling_steps = obj.DP.calc_rolling_steps(obj.RP.pdf_yrstep);
                nsteps = length(obj.DP.rolling_steps);
                rolling_range = zeros(nsteps,2);
                for j=1:nsteps
                    [bstart,bend] = obj.using_range("rolling_yrs",j);
                    rolling_range(j,:) = [bstart, bend];
                end
                
                if (obj.DP.isPrecipRun)
                    na_map = isnan(obj.anoms);
                    [obj.rolling_hist, obj.lohi_all] = histogram_anomalies(obj.anoms,     na_map,    yr_range, obj.DP.rolling_steps, obj.RP.pdf_yrs, obj.RP.edges, false, false, obj.RP.yrlen, rolling_range);
%                   [obj.rolling_hist, obj.lohi_all] = histogram_anomalies(obj.anoms,     na_map,    yr_range, obj.DP.rolling_steps, obj.RP.pdf_yrs, obj.RP.edges, true,  true,  obj.RP.yrlen);
                else
                    [obj.rolling_hist, obj.lohi_all] = histogram_anomalies(obj.anoms, obj.na_map_in, yr_range, obj.DP.rolling_steps, obj.RP.pdf_yrs, obj.RP.edges, true, true, obj.RP.yrlen);
                end
%               if (isempty(obj.rolling_hist)), return; end  % No check on valid count for rolling period...assume the  model has enough data.
            end
            
%--------------  code to mirror precip to create approx. gaussian.  This is currently not used, but left in in case it
%is needed again later.
                % for precip, all histogram values will be on positive size of hist.  We flip the hist to make it
                % symmetrical.  Edges are already flipped, going as far negative as positive.
                
                % could probably just append zeros, do circular convolution, and then add the two halves back
                % together...
%             if (obj.DP.isPrecipRun)
% %               obj.base_hist = [flipud(obj.base_hist);obj.base_hist];
%                 nr  = size(obj.base_hist,1);                
%                 mid = floor(nr/2);
%                 last = nr - mod(nr,2);
%                 obj.base_hist(1:mid,:) = flipud(obj.base_hist((mid+1):last,:));
%                 
%                 
%                 if (~isempty(obj.DP.rolling_yrs))
% %                   obj.rolling_hist = [flipud(obj.rolling_hist);obj.rolling_hist];
%                     obj.rolling_hist(1:mid,:,:) = flipud(obj.rolling_hist((mid+1):last,:,:));
%                 end
%                 
%  %              obj.RP.edges = [-fliplr(obj.RP.edges(1:end-1)), obj.RP.edges];  % mirroring the edges is now done in ARRM_V2_run_gridded
%             end
%                         
%    92.0000   -0.0473
%    23.0000   -0.0315
%          0   -0.0158
%          0         0
%          0    0.0158
%          0    0.0315
%    23.0000    0.0473
%    92.0000    0.0631
%    55.0000    0.0788                        
%             figure(99);       debugging...
%             rawfix = reshape(rawfix,365,57);
%             rd=obj.anoms+obj.trend+obj.avg_base+obj.moving_clim;
%             anms=obj.anoms;
%             anms(obj.na_map_in)=obj.moving_clim(obj.na_map_in);
%             rd2=anms+obj.trend+obj.avg_base+obj.moving_clim;
%             rd=reshape(rd,365,57);
%             rd2=reshape(rd2,365,57);
%             figure(99);
%             plot(1:365,mean(rawfix,2),'c:','linewidth',3);
%             hold on;
%             plot(1:365,obj.base_clim);
%             plot(1:365, nanmean(rd,2));
%             plot(1:365, mean(rd2,2));
%             hold off;
%-------------------
        end
        
        function calc_pdfs(obj, do_base, do_rolling)
            
            % Get pdf's and CDFs of anomalies
            
            if (obj.RP.pdf_yrlen < 1)
                obj.pdf_base = [];
                obj.CDF_base = [];
                obj.pdf_rolling = [];
                obj.CDF_rolling = [];
                return;
            end
                                        
            if (~exist('do_base','var')),    do_base    = true; end
            if (length(do_base) > 1)
                do_rolling = do_base(2);
                do_base=do_base(1);
            end
            if (~exist('do_rolling','var')), do_rolling = true; end

                % make sure we've got histograms before we do pdfs...
            if (do_base && ~isempty(obj.DP.base_yrs) && isempty(obj.base_hist))
                obj.calc_hists();
            end
            if (do_rolling && ~isempty(obj.DP.rolling_yrs) && isempty(obj.rolling_hist))
                obj.calc_hists();
            end
            if (obj.insufficient_data), return; end
            
            if (do_base && ~isempty(obj.DP.base_yrs))
                
                    % Ian:  might need to sum and repmat if pdf_yrlen is 1.
                
                if (obj.RP.pdf_yrlen>1 && (obj.RP.pdf_yrlen ~= obj.RP.yrlen || ~strcmp(obj.RP.pdf_map_method,'linear')))        % change time dimension?
%                    fprintf(obj.DP.fidlog, 'mapping basic\n');
                        % Yes.  remap to alternate time (day-of-year) mapping
%                     if (obj.DP.isPrecipRun)
%                         distrib = obj.wet_day_clim;
%                     else
   %                    distrib = calc_day_mapping_distribution(obj);
                        distrib = obj.calc_day_mapping_distribution();
%                     end
                    [mapped_hist, dayPos] = remap_time_dimension(obj.base_hist, obj.RP.pdf_yrlen, distrib);

%                   try
                        midpdf = obj.create_kde_pdfs_lpf(mapped_hist, dayPos, false);      % dayPos added so we can plot the final pdf back as 365-day.
 %                   catch
 %                       fprintf("oops!  error running create_kde_pdfs_lpf(...)\n");
 %                   end
                    
                        % map back to original time (day-of-year) scale
                    obj.pdf_base = remap_time_dimension(midpdf, obj.RP.yrlen, [], dayPos);
                    obj.CDF_base = cumsum(obj.pdf_base,'omitnan');
                    obj.CDF_base = obj.CDF_base ./ obj.CDF_base(end,:);     % make sure CDF's sum to 1.
                    obj.pdf_base = diff([zeros(1,obj.RP.yrlen); obj.CDF_base],1,1);    % to stay consistent with the CDF, in case of rounding near the top.
                else    % no time dimension warping.  Just create the pdfs & CDFs
                    if (obj.RP.pdf_yrlen == 1)
                        obj.RP.pdf_yrlen = obj.RP.yrlen;
                        obj.base_hist = repmat(mean(obj.base_hist,2),1,obj.RP.yrlen);    % kludge.  Should sum and do a single pdf/CDF...
                    end
                    [obj.pdf_base, obj.CDF_base] = obj.create_kde_pdfs_lpf(obj.base_hist, [], false);
                end
                
                if (any(obj.DP.cdf_append_pts(:) ~= 0))
                            % calculate means, std devs, skews & kurtosis
%                   [obj.means_base, obj.sigmas_base, obj.skewness_base, obj.kurtosis_base] = pdf_stats(setnans(obj.pdf_base,0), obj.RP.bins,     [], false);
                    [obj.means_base, obj.sigmas_base, obj.skewness_base, obj.kurtosis_base] = pdf_stats(setnans(obj.pdf_base,0), obj.RP.out_bins, [], false);
                    if (strcmpi(obj.DP.cdf_append_type,'empirical'))
                        [obj.pdf_base, obj.CDF_base, obj.CDF_okflags_base] = obj.append_basic_cdf(obj.pdf_base, obj.CDF_base, obj.means_base'); %, obj.sigmas_base);
                    elseif (strcmpi(obj.DP.cdf_append_type,'sk_normal'))
                        [obj.pdf_base, obj.CDF_base, obj.CDF_okflags_base] = obj.append_sk_normal_cdf(obj.pdf_base, obj.CDF_base, obj.means_base'); %, obj.sigmas_base');
                    else
                        [obj.pdf_base, obj.CDF_base, obj.CDF_okflags_base] = obj.append_normal_cdf(obj.pdf_base, obj.CDF_base, obj.means_base');    % transpose because append_normal needs year running down the column...
                    end
                else                
                    obj.CDF_okflags_base = obj.CDF_base >= obj.RP.cdf_thresh  &  obj.CDF_base <= 1.0-obj.RP.cdf_thresh;
                end
                
                    %calculate stats (or recalculate...should have changed slightly after appending corrected tails.)
%               [obj.means_base, obj.sigmas_base, obj.skewness_base, obj.kurtosis_base] = pdf_stats(setnans(obj.pdf_base,0), obj.RP.bins,     [], false);
   %            if (~isempty(obj.RP.prcp_distrib) && strlength(obj.RP.prcp_distrib)>0)
                if (obj.RP.isPrecipRun && any(strcmp(obj.RP.prcp_distrib,["generalizedpareto", "inversegaussian", "loglogistic", "lognormal"])))
                    [obj.means_base, obj.sigmas_base, obj.skewness_base, obj.kurtosis_base] = pdf_stats(setnans(obj.pdf_base,0), obj.RP.zbins,     [], false);
                else
                    [obj.means_base, obj.sigmas_base, obj.skewness_base, obj.kurtosis_base] = pdf_stats(setnans(obj.pdf_base,0), obj.RP.out_bins, [], false);
                end            
                pthresh_probs = sort([obj.RP.far_outlier_anchor_pt, .5, 1-obj.RP.far_outlier_anchor_pt]);  % sort in case user specified anchor point > .5
                obj.pthresh_base = find_problines(obj.pdf_base, obj.CDF_base, pthresh_probs, obj.RP.out_edges, false, 1, obj.RP.cdf_thresh, false); % NOTE:  no need to filter using nterms & sig_terms;  we've already done that.      
                
            end

                    % get pdf's and CDFs for moving windows
            if (do_rolling && ~isempty(obj.DP.rolling_yrs))
                
                if (obj.RP.pdf_yrlen>1 && (obj.RP.pdf_yrlen ~= obj.RP.yrlen || ~strcmp(obj.RP.pdf_map_method,'linear')))        % change time dimension?
%                   fprintf(obj.DP.fidlog, 'mapping rolling\n');
                        % Yes.  remap to alternate time (day-of-year) mapping
                        
                    nsets = size(obj.rolling_hist, 3);
%                     if (obj.DP.isPrecipRun)
%                         distrib = repmat(obj.wet_day_clim, nsets, 1);
%                     else
%                       distrib = calc_day_mapping_distribution(obj, true);
                        distrib = obj.calc_day_mapping_distribution(true);
%                     end
                    try 
                        [mapped_hist, dayPos] = remap_time_dimension(obj.rolling_hist, obj.RP.pdf_yrlen, distrib);
                    catch me
                        report_me_error(me);
                        fprintf("oops\n");
                        rethrow(me);
                    end
                    
                    midpdf  = obj.create_kde_pdfs_lpf(mapped_hist, dayPos, true);

                        % map back to original time (day-of-year) scale

                    obj.pdf_rolling = remap_time_dimension(midpdf, obj.RP.yrlen, [], dayPos);
                    obj.CDF_rolling = zeros(size(obj.pdf_rolling));
                    for j=1:nsets
                        obj.CDF_rolling(:,:,j) = cumsum(obj.pdf_rolling(:,:,j),'omitnan');
                        obj.CDF_rolling(:,:,j) = obj.CDF_rolling(:,:,j) ./ obj.CDF_rolling(end,:,j);     % make sure CDF's sum to 1.
                    end
                    obj.pdf_rolling = diff([zeros(1,obj.RP.yrlen,nsets); obj.CDF_rolling],1,1);    % to stay consistent with the CDF, in case of rounding near the top.
                else
                    if (obj.RP.pdf_yrlen == 1)
                        nsets = size(obj.rolling_hist, 3);
                        obj.RP.pdf_yrlen = obj.RP.yrlen;
                        for i=1:nsets
                            obj.rolling_hist(:,:,i) = repmat(mean(obj.rolling_hist(:,:,i),2),1,obj.RP.yrlen);    % kludge.  Should sum and do a single pdf/CDF...
                        end
                    end
                    [obj.pdf_rolling, obj.CDF_rolling] = obj.create_kde_pdfs_lpf(obj.rolling_hist,[], true);
                end
                
                if (any(obj.DP.cdf_append_pts(:) ~= 0))
                            % calculate stats
%                   [obj.means_rolling, obj.sigmas_rolling, obj.skewness_rolling, obj.kurtosis_rolling] = pdf_stats(setnans(obj.pdf_rolling,0), obj.RP.bins,     [], false);
                    [obj.means_rolling, obj.sigmas_rolling, obj.skewness_rolling, obj.kurtosis_rolling] = pdf_stats(setnans(obj.pdf_rolling,0), obj.RP.out_bins, [], false);
                    if (strcmp(obj.DP.cdf_append_type,'empirical'))
                        [obj.pdf_rolling, obj.CDF_rolling, obj.CDF_okflags_rolling] = obj.append_basic_cdf(obj.pdf_rolling, obj.CDF_rolling, obj.means_rolling); %, obj.sigmas_rolling);
                    elseif (strcmpi(obj.DP.cdf_append_type,'sk_normal'))
                        [obj.pdf_rolling, obj.CDF_rolling, obj.CDF_okflags_rolling] = obj.append_sk_normal_cdf(obj.pdf_rolling, obj.CDF_rolling, obj.means_rolling); %, obj.sigmas_rolling);
                    else
                        try
                            [obj.pdf_rolling, obj.CDF_rolling, obj.CDF_okflags_rolling] = obj.append_normal_cdf(obj.pdf_rolling, obj.CDF_rolling, obj.means_rolling);
                        catch
                            fprintf("oops appending normal cdf");
                        end
                    end
                else
                    obj.CDF_okflags_rolling = obj.CDF_rolling > obj.RP.cdf_thresh  &  obj.CDF_rolling < 1.0-obj.RP.cdf_thresh;
                end
                
                    % recalculate stats;  should have changed slightly after appending corrected tails.
%               [obj.means_rolling, obj.sigmas_rolling, obj.skewness_rolling, obj.kurtosis_rolling] = pdf_stats(setnans(obj.pdf_rolling,0), obj.RP.bins,     [], false);
                
%               if (~isempty(obj.RP.prcp_distrib) && strlength(obj.RP.prcp_distrib)>0)
                if (obj.RP.isPrecipRun && any(strcmp(obj.RP.prcp_distrib,["generalizedpareto", "inversegaussian", "loglogistic", "lognormal"])))
                    [obj.means_rolling, obj.sigmas_rolling, obj.skewness_rolling, obj.kurtosis_rolling] = pdf_stats(setnans(obj.pdf_rolling,0), obj.RP.zbins,    [], false);
                else
                    [obj.means_rolling, obj.sigmas_rolling, obj.skewness_rolling, obj.kurtosis_rolling] = pdf_stats(setnans(obj.pdf_rolling,0), obj.RP.out_bins, [], false);
                end
                
                pthresh_probs = sort([obj.RP.far_outlier_anchor_pt, .5, 1-obj.RP.far_outlier_anchor_pt]); % sort in case user specified anchor point > .5
                obj.pthresh_rolling = zeros(obj.RP.yrlen,3,length(obj.DP.rolling_steps));
                for j=1:length(obj.DP.rolling_steps)
                    obj.pthresh_rolling(:,:,j) = find_problines(obj.pdf_rolling(:,:,j), obj.CDF_rolling(:,:,j), pthresh_probs, obj.RP.out_edges, false, 1, obj.RP.cdf_thresh); % NOTE:  no need to filter using nterms & sig_terms;  we've already done that.

                end
            end
        end
        
        function calc_problines(obj, do_base, do_rolling)
            
            % Get pdf's and CDFs of anomalies
                        
            if (~exist('do_base','var')),    do_base = true; end
            if (length(do_base) > 1)
                do_rolling = do_base(2);
                do_base=do_base(1);
            end
            if (~exist('do_rolling','var')), do_rolling = true; end

                % make sure we've got pdfs before we do problines...
            if (do_base && ~isempty(obj.DP.base_yrs) && isempty(obj.pdf_base))
                obj.calc_pdfs(do_base, do_rolling);
            end
            if (do_rolling && ~isempty(obj.DP.rolling_yrs) && isempty(obj.pdf_rolling))
                obj.calc_pdfs(do_base, do_rolling);
            end
            if (obj.insufficient_data), return; end
            
            if (do_base && ~isempty(obj.DP.base_yrs))
                
                [obj.problines_base,  obj.pdf_xlines_base,  obj.CDF_xlines_base] = find_problines(obj.pdf_base,  obj.CDF_base,  obj.RP.probs, obj.RP.edges, false,  1, obj.RP.cdf_thresh, false); % NOTE:  no need to filter using nterms & sig_terms;  we've already done that.
                [bstart,bend] = obj.using_range('base_yrs');
                [obj.probcounts_base, obj.totcount_base] = find_probcounts(obj.problines_base,obj.anoms(bstart:bend,:));
            end

                    % get pdf's and CDFs for moving windows
            if (do_rolling && ~isempty(obj.DP.rolling_yrs))
                
                
                    % Calculate problines on rolling pdfs, CDFs if requested
                nsets  = size(obj.pdf_rolling, 3);
                nprobs = length(obj.RP.probs);
                obj.problines_rolling  = nan(obj.RP.yrlen, nprobs, nsets);
                obj.pdf_xlines_rolling = nan(obj.RP.yrlen, nprobs, nsets);
                obj.CDF_xlines_rolling = nan(obj.RP.yrlen, nprobs, nsets);
                obj.probcounts_rolling = nan(nprobs,nsets);
                obj.totcounts_rolling  = nan(nsets,1);
                for j=1:nsets
                    [obj.problines_rolling(:,:,j),  obj.pdf_xlines_rolling(:,:,j),  obj.CDF_xlines_rolling(:,:,j)] = find_problines(obj.pdf_rolling(:,:,j), obj.CDF_rolling(:,:,j),  obj.RP.probs, obj.RP.edges, false, 1, obj.RP.cdf_thresh); % again, no need to filter.
%                       [obj.probcounts_rolling(:,i), obj.totcounts_rolling(i)]  = find_probcounts(obj.problines_rolling(:,:,i), obj.RP.probs, squeeze(obj.rolling_hist(:,:,i))', obj.RP.bins);
                    [bstart,bend] = obj.using_range('rolling_yrs', j);
                    [obj.probcounts_rolling(:,j), obj.totcounts_rolling(j)]  = find_probcounts(obj.problines_rolling(:,:,j), obj.anoms(bstart:bend,:));
                end
                
            end
        end
        
        function [probs, far_outliers, outliers] = calc_base_probs(obj)
            
            % Get pdf's and CDFs of anomalies
            
            if (isempty(obj.pdf_base))
                obj.calc_pdfs(true, false);
            end
            if (obj.insufficient_data), probs=[]; outliers=[];  far_outliers=[]; return; end

            prcp_min = obj.DP.prcp_min;
            prcp_medians = obj.pthresh_base(:,2);
            
            [bstart,bend] = obj.using_range('base_yrs');
            [obj.anom_probs_base, mapped_nas, far_outliers, outliers]    = calc_anom_probs(obj.anoms(bstart:bend), obj.RP.out_bins, obj.CDF_base, obj.na_map_in(bstart:bend), obj.RP.far_outlier_thresh, obj.RP.outlier_thresh, obj.DP.isPrecipRun, prcp_min,prcp_medians);
            obj.anom_sdevs_base = norminv(obj.anom_probs_base);     % std. dev equivalents for each anomaly
            obj.insert_data(far_outliers, 'far_outlier_map','base_yrs');
            obj.insert_data(outliers,'outlier_map','base_yrs');
            
            obj.anom_probs_base(mapped_nas) = nan;       % for now.  may want to change this, ian.
            probs = obj.anom_probs_base;
        end
       
        function [probs, far_outliers, outliers] = calc_rolling_probs(obj)
           
            if (isempty(obj.pdf_rolling))
                obj.calc_pdfs(false, true);
            end
            if (obj.insufficient_data), probs=[];  outliers=[];  far_outliers=[]; return; end
            
            prcp_min = obj.DP.prcp_min;
            prcp_medians = obj.pthresh_base(:,2);
            
            for j=1:length(obj.DP.rolling_steps)
                [bstart,bend] = obj.using_range('rolling_yrs',j);
                if (obj.RP.sigma_normalize)
                    mybins = obj.RP.out_bins;
                else
                    mybins = obj.RP.bins;
                end
%               [probs_rolling, mapped_nas, far_outliers, outliers] = calc_anom_probs(obj.anoms(bstart:bend), obj.RP.out_bins, obj.CDF_rolling(:,:,j), obj.na_map_in(bstart:bend), obj.RP.far_outlier_thresh, obj.RP.outlier_thresh,  obj.DP.isPrecipRun,  prcp_min, prcp_medians);
%               [probs_rolling, mapped_nas, far_outliers, outliers] = calc_anom_probs(obj.anoms(bstart:bend), obj.RP.bins,     obj.CDF_rolling(:,:,j), obj.na_map_in(bstart:bend), obj.RP.far_outlier_thresh, obj.RP.outlier_thresh,  obj.DP.isPrecipRun,  prcp_min, prcp_medians);
                [probs_rolling, mapped_nas, far_outliers, outliers] = calc_anom_probs(obj.anoms(bstart:bend), mybins,          obj.CDF_rolling(:,:,j), obj.na_map_in(bstart:bend), obj.RP.far_outlier_thresh, obj.RP.outlier_thresh,  obj.DP.isPrecipRun,  prcp_min, prcp_medians);
                obj.insert_data(far_outliers, 'far_outlier_map','rolling_yrs',j);
                obj.insert_data(outliers,'outlier_map','rolling_yrs',j);
                probs_rolling(mapped_nas) = nan;       % for now.  may want to change this, ian.
                obj.insert_data(probs_rolling,'anom_probs_rolling','rolling_yrs',j);
%               obj.anom_probs_rolling(bstart:bend) = probs_rolling;
            end
            obj.anom_sdevs_rolling = norminv(obj.anom_probs_rolling);     % std. dev equivalents for each anomaly
            probs = obj.anom_probs_rolling;
       end        
        
        function assembled = reassemble(obj, anoms, delta_hist_mdl)
            if (~exist('anoms','var')), anoms = obj.anoms; end
            if (~exist('delta_hist_mdl','var')), delta_hist_mdl = 0; end    % difference between hist avg_base and mdl_avg_base.  
                                                                            % this will only be non-zero for separate_hist runs.
            if (isempty(obj.trend))
                if (obj.DP.isPrecipRun)
                    assembled = anoms(:);
                else
                    assembled =  obj.moving_clim(:) + obj.avg_base(:) + anoms(:) + delta_hist_mdl; 
                end
            else
                assembled =  obj.moving_clim(:) + obj.trend(:) + obj.avg_base + anoms(:) + delta_hist_mdl;   % (:) needed here in case we're adding both rows and columns.
            end
        end
        
        function insert_data(obj, data, field, yr_typ, rolling_set)     % removed "do_clear" flag...
            if (~exist('rolling_set','var')), rolling_set = []; end

            if (~isempty(data))
                [ixstart,ixend] = obj.using_range(yr_typ, rolling_set);
                if (isempty(obj.(field)))         % if empty, create column vector so we get a column vector as output.
                    if (~islogical(obj.(field)))
                        obj.(field) = zeros(ixend,1); 
                    else
                        obj.(field) = false(ixend,1); 
                    end
                end
                obj.(field)(ixstart:ixend) = data;
            end
        end           
        
        function [edges, edgerange, nedges, dx] = calc_binning(obj)
            % Purpose here is to find a reasonable range for our PDFs, and sufficient binning over that range to make
            % sure our convolutions go to zero in the wings, and that we don't have any discontinuities when we move
            % into quantile space.  
            % One-size won't fit all here.  Narrow distributions need a smaller range with closer edges to avoid
            % discontinuities in quantile space, while wider distributions need wider range, but if edges are too close, 
            % then there are so many edges that the processing slows down.
            
            if (obj.DP.isPrecipRun)
%               [edges, edgerange, nedges, dx] = obj.calc_precip_binning();
                edges = obj.prcp_scaling.edges(2:end);
                obj.RP.edges = edges;
                edgerange = [edges(1),edges(end)];
                nedges = length(edges);
%               dx = obj.prcp_scaling.dx;                
                obj.binning_checked = true;
                return;
            end

            
            if (~isempty(obj.DP.base_yrs))
                [edgemin, edgemax, nedges, ssig_minsig] = calc_binning_sub(obj, 'base_yrs');
            else
                [edgemin, edgemax, nedges, ssig_minsig] = calc_binning_sub(obj, 'rolling_yrs');
            end
                
      %     obj.DP.print_log("base_period: done calc_binning_sub  edgemin %.4f edgemax %.4f nedges %d ssig_minsig %.4f\n", edgemin, edgemax, nedges, ssig_minsig);
            
            nsets = length(obj.DP.rolling_steps);
            for i=1:nsets
                [emin, emax, nes, ssmin] = calc_binning_sub(obj, 'rolling_yrs', i);
                edgemin  = min(edgemin, emin);
                edgemax  = max(edgemax, emax);
                nedges   = max(nedges, nes);
                ssig_minsig = min(ssig_minsig, ssmin);
      %         obj.DP.print_log("rolling period: set %d: done calc_binning_sub  edgemin %.4f edgemax %.4f nedges %d ssig_minsig %.4f\n", i, edgemin, edgemax, nedges, ssig_minsig);
            end
            
            edgerange = [edgemin,edgemax];
            edgerng = range(edgerange);
            
            dx = edgerng/(nedges-1);
            ssig_min =   silversig_min(nedges-1);
            
            n=0;
            while(ssig_min > .5*ssig_minsig/dx)     % Ian - you should make this a binary search to speed it up.
                nedges = nedges+50;
                dx = edgerng/(nedges-1);
                ssig_min =  silversig_min(nedges-1);
%               obj.DP.print_log( '%d %4d %8.4f > %8.4f \n', n, nedges, ssig_min, ssig_minsig/dx);
                n=n+1;
                if (mod(n,50)==0), fprintf("calc_binning:  n  %d  nedges  %d  dx  %.6f   ssig_min %.4f  ssigm_minsig val:  %.4f\n", n, nedges, dx, ssig_min, .75*ssig_minsig/dx); end
                if (n>500)
                    obj.DP.warn_log("calc_binning:  range calculating oops  %s %s %s %s : dx %.6f ssig_min: %.4f  ssig_minsig val: %.4f\n", ...
                                    obj.DP.model, obj.DP.ensemble, obj.DP.varname, obj.DP.scenario, dx, ssig_min, .75*ssig_minsig/dx);
                    obj.DP.error_log("calc_binning:  range calculating oops  %s %s %s %s dx %.6f ssig_min: %.4f  ssig_minsig val: %.4f\n", ...
                                    obj.DP.model, obj.DP.ensemble, obj.DP.varname, obj.DP.scenario, dx, ssig_min, .75*ssig_minsig/dx);
                end
            end
            
            edges = linspace(edgemin, edgemax, nedges);
            
            obj.RP.edges = edges;
            obj.binning_checked = true;
            
            edgerange = [edgemin, edgemax];
%           obj.DP.print_log('binning range:     %9.5f to %9.5f, %6d bins, ssig_min: %9.5f   %s\n', edgerange, nedges, ssig_min, obj.DP.runLbl);
            
        end   
        
        function [edgemin, edgemax, nedges, ssig_minsig] = calc_binning_sub(obj, yrlbl, set)
            % Purpose here is to find a reasonable range for our PDFs, and sufficient binning over that range to make
            % sure our convolutions go to zero in the wings, and that we don't have any discontinuities when we move
            % into quantile space.  Some of the numbers below are empirically determined...

            if (~exist('set','var'))
                [bstart,bend, nyrs] = obj.using_range(yrlbl);
            else
                [bstart,bend, nyrs] = obj.using_range(yrlbl, set);
            end
            my_anoms = obj.anoms(bstart:bend); 
            
            my_anoms = reshape(my_anoms, obj.RP.yrlen, nyrs);
            sigs = std(my_anoms, 1,2, "omitnan");
            sigs  = fillmissing(sigs, 'nearest','endvalues','nearest');
            clim_sigs = climatology(sigs, obj.RP.anom_nterms(1), obj.RP.anom_sig_terms(1), obj.RP.yrlen);
            
            minv = min(my_anoms(:),[],"omitnan");
%             minix = find(my_anoms(:)==minv,1);
            maxv = max(my_anoms(:),[],"omitnan");
%             maxix = find(my_anoms(:) == maxv,1);
%           obj.DP.print_log("\t\t\t\t\tcalc_binning_sub:  minv = %10.6f   maxv = %10.6f\n", minv, maxv);
              
%           npts = obj.RP.pdf_yrs * 365/obj.RP.pdf_yrlen;      % will need to fix this for precip, Ian!
            npts = nyrs * obj.RP.yrlen/obj.RP.pdf_yrlen;      % will need to fix this for precip, Ian!
%             mns_doy = find(clim_sigs == min(clim_sigs),1);
%             min_doy = mod(minix-1, obj.RP.yrlen)+1;
%             max_doy = mod(maxix-1, obj.RP.yrlen)+1;
%                 ssig_minv   = 1.06 * clim_sigs(min_doy) .* (npts.^(-.2));
%                 ssig_maxv   = 1.06 * clim_sigs(max_doy) .* (npts.^(-.2));
%                 ssig_minsig = 1.06 * clim_sigs(mns_doy) .* (npts.^(-.2));
            ssig_max    = 1.06 * max(clim_sigs) .* (npts.^(-.2));
            ssig_minsig = 1.06 * min(clim_sigs) .* (npts.^(-.2));
                       
            if(obj.RP.do_calc_binning)
                edgemin = minv;
                edgemax = maxv;
%               nedges  = 280;
                nedges  = 501;
%               obj.DP.print_log("do_calc_binning is true.  starting nedges = %d\n", nedges);
            else
                edgemin = obj.RP.edges(1);
                edgemax = obj.RP.edges(end);
                nedges  = length(obj.RP.edges);
%               obj.DP.print_log("do_calc_binning is false.  starting nedges = %d\n", nedges);
            end
            
                % make sure we go at least 9.5 sigmas past min & max values in the data
                % so convolution will go to essentially 0.
%           edgemin = min(edgemin, 1.5*((-9.5)*ssig_max+minv));
%           edgemax = max(edgemax, 1.5*(( 9.5)*ssig_max+maxv));
            edgemin = min(edgemin, 1.5*((-7.5)*ssig_max+minv));
            edgemax = max(edgemax, 1.5*(( 7.5)*ssig_max+maxv));
            edgerange = [edgemin,edgemax];
            edgerng = range(edgerange);
            
            dx = edgerng/(nedges-1);
            ssig_min =   silversig_min(nedges-1);
            
            n=0;
                % now increase the number of edges until our silverman's rule-of-thumb sigma is small enough to
                % guarantee we don't have discontinuities in our convolved surfaces.
            while(ssig_min > .75*ssig_minsig/dx)
                nedges = nedges+20;
                dx = edgerng/(nedges-1);
                ssig_min =  silversig_min(nedges-1);
 %              fprintf(obj.DP.fidlog, '%d %4d %8.4f > %8.4f \n', n, nedges, ssig_min, ssig_minsig/dx);
                n=n+1;
%               if (mod(n,50)==0), fprintf("calc_binning_sub:  n  %d  nedges  %d  dx  %.6f   ssig_min %.4f  ssigm_minsig val:  %.4f\n", n, nedges, dx, ssig_min, .75*ssig_minsig/dx); end
                if (n>500)
%                     obj.DP.warn_log("calc_binning_sub:  range calculating oops  %s %s %s %s\n", obj.DP.model, obj.DP.ensemble, obj.DP.varname, obj.DP.scenario);
%                     error("calc_binning_sub:  range calculating oops  %s %s %s %s\n", obj.DP.model, obj.DP.ensemble, obj.DP.varname, obj.DP.scenario);
                    obj.DP.warn_log("calc_binning_sub:  %s : too many iterations calculating ssig_min and nbins  %s %s %s %s : dx %.6f ssig_min: %.4f  ssig_minsig val: %.4f\n", ...
                                    obj.DP.stnName, obj.DP.model, obj.DP.ensemble, obj.DP.varname, obj.DP.scenario, dx, ssig_min, .75*ssig_minsig/dx);
                    obj.DP.error_log("calc_binning:  %s : too many iterations calculating ssig_min and nbins  %s %s %s %s dx %.6f ssig_min: %.4f  ssig_minsig val: %.4f\n", ...
                                    obj.DP.stnName, obj.DP.model, obj.DP.ensemble, obj.DP.varname, obj.DP.scenario, dx, ssig_min, .75*ssig_minsig/dx);
                end
            end
            
%             if (~exist('set','var'))
%                 setlbl = 'base   ';
%             else
%                 setlbl = sprintf('step %2d', set);
%             end
%           obj.DP.print_log('calc_binning_sub:  calculated range: %s %s:  %.4f to %.4f, %d bins, ssig_min: %.4f\n', obj.DP.runLbl, setlbl, edgerange, nedges, ssig_min)
           
        end
        
        function distrib = calc_day_mapping_distribution(obj, do_rolling)
            
            if (~exist('do_rolling','var') || ~do_rolling)        % calculate for basic 

                if (obj.DP.isPrecipRun)
%                   clim = obj.wet_day_clim;
                    clim = obj.wet_day_clim_raw;
                else
                    clim = obj.base_clim;
                end
                clim = max(clim,0);     % in really dry locations, the filtered climatology can go negative...
                if (strncmpi(obj.RP.pdf_map_method, 'cos',3))
                    slopes = climatology(circular_slopes(clim),6,2);
                    distrib = 1+sin(abs(slopes)*pi/4);
                elseif (strncmpi(obj.RP.pdf_map_method, 'lin',3))
                    distrib = ones(1,obj.RP.yrlen);
                elseif (strncmpi(obj.RP.pdf_map_method, 'clim',4))
                    distrib = clim;
                elseif (strcmpi(obj.RP.pdf_map_method,"wetclim"))
%                   slopes = climatology(circular_slopes(clim),6,2);
                    slopes = circular_slopes(clim);
                    distrib = clim .*(1+sin(abs(slopes)*pi/4))';
                end
               
            else
                nsets = length(obj.DP.rolling_steps);
                if (strncmpi(obj.RP.pdf_map_method, 'lin',3))
                    distrib = ones(nsets,obj.RP.yrlen);
                else                  
                    distrib = zeros(nsets,obj.DP.yrlen);                                       
                    for j=1:nsets
                        if (obj.DP.isPrecipRun)
                            if (isempty(obj.wet_day_rolling_clim_raw))
%                               clim = obj.wet_day_clim;
                                clim = obj.wet_day_clim_raw;
                            else
                                clim = obj.wet_day_rolling_clim_raw(j,:);
                            end
                        else
                            if (isempty(obj.rolling_clim))
                                clim = obj.base_clim;
                            else
                                clim = obj.rolling_clim(j,:);
                            end
                        end
                        
                        if (strncmpi(obj.RP.pdf_map_method, 'cos',3))
                            slopes = climatology(circular_slopes(clim),6,2);                           
                            distrib(j,:) = 1+sin(abs(slopes)*pi/4);
                        elseif (strncmpi(obj.RP.pdf_map_method, 'clim',4))
                            distrib(j,:) = clim;
                        elseif (strcmpi(obj.RP.pdf_map_method,"wetclim"))                            
 %                          slopes = climatology(circular_slopes(clim),6,2);                           
                            slopes = circular_slopes(clim);                           
                            distrib = clim .*(1+sin(abs(slopes)*pi/4))';                            
                        end
                    end 
                end
            end            
        end
        
        function trim_data_yrs(obj, trim_yrs, keep_all, trim_fields)
            % trims all years to trim_yrs.  
            % if missing or empty, trims to DA.data_yrs.
            % updates DP.data_final_yrs.
            
            if (~exist('trim_yrs','var') || isempty(trim_yrs)), trim_yrs = obj.DP.data_final_yrs; end
            if (~exist('keep_all','var') || isempty(keep_all)), keep_all = obj.DP.keep_all; end
            if (isempty(trim_yrs))
                trim_yrs = obj.DP.data_yrs;
            else
                trim_yrs = min_yr_range(trim_yrs, obj.DP.data_yrs);
            end
            [d_start, d_end] = obj.using_range(trim_yrs);
            
            if (keep_all)       % this needs updating, Ian!
                obj.anoms_all = obj.anoms;
                obj.trend_all = obj.trend;
                obj.moving_clim_all = obj.moving_clim;
                obj.moving_clim_1D_all = obj.moving_clim_1D;
            end
            
            if (~exist('trim_fields','var') || isempty(trim_fields))
                trim_fields = ["anoms","na_map_in","outlier_map","far_outlier_map","na_map_mapped","trend", ... 
                               "moving_clim", "moving_clim_1D", "raw_data", "raw_fix", ...
                               "mapped_anoms_im","outlier_mapped_anoms","mapped_output", "anom_probs_rolling","anom_sdevs_rolling"];
            end
            
            for i=1:length(trim_fields)
                prop = trim_fields{i};
                if (~isempty(obj.(prop)))
                    obj.(prop) = obj.(prop)(d_start:d_end);
                end
            end
                       
            obj.DP.data_final_yrs = trim_yrs;
            
            obj.DP = obj.DP.set_yr_limits(trim_yrs);
        end
        
        function DA_struct = toStruct(obj)
            pdf_props  = ["pdf_base","CDF_base","CDF_okflags_base","pdf_rolling","CDF_rolling","CDF_okflags_rolling"];
            prob_props = ["problines_base","pdf_xlines_base","CDF_xlines_base","probcounts_base","problines_rolling","pdf_xlines_rolling","CDF_xlines_rolling","probcounts_rolling"];
            keep_all_props = ["anoms_all","trend_all","moving_clim_all","moving_clim_1D_all"];
            props = properties(obj);
            DA_struct = struct();
            for i=1:length(props)
                prop = props{i};
                if (isa(obj.(prop),'ARRM_V2_RunParams') || isa(obj.(prop),'ARRM_V2_DataParams'))
                    DA_struct.(prop) = obj.(prop).toStruct();
                elseif (contains(pdf_props, prop))
                    if (obj.RP.do_pdfs)
                        DA_struct.(prop) = obj.(prop);
                    end
                elseif (contains(prob_props,prop))
                    if (obj.RP.do_probs)
                        DA_struct.(prop) = obj.(prop);
                    end
                elseif (any(strcmpi(keep_all_props,prop)) && obj.DP.keep_all)
                    DA_struct.(prop) = obj.(prop);
                else
                    DA_struct.(prop) = obj.(prop);
                end
            end
        end
        
        function obj = clean(obj)
            mc=?ARRM_V2_DisaggregateSignal;  % get metadata about DA class.
            props = properties(obj);
            for i=1:length(props)
                if (~mc.PropertyList(i).NonCopyable)
                    obj.(props{i}) = [];
                end        
            end
        end
        
        function [pdfs, CDFs] = create_kde_pdfs_lpf(obj, hist_tbl, dayPos, using_rolling)   
        %
        %   This creates pdf's from the histograms, using low-pass filtering.  Similar to quantile regression, but
        %   instead of using regression to calulate the shape of the quantile line, it uses low-pass filtering,
        %   essentially finding the 'climatology' of each quantile line.
        %   It warps the CDF into probability space (histogram of anomalies spaced evenly in probability space (std.
        %   deviations) rather than evenly spaced in temperature).  Then low-pass-filters along each probability line,
        %   and then warps the CDF back to anomaly vs day of year.
        %   outputs the smoothed PDF surface.  
        %   Calculates CDFs from pdfs only if nargout > 1.
        %   
%            runLbl = obj.DP.runLbl;
             runLbl="";
            
            [~,pdf_yrlen,nsets] = size(hist_tbl);

                % this creates a 3D matrix with dimensions (yvals, day-of-year, #years) , 
                % where #years is total years / yr_step -- or think of it as 1 2D matrix (nyvals x pdf_yrlen) for each 5 year period.
%           if (~isempty(obj.RP.prcp_distrib) && strlength(obj.RP.prcp_distrib)>0)
            if (obj.RP.isPrecipRun && any(strcmp(obj.RP.prcp_distrib,["generalizedpareto", "inversegaussian", "loglogistic", "lognormal"])))
                bins = to_column(obj.RP.zbins);
            else
                bins = to_column(obj.RP.bins);
            end

%             nterms = obj.RP.anom_nterms(end);
%             sig_terms = obj.RP.anom_sig_terms(end);
%             nterms = obj.RP.clim_nterms;
%             sig_terms = obj.RP.clim_sig_terms;
%             nterms = obj.RP.anom_nterms;
%             sig_terms = obj.RP.anom_sig_terms;

            sdrange=[-7.5,7.5];
            sdsteps = 50*range(sdrange)+1;
            sdvals = linspace(sdrange(1), sdrange(2), sdsteps);
            sdprbs = normcdf(sdvals);
            sdlines = -7:7; %sdvals(11:50:end-10);         % for plotting, plots lines at each std dev.
            sdprbl = normcdf(sdlines);
            
%             sdprbs = normcdf(zbins);
%             sdlines = -7:7;
%             sdprbl = normcdf(sdlines);
            
            pdfs = [];
            CDFs = [];
            
            if (any(strcmp(obj.RP.runType,'temp')))     % for labeling figures.
                ylbl = 'Deg C';
            elseif (obj.DP.isPrecipRun)
                probs_flag  = 1;
                ylbl = 'Precip, scaled';
            else
                probs_flag = 3;                 
                ylbl = '';
            end
            for j=1:nsets
                if (obj.DP.intern_plotflag(j))                 % see notes in ARRM_V2_DataParams on creating internal plots.
                    fignum = 2*(j-1)+1+5*obj.DP.figbase; 
                end
                pdfs_t = squeeze(hist_tbl(:,:,j));
                
                nnegs = sum(pdfs_t(:)<0);
                if ( nnegs > 0)
                    obj.DP.warn_log(sprintf("negatives in raw histogram:  %3d  min %.5g\n", nnegs, min(pdfs_t(:))));
                end
                
                    % see notes in ARRM_V2_DataParams on creating internal plots.
                if (obj.DP.intern_plotflag(j))
                    if (~using_rolling)
                        raw_hist = obj.base_hist;
                    else
                        raw_hist = squeeze(obj.rolling_hist(:,:,j));
                    end
                    lbl = sprintf("%s %s raw histogram", obj.DAType, obj.DA_yrlbl);
                    if (nsets > 1)
                        lbl = sprintf("%s, %d", lbl, ceil(mean(obj.rolling_step_years(j))));
                    end
                    plot_surface(fignum, raw_hist, 1:obj.RP.yrlen, bins, false, [3,3,1],lbl,                       true, 'time', ylbl, 'counts'); 
                    plot_surface(fignum,   pdfs_t, 1:pdf_yrlen,    bins, false, [3,3,4],'time-remapped histogram',false, 'time', ylbl, 'counts');
                end
                
%                myplines = find_problines(pdfs_t, [], obj.RP.probs, obj.RP.edges, false, 1, 1e-14,false,3,1);
                
                    % first, do KDE on the histogram surface
                [pdfs_t,~,clim_sigs_orig, silversigs]=kde_surf(obj, pdfs_t, bins);
                
                nnegs = sum(pdfs_t(:)<0);
                if ( nnegs > 0)
                    fprintf("negatives after kde_surf:  %3d  min %.5g\n", nnegs, min(pdfs_t(:)));
                end
                
%                mykplines = find_problines(pdfs_t, [], obj.RP.probs, obj.RP.edges, false, 1, 1e-14,false,3,1);
                
%                plot_kdelines(99, myplines, mykplines, obj.RP.probs);

                    % see notes in ARRM_V2_DataParams on creating internal plots.
%               if (obj.DP.intern_plotflag(j)), plot_surface(fignum, pdfs_t, 1:pdf_yrlen, bins, false, [3,3,4],'kde_surface',false, 'time', ylbl, 'counts'); end
                
%               x = to_column(obj.RP.bins);
%                 mymus = sum(x .* pdfs_t);                       % mean daily value, s/b around midpoint of bins.
%                 myclim_mus = obj.DA_climatology(mymus(1,:),nterms,sig_terms, pdf_yrlen)';
                
%                 if (obj.RP.do_normalize_sigmas)
% %                   mu0 = mean(mus(1,:));                               % center point of distributions
% %                   bin_ctr = interp1(obj.RP.bins, 1:nbins, mu0);                  % equivalent bin (fractional).
%                 
%                     mysigs = sqrt(sum(x.^2 .* pdfs_t) - mymus.^2); 
%                     myclim_sigs = obj.DA_climatology(mysigs,nterms,sig_terms, pdf_yrlen)';   %<---- should be changed if we keep this, Ian.
% %                   sigma_1 = mean(clim_sigs(1,:));                    % get silverman's sigmas
% %                   pdfs_u1 = pdfs_t;
%                     pdfs_t = obj.normalize_pdfs_by_sigmas(x, myclim_mus, myclim_sigs, pdfs_t, max(abs(myclim_sigs)), true);
%                                             nnegs = sum(pdfs_t(:)<-1e-18);
%                                             if ( nnegs > 0)
%                                                 fprintf("negatives:  %3d  min %.5g\n", nnegs, min(pdfs_t(:)));
%                                             end
% %                   pdfs_u2 = pdfs_t;
%                     lbl = "sigma-adjusted kde_surface";
%                 else
                if (obj.RP.median_normalize)
                    x = to_column(obj.RP.bins);
 %                  obj.pthresh_base = find_problines(obj.pdf_base, obj.CDF_base, pthresh_probs, obj.RP.edges, false, 1, obj.RP.cdf_thresh, false); % NOTE:  no need to filter using nterms & sig_terms;  we've already done that.                    
                
%                   pdf_medians = find_problines(pdfs_t, [], 0.5, obj.RP.edges, false, 1, obj.RP.cdf_thresh, false, 6, 2);
                    pdf_medians = find_problines(pdfs_t, [], 0.5, obj.RP.edges, false, 1, 1e-12, false, obj.RP.clim_nterms, obj.RP.clim_sig_terms);
%                   pdfs_t = obj.normalize_pdfs_means(x, myclim_mus, pdfs_t, true);
                    pdfs_t = obj.normalize_pdfs_medians(x, pdf_medians, pdfs_t, true);
                    lbl = "median-centered kde_surface";
                end
  
                    % see notes in ARRM_V2_DataParams on creating internal plots.
                if (obj.DP.intern_plotflag(j)), plot_surface(fignum, pdfs_t, 1:pdf_yrlen, bins, false, [3,3,7],lbl,false, 'time', ylbl, 'counts'); end
                
                if (obj.RP.sigma_normalize)
                    [my_norm_sigmas, ~] = remap_time_dimension(obj.norm_sigmas', obj.RP.pdf_yrlen, []);  % compress norm_sigmas the same as we have with the histogram table.
                    my_norm_sigmas = my_norm_sigmas*pdf_yrlen/obj.RP.yrlen;
                else
%                   my_norm_sigmas = ones(1,obj.RP.pdf_yrlen);
                    my_norm_sigmas = ones(1,size(pdfs_t,2));
                end
                
                        % next, map to probability space, and smooth.
                try
                [prb_noisy,~,prb_smooth] = probline_surf(pdfs_t,    bins,sdvals,1, 'fwd',     false, obj.RP.anom_nterms(1), obj.RP.anom_sig_terms(1), my_norm_sigmas);
                catch me
                    report_me_error(me);
                    fprintf("oops create_kde_pdfs_lpf probline_surf\n");
                    rethrow(me);
                end
                
                if (obj.DP.intern_plotflag(j))                     % see notes in ARRM_V2_DataParams on creating internal plots.
                    plot_surface(fignum,prb_noisy, 1:pdf_yrlen,sdvals,sdlines,[3,3,2], sprintf('PDF v STDev %s', obj.DA_title),     false, 'time','std devs',ylbl,[]);
                    plot_surface(fignum,prb_smooth, 1:pdf_yrlen,sdvals,sdlines,[3,3,5], 'smoothed v STDEV',false, 'time','std devs',ylbl,[]);
                    plot_surface(fignum,prb_smooth, 1:pdf_yrlen,sdprbs,sdprbl,[3,3,8], 'smoothed v Prob',false, 'time','probability',ylbl,[]);
                end
%               zzz
                    % added for sigma_normalizing.  Will need to update precip to use out_edges.
                    % if we keep this, we can do away with all the calc_binning stuff and just use sigma_normalizing to do the histogramming.
%                 if (obj.RP.sigma_normalize)
%                     min_anom = min(obj.anoms,[],"all","omitnan");
%                     max_anom = max(obj.anoms,[],"all","omitnan");
%                     myrange = max_anom - min_anom;
%                     outrange = [min_anom - .125*myrange, max_anom + .125*myrange];
%                     if (myrange < 50)
%                         out_edges = linspace(outrange(1),outrange(2),625);
%                     else
%                         out_edges = outrange(1):.1:outrange(2);
%                     end
%                     obj.RP.out_edges = out_edges;
%                     obj.RP.out_binstep = (outrange(2)-outrange(1))/(length(out_edges)-1);
%                 else
%                     obj.RP.out_edges = obj.RP.edges;
%                     obj.RP.out_binstep = obj.RP.binstep;
%                 end      
                
                    
                if (obj.RP.sigma_normalize)
                    if (isempty(obj.RP.out_edges))
                        out_min = min(prb_smooth,[], "all","omitnan");
                        out_max = max(prb_smooth,[], "all","omitnan");
                        out_lims = [floor(out_min-2), ceil(out_max+2)];
                        dx = min(0.1, range(out_lims)/500);
                        nbins = ceil(range(out_lims) - 1e-12)/dx;    % subtract fudge factor, then use ceil(...), to avoid computer math floating point errors.
%                         if (mod(npts,2)==0)
%                             npts=npts+1; 
%                         end
                        nedges = nbins+1;
                        obj.RP.out_edges = linspace(out_lims(1), out_lims(2), nedges);
                    end
                else
                    obj.RP.out_edges = obj.RP.edges;
                end      
            
%               if (isempty(obj.RP.prcp_distrib) || strlength(obj.RP.prcp_distrib)==0)
                if (isempty(obj.RP.prcp_distrib) || any(strcmp(obj.RP.prcp_distrib,["log","pwr"])))
                        % now map back into standard pdf format
                    nout_bins = length(obj.RP.out_bins);
    %               [pdfs_t,nanmap, ~]  = probline_surf(prb_smooth,            bins,sdvals,1, 'reverse', false, [], [], my_norm_sigmas);
                    [pdfs_t,nanmap, ~]  = probline_surf(prb_smooth, obj.RP.out_bins,sdvals,1, 'reverse', false, [], [], my_norm_sigmas);
                    nnegs = sum(pdfs_t(:)<0);
                    if ( nnegs > 0)
                        fprintf("negatives:  %3d  min %.5g\n", nnegs, min(pdfs_t(:)));
                    end
                    if (obj.DP.intern_plotflag(j)), plot_surface(fignum,pdfs_t,1:pdf_yrlen,obj.RP.out_bins,probs_flag,[3,3,3], lbl,false, 'time',ylbl,'Probability',[]); end
                else
                        % now map back into standard pdf format
                    nout_bins = length(obj.RP.out_bins);
                    [pdfs_t,nanmap, ~]  = probline_surf(prb_smooth,            bins,sdvals,1, 'reverse', false, [], [], my_norm_sigmas);
                    nnegs = sum(pdfs_t(:)<0);
                    if ( nnegs > 0)
                        fprintf("negatives:  %3d  min %.5g\n", nnegs, min(pdfs_t(:)));
                    end
                    if (obj.DP.intern_plotflag(j)), plot_surface(fignum,pdfs_t,1:pdf_yrlen,           bins,probs_flag,[3,3,3], lbl,false, 'time',ylbl,'Probability',[]); end
                end                
                
%               if (obj.DP.intern_plotflag(j)), plot_surface(fignum,pdfs_t,1:pdf_yrlen,    bins,probs_flag,[3,3,3], lbl,false, 'time',ylbl,'Probability',[]); end
%               if (obj.DP.intern_plotflag(j)), plot_surface(fignum,pdfs_t,1:pdf_yrlen,obj.RP.out_bins,probs_flag,[3,3,3], lbl,false, 'time',ylbl,'Probability',[]); end
                
                        % undo normalizing by sigmas.
%                 if (obj.RP.do_normalize_sigmas)
%                     pdfs_t(nanmap) = 0;
%                     pdfs_t = obj.normalize_pdfs_by_sigmas(x, myclim_mus, myclim_sigs, pdfs_t, max(abs(myclim_sigs)), false);
%                                             nnegs = sum(pdfs_t(:)<0);
%                                             if ( nnegs > 0)
%                                                 fprintf("negatives:  %3d  min %.5g\n", nnegs, min(pdfs_t(:)));
%                                             end
% %                   pdfs_t(nanmap) = nan;
%                     pdfs_t(pdfs_t == 0) = nan;
%                     lbl = 'Denormalized by Sigmas';
%                 else
                if (obj.RP.median_normalize)
                    pdfs_t(nanmap) = 0;
%                   pdfs_t = obj.normalize_pdfs_means(x, myclim_mus, pdfs_t, false);
                    pdfs_t = obj.normalize_pdfs_medians(obj.RP.out_bins, pdf_medians, pdfs_t, false);
                    pdfs_t(pdfs_t == 0) = nan;
                    lbl = "median reset";
                end
                if (obj.RP.adjust_pdfs)
 %                   if (obj.DP.intern_plotflag(j)), plot_surface(fignum,pdfs_t,1:pdf_yrlen,    bins,probs_flag,[3,3,6], sprintf('%s top view', lbl),false, 'time',ylbl,'Probability',[],[90,90]); end
                    if (obj.DP.intern_plotflag(j)), plot_surface(fignum,pdfs_t,1:pdf_yrlen,obj.RP.out_bins,probs_flag,[3,3,6], sprintf('%s top view', lbl),false, 'time',ylbl,'Probability',[],[90,90]); end
                
                    if (obj.RP.sigma_normalize)
                                           % this is simple if we sigma-normalized.  In that case, we used (almost)
                                            % identical silversigs, with the same # of values in each day, and all days
                                            % had a sigma of (approximately) 1.  So we'll simply reduce the bins by
                                            % sqrt(1+mean(silversigs).
                                            % see notes in ARRM_V2_DataParams on creating internal plots.
                        scale_factor = sqrt(1+(mean(silversigs))^2);
                        try
                             sigs = pdf_sigmas(pdfs_t, obj.RP.out_bins, [], false);        % sigs is widened from original sigmas by KDE smoothing.
                        catch
                            fprintf("oops\n");
                        end
                        sigs_orig = sigs ./scale_factor;
                        pdfs_t = adjust_kde_pdfs(pdfs_t, sigs_orig, obj.RP.out_bins);
                    else

                        try
    %                       pdfs_u = pdfs_t;
%                           if (~isempty(obj.RP.prcp_distrib) && strlength(obj.RP.prcp_distrib)>0)
                            if (obj.RP.isPrecipRun && any(strcmp(obj.RP.prcp_distrib,["generalizedpareto", "inversegaussian", "loglogistic", "lognormal"])))
                                pdfs_t = adjust_kde_pdfs(pdfs_t, clim_sigs_orig,     bins);
                            else
                                pdfs_t = adjust_kde_pdfs(pdfs_t, clim_sigs_orig, obj.RP.out_bins);
                            end
                            nnegs = sum(pdfs_t(:)<0);
                            if ( nnegs > 0)
                                fprintf("negatives:  %3d  min %.5g\n", nnegs, min(pdfs_t(:)));
                            end
                        catch
    %                       obj.DP.warn_log(2, string(oops()));
                            obj.DP.warn_log("create_kde_pdfs_lpf adjust_kde_pdfs oops\n");
                        end
                    end
                end                    % see notes in ARRM_V2_DataParams on creating internal plots.
                if (obj.DP.intern_plotflag(j))
                    if (~isempty(dayPos))
                        fin_pdf = remap_time_dimension(pdfs_t, obj.RP.yrlen, [], dayPos(j,:));
                    else
                        fin_pdf = pdfs_t;
                    end
                    plot_surface(fignum,fin_pdf,1:obj.RP.yrlen,           bins,probs_flag,[3,3,9], sprintf('final pdf %s', runLbl),false, 'time',ylbl,'Probability',[]); 
%                   plot_surface(fignum,fin_pdf,1:obj.RP.yrlen,obj.RP.out_bins,probs_flag,[3,3,9], sprintf('final pdf %s', runLbl),false, 'time',ylbl,'Probability',[], [90,90]); 
                end

                if (isempty(pdfs))
                    pdfs = zeros(nout_bins, pdf_yrlen, nsets);
                    if (nargout > 1)
                        CDFs = zeros(nout_bins, pdf_yrlen,nsets);
                    end
                end
                if (nsets > 1)
                    pdfs(:,:,j) = pdfs_t; %#ok<AGROW>
                else
                    pdfs(:,:) = pdfs_t;
                end

                    % Now integrate to create PDFs
                if (nargout > 1)
                    CDFs = zeros(size(pdfs));
                    try
                        CDFs(:,:,j) = cumsum(pdfs(:,:,j),'omitnan');
                        CDFs(:,:,j) = CDFs(:,:,j) ./ CDFs(end,:,j);   % normalize
                    catch
                        fprintf("oops integrating CDFs\n");
                    end
                end
            end

            pdfs = pdfs ./ sum(pdfs,"omitnan");   % normalize
            if (nargout > 1)
                CDFs = squeeze(CDFs);
            end
        end             % end of create_kde_pdfs_lpf
        
%         function [pdfs, CDFs] = create_kde_pdfs_norm_sigmas(obj, hist_tbl)
%         %   Not used.  Superceded by create_kde_pdfs_lpf.
%         %   This function currently not used.  It normalizes the pdfs by std deviation, smooths along the time dimension
%         %   with a gaussian, then denormalizes.  
%         %   Calculates smoothed pdfs, smoothed both along the day-of-year and the pdf  dimensions.
%         %   Calculates CDFs from pdfs only if nargout > 1.
% 
%             dx = obj.RP.bins(2) - obj.RP.bins(1);
%             
%             [nbins,yrlen,nsets] = size(hist_tbl);
% 
%                 % this creates a 3D matrix with dimensions (yvals, day-of-year, #years) , 
%                 % where #years is total years / yr_step -- or think of it as 1 2D matrix (nyvals x yrlen) for each 5 year period.
% 
%         %         % smooth histograms along day-of-year
%         %         
%         %     G = GAUSS_fft(yrlen, sigrange*yrlen/365, true);
%         %     htbl = ifft(G .* fft(hist_tbl,yrlen,2),yrlen,2);
% 
%                 % Now do KDE down the histogram (yval) dimension
% 
%                 % normalize histograms to sum to 1
% 
%             sumhist = sum(hist_tbl);
%             pdf_tbl = hist_tbl ./ sumhist;
%             sumhist = squeeze(sumhist);         % get rid of singleton 1st dimension.  (we needed it 1 x yrlen x npdfs for previous statement)
%             if (isrow(sumhist)), sumhist=sumhist'; end
%             x = obj.RP.bins;
%             if (isrow(x)), x=x'; end
% 
%             GS         = zeros(nbins,yrlen);
%             pdfs = zeros(size(pdf_tbl));
%             if (nargout > 1)
%                 CDFs = zeros(size(pdf_tbl));
%             end
%             
%             anom_nterms = obj.RP.anom_nterms(end); % or (1).  Added 2nd term for create_kde_pdfs_lpf.
%             anom_sig_terms = obj.RP.anom_sig_terms(end);  % or (1).  Added 2nd term for create_kde_pdfs_lpf.
%             
%             mus       = zeros(2,yrlen);
%             clim_mus  = zeros(2,yrlen);
%             sigs      = zeros(6,yrlen);
%             clim_sigs = zeros(6,yrlen);
%             if (any(strcmp(obj.RP.runType,'temp')))
%                 ylbl = 'Deg C';
%             elseif (obj.DP.isPrecipRun)
%                 ylbl = 'Precip, mm';
%             else
%                 ylbl = '';
%             end
% % z=0;
%             for j=1:nsets
%                 if (obj.DP.intern_plotflag(j)), fignum = 2*(j-1)+1+5*obj.DP.figbase; end                     % see notes in ARRM_V2_DataParams on creating internal plots.
%                 pdfs_t = squeeze(pdf_tbl(:,:,j));
%                 npts   = sumhist(:,j);
%         %                                     counts = nansum(y,1).*s(1,:,j);             % # of readings for each data
%         %                                     prob = y ./ s(:,j);                         % daily probabilities of each value               
%                 mus(1,:) = sum(x .* pdfs_t);                       % mean daily value, s/b around midpoint of bins.
% %               mu0 = mean(mus(1,:));                               % center point of distributions
% %               bin_ctr = interp1(obj.RP.bins, 1:nbins, mu0);                  % equivalent bin (fractional).
%                 
%                 sigs(1,:) = sqrt(sum(x.^2 .* pdfs_t) - mus(1,:).^2); % NOT now:  mus should be fairly continuous, since we smoothed already along the time axis
%                 clim_sigs(1,:) = obj.DA_climatology(sigs(1,:),anom_nterms,anom_sig_terms(end), yrlen)';   %<---- should be changed if we keep this, Ian.
%                 sigma_1 = mean(clim_sigs(1,:));                    % get silverman's sigmas
% %                 if (obj.RP.do_normalize_sigmas)
% %                     clim_sigs(1,:) = obj.DA_climatology(sigs(1,:),anom_nterms,anom_sig_terms, yrlen)';
% %                     sigma_1 = mean(clim_sigs(1,:));                             % we'll normalize all days to this, as "one standard deviation"
%                     clim_mus(1,:) = obj.DA_climatology(mus(1,:),anom_nterms,anom_sig_terms(end), yrlen)';
%                     
%                     if (obj.DP.intern_plotflag(j)), plot_surface(fignum, pdfs_t, 1:yrlen, x, false, [3,2,1],'raw histogram',true, 'time', ylbl, 'prob'); end                     % see notes in ARRM_V2_DataParams on creating internal plots.
%                     if (obj.RP.do_normalize_sigmas)                    
%                         pdfs_t = obj.normalize_pdfs_by_sigmas(x, clim_mus(1,:), clim_sigs(1,:), pdfs_t, max(abs(clim_sigs(1,:))), true);
%                     end
%                     mus(2,:) = sum(x .* pdfs_t);
%                     sigs(2,:) = sqrt(sum(x.^2 .* pdfs_t) - mus(2,:).^2); % NOT now:  mus should be fairly continuous, since we smoothed already along the time axis
%                     clim_sigs(2,:) = obj.DA_climatology(sigs(2,:),anom_nterms,anom_sig_terms, yrlen)';
%                     
%                     if (obj.DP.intern_plotflag(j)), plot_surface(fignum, pdfs_t, 1:yrlen, x, false, [3,2,3],'Normalizing sigmas',false, 'time', ylbl, 'prob'); end                    % see notes in ARRM_V2_DataParams on creating internal plots.
%                                         
%                     ss0 = 1.06 * sigma_1 .* (mean(npts)^(-.2))/dx;
%                     silversigs = repmat(ss0, 1, size(sigs,2));
%                     silversigs = max(silversigs, 1.5);              % If silversigs are too small, the filter ends up distorted because the gaussian in the fft doesn't get small enough.
%                         % make gaussians based on silversigs
%                     for i=1:yrlen, GS(:,i) = GAUSS_fft(nbins, silversigs(i), true)'; end
%                     pdfs_t = abs(ifft(GS .* fft(pdfs_t)));     % s/b purely real, but may have tiny imag. part, so take abs of inverse fft here.
%                     mus(3,:) = sum(x .* pdfs_t);
%                     sigs(3,:) = sqrt(sum(x.^2 .* pdfs_t) - mus(3,:).^2); % NOT now:  mus should be fairly continuous, since we smoothed already along the time axis
%                     clim_sigs(3,:) = obj.DA_climatology(sigs(3,:),anom_nterms,anom_sig_terms, yrlen)';
% 
%                     if (obj.DP.intern_plotflag(j)), plot_surface(fignum, pdfs_t, 1:yrlen, x, true, [3,2,5],'KDE smoothing',false, 'time', ylbl, 'prob'); end                    % see notes in ARRM_V2_DataParams on creating internal plots.
% 
%                     clim_mus(2,:) = obj.DA_climatology(mus(2,:), anom_nterms, anom_sig_terms, yrlen)';
%                  
% %                    pdfs_t = obj.smooth_along_time(yrlen, obj.RP.sigrange, pdfs_t, false, obj.RP.do_variable_gaussians, max(clim_sigs), .95);
% %                   pdfs_t = obj.smooth_along_time(yrlen, obj.RP.sigrange, pdfs_t, false, false, mean(clim_sigs(1,:)), .95);
%                                                                                                                                    %         95th, 99.9th
% %                   pdfs_t = obj.smooth_along_time(yrlen, obj.RP.sigrange, pdfs_t, false, obj.RP.do_variable_gaussians, bin_ctr, sigma_1/dx, 1.645,3.29, 3);
%  %                  pdfs_t = obj.smooth_along_time(yrlen, obj.RP.sigrange, pdfs_t, obj.sigrange_probs, obj.sig_scale, false);
%                     pdfs_t = obj.smooth_along_time(pdfs_t, yrlen, obj.RP.sigrange, obj.RP.sigrange_probs, obj.RP.sig_adj_linear);
%                     mus(4,:) = sum(x .* pdfs_t);
%                     
%                     if (obj.DP.intern_plotflag(j)), plot_surface(fignum, pdfs_t, 1:yrlen, x, true, [3,2,2],'time smoothing',false, 'time', ylbl, 'prob'); end                    % see notes in ARRM_V2_DataParams on creating internal plots.
% 
%                     sigs(4,:) = sqrt(sum(x.^2 .* pdfs_t) - mus(4,:).^2); 
%                     if (obj.RP.do_normalize_sigmas)
%                         pdfs_t = obj.normalize_pdfs_by_sigmas(x, clim_mus(1,:), clim_sigs(1,:), pdfs_t, max(abs(clim_sigs(1,:))), false);
%                     end
%                     mus(5,:) = sum(x .* pdfs_t);
%                     sigs(5,:) = sqrt(sum(x.^2 .* pdfs_t) - mus(5,:).^2); 
%                     
%                     if (obj.DP.intern_plotflag(j)), plot_surface(fignum, pdfs_t, 1:yrlen, x, true, [3,2,4],'remap back to original means & sigmas',false, 'time', ylbl, 'prob'); end                     % see notes in ARRM_V2_DataParams on creating internal plots.
% %                 else
% %                     % need to smooth along time first for this.
% %                    pdfs_t = obj.smooth_along_time(yrlen, obj.RP.sigrange, pdfs_t, false, obj.RP.do_variable_gaussians, bin_ctr, sigma_1/dx, 1.645,3.29, 3);
% %                     silversigs = 1.06 * sigs(1,:) .* (ss.^(-.2))/dx;
% %                     silversigs = obj.DA_climatology(silversigs,yrlen, 4,2);
% %                     silversigs = max(silversigs, 1.5);              % If silversigs are too small, the filter ends up distorted because the gaussian in the fft doesn't get small enough.
% %                         % make gaussians based on silversigs
% %                     for i=1:yrlen, GS(:,i) = GAUSS_fft(nbins, silversigs(i), true)'; end
% %                     pdfs_t = abs(ifft(GS .* fft(pdfs_t)));     % s/b purely real, but may have tiny imag. part, so take abs of inverse fft here.
% % %                    pdfs_t = obj.smooth_along_time(yrlen, obj.RP.sigrange, pdfs_t, false, obj.RP.do_variable_gaussians, bin_ctr, sigma_1/dx, 1.645,3.29, 5);
% %                 end
% 
%                 ss=sum(pdfs_t);
%                 pdfs_t = pdfs_t./ss;
% %                plot_surface(98, pdfs_t, 1:yrlen, x, true);
%                 
%                 if (obj.RP.adjust_pdfs)
%                         % now adjust for spreading from convolving w/ silversigs.
%                     pdfs_adj = zeros(nbins, yrlen);
%                     ssigs = silversigs*dx;
%                     clim_ssigs = obj.DA_climatology(ssigs, anom_nterms, anom_sig_terms, yrlen)';
%                     
%                     x_adj  = zeros(size(pdfs_t));
% %                   x_adj2 = zeros(size(pdfs_t));
% 
%                     for i=1:yrlen
%                 %        out_yvals_adj = mus(i) + (out_yvals-mus(i)) * sigs(i)/sqrt(sigs(i)^2+ssigs(i)^2);        % original equation...
%                 %        out_yvals_adj = (out_yvals) * sigs(i)/sqrt(sigs(i)^2+ssigs(i)^2);   % kludge to try tuning things
%                        x_adj(:,i)  = clim_mus(1,i) + (x-clim_mus(1,i)) * clim_sigs(2,i)/sqrt(clim_sigs(2,i)^2+(clim_ssigs(i))^2);
% %                      x_adj2(:,i) =      mus(1,i) + (x-     mus(1,i)) *      sigs(1,i)/sqrt(     sigs(1,i)^2+(     ssigs(i))^2);
%                         try
%                             pdfs_adj(:,i) = interp1(x_adj(:,i), pdfs_t(:,i), x,'pchip',0);
%                         catch me
%                             obj.DP.warn_log(oops("oops! jc_kdm:  problem adjusting for convolving gaussian with gaussian\n"));
%                             msgtext=getReport(me);
%                             obj.DP.warn_log('%s\n', msgtext);
%                             obj.DP.warn_log('------\n');        
%                         end
%                     end  
%                          % normalize adjusted pdfs to make them pdfs.    
%                     ss = sum(pdfs_adj);
%                     pdfs_t = pdfs_adj ./ ss;
%                     if (sum(isnan(pdfs_t(:)))>0)
%                         obj.DP.warn_log(oops('oops! jckdm:  pdf adj produced nans\n'));
%                     end
%                 end
%                 mus(6,:) = sum(x .* pdfs_t);
%                 sigs(6,:) = sqrt(sum(x.^2 .* pdfs_t) - mus(6,:).^2); 
%                 
%                 if (obj.DP.intern_plotflag(j))                     % see notes in ARRM_V2_DataParams on creating internal plots.
%                     plot_surface(fignum, pdfs_t, 1:yrlen, x, true, [3,2,6],'Adjusting for KDE spreading',false, 'time', ylbl, 'prob');             
%                     obj.plot_sigmas(fignum+1, mus, clim_mus,  sigs, clim_sigs);
%                 end
%                 
%                 pdfs(:,:,j) = pdfs_t;
%                     % Now integrate to create PDFs
%                 if (nargout > 1)
%                     CDFs(:,:,j) = cumsum(pdfs(:,:,j),'omitnan');
%                 end
%             end
% 
%             if (nsets == 1)
%                 pdfs = squeeze(pdfs);
%                 if (nargout > 1)
%                     CDFs = squeeze(CDFs);
%                 end
%             end
%         end
        
%       function [wd_clim, wd_phase] = calc_precip_climatology(obj, fldname, pdf_yrs)        % this should do something with NAs to normalize for missing data...
%       function                       calc_precip_climatology(obj, fldname, pdf_yrs)        % this should do something with NAs to normalize for missing data...
        function                       calc_precip_climatology(obj, fldname)                 % this should do something with NAs to normalize for missing data...
            
            prcp_min = obj.DP.prcp_min;
            yrlen = obj.DP.yrlen;
                
            prcp = obj.(fldname);
            
            prcp = reshape(prcp, yrlen,[]);
%           nyrs = size(prcp,2);
%             if (~exist('pdf_yrs','var') || isempty(pdf_yrs))
%                 pdf_yrs = nyrs;
%             end
            
            prcp_counts = sum(prcp>=prcp_min,2);
            
            [wd_clim, wd_avg, wd_phase, wd_clim_raw] = obj.DA_climatology(prcp_counts, 12, 4, yrlen);
                % transpose, because we need it as a row later.
            wd_clim = wd_clim';
            wd_clim_raw = wd_clim_raw';

            obj.wet_day_clim  = wd_clim; 
            obj.wet_day_avg   = wd_avg; 
            obj.wet_day_phase = wd_phase;
            obj.wet_day_clim_raw = wd_clim_raw;
% [clim, avg, phase, clim_raw] = DA_climatology(obj, y, nterms, sig_terms, yrlen)     
                        % climatology for each rolling step
            if (isempty(obj.DP.rolling_yrs))
                nsets = 1;
            else
                nsets = length(obj.DP.rolling_steps);
                obj.wet_day_rolling_clim = zeros(nsets,obj.RP.yrlen);
                obj.wet_day_rolling_clim_raw = zeros(nsets,obj.RP.yrlen);
                prcp_counts = zeros(yrlen, nsets);
                for j=1:nsets
                    [rstart, rend] = obj.using_range('rolling_yrs', j);
%                   try
                    p=reshape(prcp(rstart:rend), yrlen,[]);
                    prcp_counts(:,j) = sum(p>=prcp_min,2);

                        [wdrc,~,~, wdrc_raw] = obj.DA_climatology(prcp_counts(:,j), 12, 4, yrlen); 
                        obj.wet_day_rolling_clim(j,:) = wdrc';
                        obj.wet_day_rolling_clim_raw(j,:) = wdrc_raw';
 %                  catch
 %                      fprintf("oops wet day rolling clim\n");
 %                  end
                end
            end                        
                        
            
            
 %          tot_prcp_events = sum(prcp_counts);
            
                    % update the pdf_yrlen, based on the total number of precip events.
                    
% tot_precip here is smaller of total obs precip and total model precip. 
% per statisticians, we need a minimum of 57 events to make reasonable estimates, so we make sure we get at least 57/day.
% Here, we're using a yearlength to give us at least 75 events/day when we use a filter width of one month.
%
%       IAN:  check this for doing the filtering in z-space!
            pdf_yrlens = zeros(nsets,1);
            for j=1:nsets
                tot_prcp_events = sum(prcp_counts(:,j));
                pdf_yrlens(j) = min(obj.RP.pdf_yrlen, floor(tot_prcp_events/65));      % aim for at least 65 precip events in each compressed day  (stats folks say 57.  I'm being conservative here.)
                if (pdf_yrlens(j) < 4)
                    pdf_yrlens(j) = 1;       % if so few precips that we can't do every 3 months, collapse it to a single pdf
                end
            end
            
            obj.RP.pdf_yrlen = min(pdf_yrlens);         % would be better to have unique value for each pdf_yrstep.  But thiswould only be for very dry locations, 
                                                            % and would complicate a lot of things in the code to address those few locations.
                
            obj.RP.anom_nterms = min(obj.RP.anom_nterms, floor(obj.RP.pdf_yrlen/2));      %  shrink # cycles for filtering down if yrlen too short (<24).  2 cycles only if yrlen == 8.
            obj.RP.clim_nterms = min(obj.RP.clim_nterms, floor(obj.RP.pdf_yrlen/2));      %  shrink # cycles for filtering down if yrlen too short (<24).  2 cycles only if yrlen == 8.
            
            obj.RP.anom_sig_terms = min(obj.RP.anom_sig_terms, obj.RP.anom_nterms/2);
            obj.RP.clim_sig_terms = min(obj.RP.clim_sig_terms, obj.RP.clim_nterms/2);            
        end
        
        function [clim, avg, phase, clim_raw] = DA_climatology(obj, y, nterms, sig_terms, yrlen)        
        % clim = climatology(y, nterms, sig_term, yrlen) 
        %
        % calculates the average daily value (365x1) for the data in y  (365*nyrs x 1)  (or yrlen x 1 and yrlen*nyrs x 1)
        % then low-pass filters it with a filter defined by nterms and sig_term.
        %
        %   Note that the average here is a straight average for each day of the year (not gaussian-weighted), which is
        %   then smoothed by circular convolution along the 365 days.
        %
        % Returns a 365x1 vector representing the smoothed daily average of the data in y.
        %   Inputs:
        %       y           data to work with.  must be of multiple of 365 days long
        %       nterms      equivalent # of frequency terms to retain
        %       sig_term    sigma for gaussian to smooth rectangular filter.
        %                       good values to use for nterms & sig_term are 6 & 2.0
        %                       this gives the equivalent filtering of an ideal
        %                       filter keeping 1st 6 terms of fft (up to 6 cycles
        %                       per year)
        %                       for pure gaussian filtering, use nterms=0, sig_term set to desired gaussian sigma.  
        %       yrlen       length of year for data (usually 365).
        %frequency domain sigma to smooth with.  See "math notes" for equations to go from
        %                   either time-domain gaussian sigma or idea-filter-equivalent sigma.

            if (~exist('yrlen',    'var') || isempty(yrlen)),        yrlen  = obj.RP.yrlen;          end
            if (~exist('nterms',   'var') || isempty(nterms)),       nterms = obj.RP.clim_nterms;    end
            if (~exist('sig_terms','var') || isempty(sig_terms)), sig_terms = obj.RP.clim_sig_terms; end

            [clim, ~, phase, avg, clim_raw] = climatology(y, nterms, sig_terms, yrlen);
        end

        function [clim, avg, phase, clim_raw] = old_DA_climatology(obj, y, nterms, sig_terms, yrlen)    % previous version, before reconciling the external climatology function with this one.
            
            if (~exist('yrlen',    'var') || isempty(yrlen)),        yrlen  = obj.RP.yrlen;          end
            if (~exist('nterms',   'var') || isempty(nterms)),       nterms = obj.RP.clim_nterms;    end
            if (~exist('sig_terms','var') || isempty(sig_terms)), sig_terms = obj.RP.clim_sig_terms; end

            if (nargout > 1), avg = mean(y); end
            [nr,nc] = size(y);
            if (nr == 1 || nc == 1)        
                    % if keeping all terms, just do straight average:
                nyrs = length(y)/yrlen;    
                clim = mean(reshape(y,yrlen,nyrs),2);     % get average daily value for data range.
            else
                clim = mean(y,2);
            end
            clim_raw = clim;

                % Create the filter, and low-pass filter the daily averages
            FILT = calc_filter(nterms, sig_terms, 1, yrlen);
            clim = lpf_FILT(clim, FILT);
            
            if (nargout > 2)
                phase = calc_phase(clim, nterms, 10, yrlen);
            end
            
        end
        
        function [clim, clim_1D, phase] = climatology_3D(obj, data)     % Ian:  this should be compared with the main climatology(...) function to make sure we are consistent.  icsf 9/25/23
        %
        % function to calculate moving (dynamic) climatology.  Long term trend should be removed first.
        %   Creates 3D climatology surface from low-pass gaussian filter of specified data.
        %
        %   smoothes climatology with gaussian filter w/ same std deviation as an ideal filter of length nterms (for 365-day year).
        %   Returns a 2-D climatology surface which is the smoothed overall surface.  
        %
        %   Long term phase shift in climatology is not calculated or removed in this version.  
        %   Code for removing phase shift is included, but commented out.
        %
        %       Inputs
        %
        %           data            input vector of temperature values, 1 per day, for period of interest, usually 1950 - 2100.
        %           RP              RunParams, containing:
        %               order           polynomial order to use (3 is good for 150-year data.  2nd order, parabolic, doesn't fit data well enough.
        %                                           (1 is good for 40- or 50-year historical period)
        %                                            0 for straight average
        %                                           '-1 :  data already detrended, or do-not-detrend.
        %               trend_yr_flags  boolean flag, 1 per year, flagging which years to include in trend calculation
        %                           (years with large fraction of NAs should be excluded to avoid biasing trend)
        %                           set to empty ( [] ) to use all years.
        %               yrlen           length of year, in days (usually 365).  must be integral value.
        %           nclim_years     # of years to use for rolling average in calculating climatology.
        %                               note:    will actually use a gaussian with the same std dev. as nclim_years.
        %                               note 2:  nclim_years does not need to be odd, but must be, integral value, such as 21.
        %                               for sigma-equiv. of n-year rolling average, sigma=sqrt((n^2-1)/12).  See MATH notes below
        %           nterms          equivalent # of fft terms to keep when filtering
        %                               larger values pass more 'noise'.  nterms=6 is equivalent to keeping all events of duration 2 months or longer.  (12/2)
        %                               6 is good if sig_terms=0 (ideal filter)
        %                               5 is good with sig_terms set to 2.0.
        %           sig_terms       sigma for gaussian to LPF nterms with.  
        %                               0 = ideal lowpass filter
        %                               2.0 = good value to use with nterms=5. This gives equivalent filtering to ideal  filter of nterms=6.
        %           n_ext_yrs       data will be extended to length of 256 years using average of 1st n_ext_yrs and last n_ext_yrs.
        %            DP              DataParams, containing:
        %               data_yrs        start, end years for entire data set
        %               clim_yrs       start, end years to use for climatology calculations
        %
        %       Outputs
        %
        %           clim            long-term climatology, as linear sequence
        %                               smoothed as a 2-D surface, first along the day, then smoothed along the year, then
        %                               reshaped back to a linear sequence.  To reshape into surface:  reshape(clim,365, nyrs);
        %                               Does NOT have the base climatology removed.
        %           clim_1D         1-D climatology, smoothed only along time, not along year as well.
        %
        % Data is extended to an even power-of-two years, then after all smoothing, results clipped back to original years.
        % FFT of even power of 2 is faster than fft of 151 or 191
        %
        %       MATH Notes:
        % notes re sigmas, gaussians & fft's of gaussians:
        %   1.  sigma of uniform distribution of length n: 
        %           for continuous:  n/2/sqrt(3) = n/sqrt(12)
        %           for discrete:  sqrt((n^2-1)/12) , 
        %           which, as m -> large, -> n/sqrt(12).
        %   1.a. n, length of uniform dist for given sigma:  
        %           n = sqrt(12*(sigma^2) + 1) 
        %   2.  sigma-0 of gaussian s.t. sigma(fft(gaussian)), in freq = sigma in time (for n points): sqrt(n/(2*pi))
        %   3.  sigma (freq. domain) for a sigma=1 (time domain) is .1592 = 1/(2*pi)
        %   4.  sigma-fft for arbitrary time domain sigma:  sigma-fft = sigma0^2 / sigma(time)
        %   4.a.    and same for given sigma-fft:  sigma(time domain) = sigma0^2 / sigma(freq domain)
        %   
        %__________________________________________________________________________________________________
        
                    % Make sure our data looks OK first.
            data = data(:);  
            npts = length(data);
            nyrs = npts/obj.RP.yrlen;

                    % Extend the data either end to a power-of-2 length years so we have something to convolve with on the ends.
                    % This appends the first few years' mean climatology at the beginning, and last few years' mean climatology
                    % at the end.
                    % Because of the efficiency of the power-of-2 FFT, it is faster to extend the series to a power of 2 years than
                    % just adding 21 years either end to get 193 years total.

            len=16;                
            nadd = len - nyrs;        
            while (len - 2*obj.RP.n_ext_yrs < nyrs)    % find smallest power of 2 long enough for extended data.
                len = len*2;
                nadd = len - nyrs;
            end
            nadd_start = floor(nadd/2);     % # of years to extend at start
            nadd_end   = nadd - nadd_start; % # of years to extend at end.

            y_ext = ARRM_V2_DisaggregateSignal.extend_data(data, obj.RP.yrlen, obj.RP.n_ext_yrs, nadd_start, nadd_end);     % data, extended so circular convolution of FFT doesn't affect opposite ends of data.

            mpts = length(y_ext);
            myrs = mpts/obj.RP.yrlen;            % this should be a power of 2)

                    % now, low-pass filter entire series with gaussian along the days,
                    % then reshape to a 2-D surface, and low-pass filter along the years dimension.

                            % get Fourier domain filter for extended sequence for smoothing along the days. 
            FILT_long = calc_filter(obj.RP.clim_nterms, obj.RP.clim_sig_terms, myrs, obj.RP.yrlen);

                    % convolve original (extended) signal along the days w/ gaussian based on nterms
                    % we do the convolution here in the fourier domain.  This is OK because we've extended the
                    % signal at either end so there's no problem with the circular effect of fft.
                    % (as long as RP.n_ext_yrs*365 >> sigt_long.)

            y2 = lpf_FILT(y_ext, FILT_long);
                        % strip off extensions to get 1-D moving climatology (smoothed only along time)
            ix1 = obj.RP.yrlen * nadd_start + 1;
            ix2 = ix1 + (obj.RP.yrlen * nyrs) - 1;
            clim_1D = y2(ix1:ix2);

                    % now, reshape and lpf along the years.
                    % See Math Notes above.
            ysurf = reshape(y2, obj.RP.yrlen, myrs);
            year_sig = sqrt((obj.RP.nclim_yrs^2-1)/12);
            clim_ext = ARRM_V2_DisaggregateSignal.lpf_surf(ysurf, year_sig);

        %---------------------------------------
        %           For now:  skip removing shifts.  We can add this in later if we find it is helpful.
        %           But removing shift is computationally expensive, and probably gains very little.

        %         % calculate and remove any phase shift in the signal.
        %     shifts = calc_shifts(base_clim, clim_ext, year_sig, day_sig);   %smooths by yrsig & daysig, then calculate shifts. 
        %          
        %     clim_ext = remove_shift(clim_ext_raw, shifts);
        %             
                    % remove the data extension from start and end of shifts.
        %     shifts = shifts(RP.n_ext_yrs+(1:RP.nyrs));
        %---------------------------------------

                    % remove the data extension from start and end
            clim  = clim_ext(:,nadd_start+(1:nyrs));
            
            if (nargout > 2)
                phase = nan(nyrs,1);
                mycos = cos(2*pi*(0:obj.RP.yrlen-1)/obj.RP.yrlen);
                for i=1:nyrs
                    y = mycos * range(clim(:,i)/2) +  mean(clim(:,i));
                    phase(i) = calc_phase(clim(:,i), obj.RP.clim_nterms, 10, obj.RP.yrlen, y);
%                     phase(i) = ic_corr_1d(mycos, clim(:,i)-mean(clim(:,i))/range(clim(:,i))*2, nterms, 10);
%                     phase(i) = mod(phase(i), obj.RP.yrlen);
                end
            end
                                        
            clim=reshape(clim,npts,1);

        end
        
        
                % NOTE:  debugging internal plots for appending probability shapes is turned on by set_internal_plotflag(32)
                % or by passing in "internal_plotflag",[0,32] in initialization
                % see notes in ARRM_V2_DataParams on creating internal plots.
                
        function [pdfs, cdfs, cdf_okflags] = append_basic_cdf(obj, pdfs, cdfs, pdf_mus) %, sigmas)   % update this, Ian!
            %this version reads empirical CDF info from a file and uses it.
            
            if (obj.DP.isStationRun)
                model = "Stations";
            else
                model = obj.DP.model;
            end
            
            if (isrow(pdf_mus)), pdf_mus = pdf_mus'; end    % 1-D may be a row instead of matrix w/ days going down 1st dimension...
            
            cdf_append_pts = obj.DP.cdf_append_pts;
            if (length(cdf_append_pts) == 1)
                cdf_append_pts(1) = -abs(cdf_append_pts(1));
                cdf_append_pts(2) =  abs(cdf_append_pts(1));
            end
            if (isrow(cdf_append_pts)), cdf_append_pts(2,:) = [-4.75,4.75]; end
            
            cdf_fname = fullfile(obj.DP.ARRM_V2_dir,"basic_cdfs",sprintf('cdf_basic.%s.%s.mat', model, obj.DP.varname));
            if (~exist(cdf_fname,'file'))
                obj.DP.warn_log("warning:  basic CDF file %s does not exist\n", cdf_fname);
                return;
            end
            load(cdf_fname, "mdl_bins","mdl_cdf","mdl_stats","okflags")
            [~, yrlen, nsets] = size(pdfs);
            my_okflags = pdfs ~= 0;
            
            cdf_okflags = false(size(cdfs));
            
            if (obj.DP.intern_plotflag(32))                    % see notes in ARRM_V2_DataParams on creating internal plots.
                cdfs1 = cdfs;
                pdfs1=pdfs;
            end
            
            for j=1:nsets
                for i=1:yrlen
%                   [out_cdf, out_okflags, anchor_probs, out_anchors, out_anchors_ix] = get_cdf_probs(cdf_append_pts, true, mdl_cdf, mdl_bins, okflags, mdl_stats.mu, cdfs(:,i,j), obj.RP.bins, my_okflags(:,i,j), pdf_mus(i,j));
                    [out_cdf, out_okflags,            ~,           ~, out_anchors_ix] = get_cdf_probs(cdf_append_pts, true, mdl_cdf, mdl_bins, okflags, mdl_stats.mu, cdfs(:,i,j), obj.RP.bins, my_okflags(:,i,j), pdf_mus(i,j));

                    ix1 = out_anchors_ix(1);
                    ix2 = out_anchors_ix(2);
                    cdfs(1:ix1,   i, j) = max(cdfs(1:ix1, i, j),   out_cdf(1:ix1),'omitnan');
                    cdfs(ix2:end, i, j) = min(cdfs(ix2:end, i, j), out_cdf(ix2:end),'omitnan');
                    cdf_okflags(:,i,j)  = out_okflags; 
                end
            end
            
            pdfs = diff([zeros(1,yrlen,nsets); cdfs],1,1);
%             pdfs1=pdfs;
%             cdfs1=cdfs;
%             myflag=false;
            if (any(pdfs(:)<0))     % because we interpolated points in get_cdf_probs, we get an occasional negative  from rounding issues.
                minp = min(pdfs(:));
                if (minp < -1e-4) 
                    obj.DP.warn_log('warning:  append_basic_cdf: max negative in CDF after get_cdf_probs: %9.6e\n', minp); 
%                     myflag=true;
                end
                pdfs(pdfs < 0) = 1e-16;
            end
            cdfs = cumsum(pdfs);        % and because we've mucked around with the values, we have to make sure we have true CDFs when we're done.
            cdfs = cdfs ./ cdfs(end,:,:);
            pdfs = diff([zeros(1,yrlen,nsets); cdfs],1,1);

%             if (myflag)
%                 save("zzz.mat","pdfs","pdfs1","cdfs","cdfs1");
%                 fprintf("bailing!\n");
%             end
 
            if (obj.DP.intern_plotflag(32))                    % see notes in ARRM_V2_DataParams on creating internal plots.
                for j=1:nsets
                    figure(99);
                    clf;
                    subplot(1,3,1);
                    surf(1:365, obj.RP.bins, pdfs(:,:,j),'facecolor','g','edgecolor','none'); 
                    hold on;  
                    surf(1:365, obj.RP.bins, pdfs1(:,:,j),'facecolor','red','edgecolor','none'); 
                    hold off;            
                    set(gca,'zscale','log');
                    lights_on(true);

                    subplot(1,3,2);
                    surf(1:365, obj.RP.bins, cdfs(:,:,j),'facecolor','g','edgecolor','none'); 
                    hold on;  
                    surf(1:365, obj.RP.bins, cdfs1(:,:,j),'facecolor','red','edgecolor','none'); 
                    hold off;            
                    set(gca,'zscale','log');
                    lights_on(true);

                    subplot(1,3,3);
                    surf(1:365, obj.RP.bins, 1-cdfs(:,:,j),'facecolor','g','edgecolor','none'); 
                    hold on;  
                    surf(1:365, obj.RP.bins, 1-cdfs1(:,:,j),'facecolor','red','edgecolor','none'); 
                    hold off;            
                    set(gca,'zscale','log');
                    lights_on(true);
                    view([-135,35]);
                end
            end

        end
        
    
                % NOTE:  debugging internal plots for appending probability shapes is turned on by set_internal_plotflag(32)
                % or by passing in "internal_plotflag",[0,32] in initialization
                % see notes in ARRM_V2_DataParams on creating internal plots.
                
        function [pdfs, cdfs, cdf_okflags] = append_sk_normal_cdf(obj, pdfs, cdfs, pdf_mus) %, pdf_sigmas)
            % This version assumes a scaled normal distribution in the tails, using the
            % specified anchor point (usual 2.5 or 3 sigma-equivalent point.)
            %
            %   Replaces pdf beyond the 1st anchor point with normal distribution whose sigma is calculated from the
            %   anchor point.
            %   Sets okflags for all points inside 2nd anchor points, although it still provides a CDF out to edge of
            %   the bin range.
            
                    % load the skew-generalized_normal table, and find the coefficients for the specific model.
            sk = load("model_pdf_coefs.mat");     % table with best sk_normal coefs for each model/var.
            
            if (obj.DP.isStationRun)
                model = "Stations";
            else
                model = obj.DP.model;
            end
            
                    % find the model in the list
            m = find(strcmp(sk.models,model), 1);
            v = find(strcmp(sk.vars, obj.DP.varname), 1);
                    % if we can't find it, just use the normal_cdf.
            if (isempty(m) || isempty(v))
                [pdfs, cdfs, cdf_okflags] = append_normal_cdf(obj, pdfs, cdfs, pdf_mus);
                return;
            end
            
            if (isrow(pdf_mus)), pdf_mus = pdf_mus'; end    % 1-D may be a row instead of matrix w/ days going down 1st dimension...
            
                    % make sure we've got 4 anchor points.
            anchor_pts = obj.DP.cdf_append_pts;
            if (length(anchor_pts) == 1)
                anchor_pts(1) = -abs(anchor_pts(1));
                anchor_pts(2) =  abs(anchor_pts(1));
            end
            if (isrow(anchor_pts)), anchor_pts(2,:) = [-4.75,4.75]; end
            [~, yrlen, nsets] = size(cdfs);
            okflags = pdfs ~= 0;
            
            cdf_okflags = false(size(cdfs));
            
            if (obj.DP.intern_plotflag(32))                % see notes in ARRM_V2_DataParams on creating internal plots.
                cdfs1 = cdfs;
                pdfs1=pdfs;
            end
            
            for j=1:nsets
                ab_bins = -10:.025:10;
                [ab_pdf, ab_cdf] = calc_model_pdfs(sk.alpha(m,v), sk.beta(m,v), 1, 0, ab_bins);    % check that mean is zero, Ian!
                ab_mean = pdf_moment(ab_pdf', 1, ab_bins,[], false);
                if (abs(ab_mean) >1e-6)
                    obj.DP.warn_log('ab_mean is not zero: %.6e\n', ab_mean);
                end
                ab_pdf = diff([0;ab_cdf]);  % recalculate pdf, because one returned about may not be quite right due to rounding near 1.
                ab_flags = ab_pdf > 1e-16;

                for i=1:yrlen
%                   [out_cdf, out_okflags, anchor_probs, out_anchors, out_anchors_ix] = get_cdf_probs(cdf_append_pts, true, mdl_cdf, mdl_bins, okflags, mdl_stats.mu, cdfs(:,i,j), obj.RP.bins, my_okflags(:,i,j), pdf_mus(i,j));
                    [out_cdf, ~          , anchor_probs,           ~, out_anchors_ix] = get_cdf_probs(anchor_pts, false, ab_cdf,  ab_bins,  ab_flags, ab_mean,  cdfs(:,i,j), obj.RP.bins,   okflags(:,i,j), pdf_mus(i,j));

                    ix1 = out_anchors_ix(1);
                    ix2 = out_anchors_ix(2);
                    cdfs(1:ix1,   i, j) = max(cdfs(1:ix1, i, j),   out_cdf(1:ix1),'omitnan');
                    cdfs(ix2:end, i, j) = min(cdfs(ix2:end, i, j), out_cdf(ix2:end),'omitnan');
                    ix3 = find(cdfs(:,i,j) >= anchor_probs(2),1);           % find first bin inside left-tail 2nd anchor pt
                    ix4 = find(cdfs(:,i,j) <= anchor_probs(4),1,'last');    % find first bin inside right-tail 2nd anchor point
                    cdf_okflags(ix3:ix4,i,j)  = true;                       % flag all points in between as OK. 
                end
            end
            pdfs = diff([zeros(1,yrlen,nsets); cdfs],1,1);
%             pdfs1=pdfs;
%             cdfs1=cdfs;
%             myflag=false;
            if (any(pdfs(:)<0))     % because we interpolated points using pchip, rather than linear, we get an occasional negative  from rounding issues.
                minp = min(pdfs(:));
                if (minp < -1e-4)
                    obj.DP.warn_log( 'warning:  append_sk_normal_pdf: max negative in CDF after get_cdf_probs: %9.6e\n', minp); 
%                     myflag=true;
                end
                pdfs(pdfs < 0) = 1e-16;
            end
            cdfs = cumsum(pdfs);        % and because we've mucked around with the values, we have to make sure we have true CDFs when we're done.
            cdfs = cdfs ./ cdfs(end,:,:);
            pdfs = diff([zeros(1,yrlen,nsets); cdfs],1,1);
%             if (myflag)
%                 save("zzz.mat","pdfs","pdfs1","cdfs","cdfs1");
%                 fprintf("bailing!\n");
%             end
 
            if (obj.DP.intern_plotflag(32))                % see notes in ARRM_V2_DataParams on creating internal plots.
                for j=1:nsets
                    figure(99);
                    clf;
                    subplot(1,3,1);
                    surf(1:365, obj.RP.bins, pdfs(:,:,j),'facecolor','g','edgecolor','none'); 
                    hold on;  
                    surf(1:365, obj.RP.bins, pdfs1(:,:,j),'facecolor','red','edgecolor','none'); 
                    hold off;            
                    set(gca,'zscale','log');
                    lights_on(true);

                    subplot(1,3,2);
                    surf(1:365, obj.RP.bins, cdfs(:,:,j),'facecolor','g','edgecolor','none'); 
                    hold on;  
                    surf(1:365, obj.RP.bins, cdfs1(:,:,j),'facecolor','red','edgecolor','none'); 
                    hold off;            
                    set(gca,'zscale','log');
                    lights_on(true);

                    subplot(1,3,3);
                    surf(1:365, obj.RP.bins, 1-cdfs(:,:,j),'facecolor','g','edgecolor','none'); 
                    hold on;  
                    surf(1:365, obj.RP.bins, 1-cdfs1(:,:,j),'facecolor','red','edgecolor','none'); 
                    hold off;            
                    set(gca,'zscale','log');
                    lights_on(true);
                    view([-135,35]);
                end
            end
        end

                % NOTE:  debugging internal plots for appending probability shapes is turned on by set_internal_plotflag(32)
                % or by passing in "internal_plotflag",[0,32] in initialization
                % see notes in ARRM_V2_DataParams on creating internal plots.
                
        function [pdfs, cdfs, cdf_okflags] = append_normal_cdf(obj, pdfs, cdfs, pdf_mus)
            % This version assumes a scaled normal distribution in the tails, using the
            % specified anchor point (usual 2.5 or 3 sigma-equivalent point.)
            %
            %   Replaces pdf beyond the 1st anchor point with normal distribution whose sigma is calculated from the
            %   anchor point.
            %   Sets okflags for all points inside 2nd anchor points, although it still provides a CDF out to edge of
            %   the bin range.
            
            if (isrow(pdf_mus)), pdf_mus = pdf_mus'; end    % 1-D may be a row instead of matrix w/ days going down 1st dimension...
            
            anchor_pts = obj.DP.cdf_append_pts;         % std.dev's equivalents (sigma-equivalents) for where to start replacing pdf w/ scaled normal.
            if (length(anchor_pts) == 1)
                anchor_pts(1) = -abs(anchor_pts(1));
                anchor_pts(2) =  abs(anchor_pts(1));
            end
            if (isrow(anchor_pts)), anchor_pts(2,:) = [-4.75,4.75]; end
            
            [~, yrlen, nsets] = size(cdfs);
            okflags = pdfs ~= 0;
            
            cdf_okflags = false(size(cdfs));
            
            if (obj.DP.intern_plotflag(32))                % see notes in ARRM_V2_DataParams on creating internal plots.
                cdfs1 = cdfs;
                pdfs1=pdfs;
            end
            
            anchor_probs = normcdf(anchor_pts(:));      % probabilities at anchor points.
            
            for j=1:nsets
                for i=1:yrlen
                    try
                    bins = obj.RP.out_bins - pdf_mus(i,j);    % shift bins so they represent the distance from the mean, not the zero-point.
                    catch
                        obj.DP.warn_log(oops());
                    end
                    flags = okflags(:,i,j);

                   try
                    anchors = interp1(cdfs(flags,i,j), bins(flags), anchor_probs, 'pchip');        % adjust bins to represent distance from mean for anchor points.
                   catch
                       obj.DP.warn_log(oops("append_normal_cdfs\n"));
                   end
                    ix1 = find(bins <  anchors(1),1,'last');    % ix of bin right before left anchor pt.
                    ix2 = find(bins >  anchors(3),1);           % ix of bin right after right-tail anchor pt
                    sig1_neg = anchors(1,1)/anchor_pts(1,1);    % sigma-equiv of left-tail anchor pt, distance/sigmas
                    sig1_pos = anchors(3)/anchor_pts(1,2);      % sigma-equiv or right-tail anchor pt
                    sig_equiv_neg = bins(1:ix1)  /sig1_neg;     % sigma-equiv of left-tail bins (adjusted for mean not exactly at 0).  Units are sigmas.
                    sig_equiv_pos = bins(ix2:end)/sig1_pos;     % sigma-equiv or right-tail bins
                    cdfs(1:ix1,   i, j) = max(cdfs(1:ix1,   i, j), normcdf(sig_equiv_neg)');    % replace tails with normal probability for sigma-equivalents, 
                    cdfs(ix2:end, i, j) = min(cdfs(ix2:end, i, j), normcdf(sig_equiv_pos)');    %       but only where normal curve is above calculated pdf
                    ix3 = find(cdfs(:,i,j) >= anchor_probs(2),1);           % find first bin inside left-tail 2nd anchor pt
                    ix4 = find(cdfs(:,i,j) <= anchor_probs(4),1,'last');    % find first bin inside right-tail 2nd anchor point
                    cdf_okflags(ix3:ix4,i,j)  = true;                       % flag all points in between as OK.
                end
            end
            
            pdfs = diff([zeros(1,yrlen,nsets); cdfs],1,1);
%             pdfs1=pdfs;
%             cdfs1=cdfs;
%             myflag=false;
            if (any(pdfs(:)<0))     % because we interpolated points in get_cdf_probs, we get an occasional negative  from rounding issues.
                minp = min(pdfs(:));
                if (minp < -1e-4)
                    obj.DP.warn_log('warning:  append_normal_cdf: max negative in CDF after get_cdf_probs: %9.6e\n', minp); 
%                     myflag = true;
                 end
                pdfs(pdfs < 0) = 1e-16;
            end
            cdfs = cumsum(pdfs);        % and because we've mucked around with the values, we have to make sure we have true CDFs when we're done.
            cdfs = cdfs ./ cdfs(end,:,:);
            pdfs = diff([zeros(1,yrlen,nsets); cdfs],1,1);

%             if (myflag)
%                 save("zzz.mat","pdfs","pdfs1","cdfs","cdfs1");
%                 fprintf("bailing!\n");
%             end
%  
            if (obj.DP.intern_plotflag(32))                % see notes in ARRM_V2_DataParams on creating internal plots.
                for j=1:nsets
                    figure(99);
                    clf;
                    subplot(1,3,1);
                    surf(1:365, obj.RP.out_bins, pdfs(:,:,j),'facecolor','g','edgecolor','none'); 
                    hold on;  
                    surf(1:365, obj.RP.out_bins, pdfs1(:,:,j),'facecolor','red','edgecolor','none'); 
                    hold off;            
                    set(gca,'zscale','log');
                    lights_on(true);

                    subplot(1,3,2);
                    surf(1:365, obj.RP.out_bins, cdfs(:,:,j),'facecolor','g','edgecolor','none'); 
                    hold on;  
                    surf(1:365, obj.RP.out_bins, cdfs1(:,:,j),'facecolor','red','edgecolor','none'); 
                    hold off;            
                    set(gca,'zscale','log');
                    lights_on(true);

                    subplot(1,3,3);
                    surf(1:365, obj.RP.out_bins, 1-cdfs(:,:,j),'facecolor','g','edgecolor','none'); 
                    hold on;  
                    surf(1:365, obj.RP.out_bins, 1-cdfs1(:,:,j),'facecolor','red','edgecolor','none'); 
                    hold off;            
                    set(gca,'zscale','log');
                    lights_on(true);
                    view([-135,35]);
                end
            end
        end
                
        function report_probcounts(probs, probcounts_base, probcounts_rolling, totcount_base, totcounts_rolling, lbl, lbl2, rolling_yrs)
            
            tcount_base = totcount_base;            % non-rounded
            tcounts_rolling = totcounts_rolling;    % non-rounded

            nsets = length(totcounts_rolling);
            totcount_base = round(totcount_base);
            totcounts_rolling = round(totcounts_rolling);
            totcounts_tot = round(sum(tcounts_rolling));
            fprintf('\n%s Probability Counts %s\n', lbl, lbl2);
            nprobs = length(probs);
            pcountsb = round(probcounts_base);
            pcountsr = round(probcounts_rolling);
            mid = ceil(nprobs/2);
            pcountsb(mid+1:end) = round(tcount_base-probcounts_base(mid+1:end));
            probFracsb = probcounts_base/tcount_base;
            probpctsb = probFracsb*100.0;
            if (nsets > 0)
                probFracsr=nan(nprobs,nsets);
                for i=1:nsets
                    pcountsr(mid+1:end,i) = round(tcounts_rolling(i)-probcounts_rolling(mid+1:end,i));
                    probFracsr(:,i) = probcounts_rolling(:,i)/tcounts_rolling(i);
                end
                probpctsr = probFracsr*100.0;
            end

            fprintf('          ');
            fprintf('%8.4f ', probs*100.0);
            fprintf('%8.4f ', 100);
            fprintf('\n');
            fprintf('          ');
            fprintf('%8d ', pcountsb);
            fprintf('%8d ', totcount_base);
            fprintf('\n');
            fprintf('          ');
            fprintf(' %7.3f%%', probpctsb);
            fprintf('  100.00');
            fprintf('\n\n')

            if (nsets > 0)
                fprintf('          ');
                fprintf('%8.4f ', probs*100.0);
                fprintf('%8.4f ', 100);
                fprintf('\n');
                for i=1:nsets
                    fprintf('   %4d   ', rolling_yrs(i));
                    fprintf('%8d ', pcountsr(:,i));
                    fprintf('%8d ', totcounts_rolling(i));
                    fprintf('\n');
                    fprintf('          ');
                    fprintf('%7.3f%% ', probpctsr(:,i));
                    fprintf('  100.00');
                    fprintf('\n')
                end

                pcntst     = sum(probcounts_rolling,2);
                pcountstot = round(pcntst);
                pcountstot(mid+1:end) = max(0, round(sum(tcounts_rolling)-pcntst(mid+1:end)));
                probFracstot = sum(probcounts_rolling,2)/sum(tcounts_rolling);
                probpctstot = probFracstot*100.0;

                fprintf('\n   total  ');
                fprintf('%8.4f ', probs*100.0);
                fprintf('%8.4f ', 100);
                fprintf('\n');
                fprintf('          ');
                fprintf('%8d ', pcountstot);
                fprintf('%8d ', totcounts_tot);
                fprintf('\n');
                fprintf('          ');
                fprintf('%7.3f%% ', probpctstot);
                fprintf('  100.00');
                fprintf('\n')
            end    

        end
        function yrs = rolling_step_years(obj, n)
            % returns year range ([start_yr, end_yr]) for rolling step n.            
            yrs = obj.DP.rolling_yrs;
             yrs(1) = max(obj.DP.data_yrs(1), yrs(1) + obj.DP.rolling_steps(n) - floor(obj.RP.pdf_yrstep/2));
            if (n ~=length(obj.DP.rolling_steps))
                yrs(2) = yrs(1) + obj.RP.pdf_yrstep-1;
            end
        end
%         
%         function [start_ix, end_ix, nyrs] = using_range(obj, yr_type, rolling_set)
%            % returns start and end indices of data to use for the specified year range (& rolling set, if given), and number of years.  
%            % yr_type can be string 'data_yrs','base_yrs','rolling_yrs', etc., or [start_yr, end_yr];
% 
%             if (exist('rolling_set','var') && ~isempty(rolling_set))
%                 yrs = obj.rolling_step_years(rolling_set);
%             elseif (isnumeric(yr_type))
%                 yrs = obj.DP.yr_type;
%             else
%                 yrs = obj.DP.(yr_type);
%             end
%             
%             nyrs = yrs(2)-yrs(1)+1;
%             start_ix = (yrs(1) - obj.DP.data_yrs(1)  )*obj.DP.yrlen+1;
%             end_ix   = (yrs(2) - obj.DP.data_yrs(1)+1)*obj.DP.yrlen;
%         end
%

        function [start_ix, end_ix, nyrs] = using_range(obj, yr_type, rolling_set)
%            % returns start and end indices of data to use for the specified year range (& rolling set, if given), and number of years.  
%            % yr_type can be string 'data_yrs','base_yrs','rolling_yrs', etc., or [start_yr, end_yr];
% 
            if (~exist('rolling_set','var')), rolling_set = []; end
            [start_ix, end_ix, nyrs] = obj.DP.using_range(yr_type, rolling_set, obj.RP.pdf_yrstep);
        end
        
        function det_y = detrend_valid_only(obj)

        % calculates trend from original data, and returns detrended signal.  Uses only years flagged as valid by trend_yr_flags
        %
        % avg value and trend are  calculated abd subtracted from signal, and also stored in the DA object for later use.  
        %
        % NOTE:  trend has base_avg subtracted, so it is zero-mean during base years only. That way when we add the
        %        the downscaled data back together, it matches the base observations.
        %        For Separate-Hist runs, the base_avg is actually the average of the model years.
        %       
        %
        %   Inputs:
        %       y               data to detrend.  Will detrend using only the trend_yrs data, but remove trend from all the data points. 
        %       RP              RunParams, containing:
        %         order           polynomial order to use (3 is good for 150-year data.  2nd order, parabolic, doesn't fit data well enough.
        %                                           (1 is good for 40- or 50-year historical period)
        %                                            0 for straight average
        %                                           '-1 :  data already detrended, or do-not-detrend.
        %         trend_yr_flags  boolean flag, 1 per year, flagging which years to include in trend calculation
        %                           (years with large fraction of NAs should be excluded to avoid biasing trend)
        %                           set to empty ( [] ) to use all years.
        %         yrlen           length of year, in days (usually 365).  must be integral value.
        %         nyrs            # of years of data over entire length of data
        %       DP              DataParams, containing:
        %         data_yrs        start, end years for entire data set
        %         trend_yrs       start, end years to use for trend calculations
        %         final_yrs       start, end years to retain in results.  (may pull more years for boundary issues...)
        %
        %   Outputs:
        %       det_y         original data, with trend subtracted, as a column vector.
        %       trend         daily (residual) trend over entire data period (not just trend_yrs) afte subtracting avg_base
        %                       Note:  trend is zero-mean during base years -- i.e., has avg_base removed.
        %       trend_params  polyfit params to recreate the long term trend.  trend = polyval(p,(0:(len-1))/yrlen,[],avg_trend);
        %       avg_trend     mean value of (residual) trend.  avg_trend+avg_base = average over entire trend period.

                % do polynomial fit to data of specified order,
                % then subtract the trend to separate out the long term trend and detrended signal

            data   = obj.raw_fix - obj.avg_base;
            data    = reshape(data, numel(data),1);
            [tr_start, tr_end] = obj.using_range('trend_yrs');
            trdata  = data(tr_start:tr_end);
            trlen   = length(trdata);

            obj.trend_params = obj.calc_trend_params(trdata);
            if (obj.insufficient_data)
                det_y = nan;
                return;
            end
            
            if (~isnan(obj.RP.trend_order) && obj.RP.trend_order >= 0)
                obj.trend = obj.recalc_trend();
                if (obj.insufficient_data)
                    det_y = nan;
                    return;
                end
                det_y = obj.raw_fix - obj.avg_base - obj.trend;
            else                        
                det_y = data - obj.avg_base;
                obj.trend = zeros(trlen(data),1);        
            end
        end 

        function trend_params = calc_trend_params(obj, trdata)
            % returns trend params, or nan if insufficient data to process (< 10 yrs with sufficient data, and not forced by manually setting trend_order)
            
                    % get the trend order and trend_yr_flags, if either is set to 'calc'.
            if (strncmpi(obj.RP.trend_order,'calc',4) || strncmpi(obj.RP.trend_yr_flags,'calc',4))
                [tr_ord, tr_flags] = ARRM_V2_DisaggregateSignal.calc_order(trdata, obj.RP.trend_thresh, obj.RP.yrlen);
                if (isnan(tr_ord))
                    trend_params = nan;
                    obj.insufficient_data = true;                    
                    return;
                end
                if (strncmpi(obj.RP.trend_order,'calc',4))
                    obj.RP.trend_order = tr_ord;
                end
                if (strncmpi(obj.RP.trend_yr_flags,'calc',4))
                    obj.RP.trend_yr_flags = tr_flags;
                end
            end    
            
            if (~isnan(obj.RP.trend_order) && obj.RP.trend_order >= 0)
                if (isempty(obj.RP.trend_yr_flags))
                    obj.RP.trend_yr_flags = true(obj.RP.nyrs,1);            
                elseif (isrow(obj.RP.trend_yr_flags))
                    obj.RP.trend_yr_flags = obj.RP.trend_yr_flags';
                end

                trlen   = length(trdata);
                x=(0:(trlen-1))'/obj.RP.yrlen;       % x is in (fractions of) years
                valid_flags = repmat(obj.RP.trend_yr_flags,obj.RP.yrlen,1);
                avg_trend=mean(trdata(valid_flags));

                p = polyfit(x(valid_flags),trdata(valid_flags),obj.RP.trend_order);

            else                        
                p = 0;
                avg_trend = mean(trdata);
            end
            len = length(obj.raw_fix);
            tr_start = obj.using_range('trend_yrs');
            trend_params = struct('p',p,'len',len,'yrlen',obj.RP.yrlen,'avg',avg_trend, 'start_yr', obj.DP.data_yrs(1), 'start_offset', tr_start, 'x','x=((1:len)-start_offset)/yrlen; % x is in (fractions of) years);  to reconstruct trend:  trend = polyval(p,x); % avg is average of trend during trend-years (it is zero-mean for base_yrs only)');
        end
        
        function trend = recalc_trend(obj)
        % calculates trend for all data, using trend params calculated on the valid trend data, using tr_start as the
        % zero-point.  So when we subtract this trend from the data, the base period will be zero-mean.
            if (isempty(obj.trend_params))
                trend = zeros(size(obj.raw_fix));
                return;
            end
            tr_start = obj.using_range('trend_yrs');
            len      = length(obj.raw_fix);
           
            xd = ((1:len)-tr_start)/obj.RP.yrlen;        % calc trend for entire data, not just the trend_yrs, but zero-point is at day tr_start
            trend = polyval(obj.trend_params.p,xd)';     % note:  trend has mean of zero for the base period.
        end

        function ok = check_valid_counts(obj, hist, valid_count)
            if (obj.DP.isPrecipRun)     % skip this for     
                    % For precip, this should check the counts of NAs in na_map and bail if too many NAs in any given
                    % month.  Skipping that for now...
                ok = true;
            else
                if (~exist('valid_count','var')), valid_count = obj.DP.monthly_valid_count; end
                if (~exist('hist','var')), hist = obj.base_hist; end

                nsets = size(hist,3);
                yrlen = size(hist,2);
                ndays = ceil(yrlen/12);
                nd    = yrlen/12;
                for i=1:nsets
                    d1 = floor((i-1)*nd) + 1;
                    d2 = min(yrlen, d1 + ndays);
                    npts = sum(sum(hist(:,d1:d2)));
                    if (~(npts >= valid_count)), ok = false; return; end
                end
                ok = true;
            end
        end
        
        function trend_params_daily = calc_daily_trend_params(obj)
           % calculate daily trend params.  This calculates the daily trend of the raw data.
           % You can create a surface from the combined set;  subtracting the basic climatology shows how the trend
           % varies over the entire year and over the time period.  Unlike the single trend-line, this uses all the 
           % data in the trend_yrs period, since we don't have to worry about biasing the trend due to missing data 
           % during one part of the year.
           
           yrlen = obj.RP.yrlen;
           if (isnan(obj.RP.trend_order) || isempty(obj.RP.trend_order) || obj.RP.trend_order == 0)
               p = zeros(yrlen,1);
               mu = mean(obj.raw_fix) - obj.avg_base;
           else
               % remove either the base climatology (if available) or the average climatology.
                [tr_start,tr_end] = obj.using_range('trend_yrs');
                data = reshape(obj.raw_fix(tr_start:tr_end), yrlen, []) - obj.avg_base;
                mu = mean(data(:));
                nyrs = size(data,2);

                nparms = obj.RP.trend_order+1;
                p = zeros(yrlen, nparms);
                x=0:(nyrs-1);
                for i=1:yrlen
                    p(i,:) = polyfit(x,data(i,:),obj.RP.trend_order);              
                end
                        % smooth any noisiness in the parameters, which will produce a smooth surface.
                        % faster than creating the surface, smoothing it, and calculating new parameters.

                for i=1:nparms      
                    p(:,i) = obj.DA_climatology(p(:,i), obj.RP.clim_nterms, obj.RP.clim_sig_terms, yrlen);
                end             
           end
           
           len = length(obj.raw_fix)/obj.RP.yrlen;
            yr_start = (obj.using_range('trend_yrs')-1)/obj.RP.yrlen;
            trend_params_daily = struct('p',p,'len',len,'yrlen',obj.RP.yrlen,'avg',mu, 'start_offset', yr_start, 'start_yr', obj.DP.data_yrs(1), 'x','x=((1:yrlen)-start_offset); % x is in (fractions of) years);  to reconstruct trend:  trend = polyval(p,x); % avg is average of trend during trend-years (it is zero-mean for base_yrs only)');
        end
        
        %function [ok, valid_counts, na_counts] = check_nans(DAs)
        function [ok, month_ok, npts, nas_in_counts, outlier_counts, far_outlier_counts, valid_counts, yr_range] = check_nans(obj, yr_type)
        %   counts the # of valid points and nas in the data for each of the DAs.
        %   Does NOT set insufficient_data flag.  That's up to whoever calls it.

            st_day_365 = [1, 32, 62, 92, 122, 153, 183, 214, 244, 274, 305, 335, 366];  % evenly divided, not actual months
            st_day_360 = [1, 31, 61, 91, 121, 151, 181, 211, 241, 271, 301, 331, 361];

%           yr_type = obj.DA_yrlbl;
            if (~strcmpi(yr_type,'rolling'))
                nsets = 1;
                [ix1, ix2, nyrs] = obj.using_range(yr_type);
                rolling=false;
            else
                nsets = length(obj.DP.rolling_steps);
                rolling=true;
            end
            month_ok            = false(12,nsets);
            npts                = zeros(12,nsets);
            nas_in_counts       = zeros(12,nsets);
            outlier_counts      = zeros(12,nsets);
            far_outlier_counts  = zeros(12,nsets);
            valid_counts        = zeros(12,nsets);
            yr_range            = zeros(2,nsets);
            
            if (obj.na_map_unknown)
                obj.na_map_in = isnan(obj.raw_data);
                obj.na_map_unknown = false;
            end

            for j=1:nsets   % check on model & downscaled is not robust.  Only checking for at least 600 points/month,
                            % not 600 points every x years.  (should really check for 600 every 25 years...)
        %        fprintf(obj.DP.fidlog, "%-5s ", lbls(i));
                if (rolling)
                    [ix1,ix2, nyrs] = obj.using_range(yr_type, j, obj.RP.pdf_yrstep);
                end        
                yrlen = obj.RP.yrlen;
                if (yrlen == 360)
                    st_day = st_day_360;
                else
                    st_day = st_day_365;
                end
                yr_range(:,j) = obj.DP.data_yrs(1) + [floor((ix1-1)/yrlen),ceil(ix2/yrlen)-1];
                na_map       = reshape(obj.na_map_in(ix1:ix2),       yrlen, nyrs);
                outliers     = reshape(obj.outlier_map(ix1:ix2),     yrlen, nyrs);
                far_outliers = reshape(obj.far_outlier_map(ix1:ix2), yrlen, nyrs);
                for i=1:12
                    s = st_day(i);
                    e = st_day(i+1)-1;
                    valid_counts(i,j)       = sum(sum(~na_map(s:e,:),2));
                    nas_in_counts(i,j)      = sum(sum(na_map(s:e,:),2));
                    outlier_counts(i,j)     = sum(sum(outliers(s:e,:),2));
                    far_outlier_counts(i,j) = sum(sum(far_outliers(s:e,:),2));
                    npts(i,j)               = (e-s+1)*nyrs;
                    if (nyrs < 25)
                        valid_thresh = valid_counts(i,j) / (.8*nyrs*(s-e+1));
                        month_ok(i,j) = nas_in_counts(i,j)/nyrs >= valid_thresh;
                    else
                        month_ok(i,j) = valid_counts(i,j) >= obj.DP.monthly_valid_count;
                    end
                end
            end        
            ok = ~(any(~month_ok(:)));
        end
        
                % getter functions for dependent data fields
        function typ = get.DAType(obj)
            typ = obj.DP.dspType;
        end
        
        function field = get.DA_yrlbl(obj)
            
            da_ix = find(strcmpi(obj.DAType,obj.DA_types),1);
            if (isempty(da_ix))
                field = 'data_yrs';
            else
                field = obj.yr_types{da_ix};
            end
        end
        
        function yrs = get.DA_yrs(obj)                     
            yrs= obj.DP.(obj.DA_yrlbl);
        end
        
        function  set_zeros_to_NAs(obj, fldnames, zero_thresh)
            %   for precip, where we want to exclude 0's from the analysis
            %   fldnames:   pass in array of fieldnames to have zero's set to NAs.  Default:  raw_data.
            %                   use "all" to set zeros in raw_data and all mapped arrays (see list below).
            %   zeros_map:  if missing, uses obj's zeros_map array.  Otherwise use zeros_map passed in.
            %                   if that's empty, it creates zeros_map array and stores map in obj.
            %
            %   returns the zeros_map used.
            
            if (~exist("zero_thresh","var") || isempty(zero_thresh))
                zero_thresh = 0;
            end
            if (isempty(obj.zeros_map))
                obj.zeros_map = (obj.raw_data <= zero_thresh);
            end
            zero_map = obj.zeros_map;
            
            if (strcmp(fldnames(1),"all"))
                fldnames = ["raw_data", "mapped_anoms", "mapped_anoms_im","outlier_mapped_anoms","mapped_output"];
            end
            for i=1:length(fldnames)
                if (~isempty(obj.(fldnames(i))))
                    obj.(fldnames(i))(zero_map) = nan;
                end
            end
        end
        
        function reset_zeros(obj, fldnames, zeros_map)
            %  This is for precip processing, where we set dry days to NAs throughout the processing.
            %  This resets original zeros back to zero from NAs, using the zeros_map.
            %  fldnames:    optional list of fldnames to reset.  default:  raw_data and mapped_output
            %                   use "all" to reset zeros in raw_data and all mapped arrays (see list below).
            %   zeros_map:  optional map to use.  If missing or empty, uses obj's zeros_map array.
            
            if (~exist("zeros_map", "var") || isempty(zeros_map)), zeros_map = obj.zeros_map; end
            
            if (~exist("fldnames", "var") || isempty(fldnames))
                fldnames = ["raw_data","mapped_output"];
            elseif (strcmp(fldnames(1),"all"))
                fldnames = ["raw_data", "mapped_anoms", "mapped_anoms_im","outlier_mapped_anoms","mapped_output"];
            end
            for i=1:length(fldnames)
                if (~isempty(obj.(fldnames(i))))
                    obj.(fldnames(i))(zeros_map) = 0;
                end
            end
        end            
        
        
    end
        
    methods(Static)
        
        function [data, na_map] = fix_nans(data, yrlen, nyrs)
        %       Replaces nans with daily average, then low=pass-filters data w/ fft filtering,
        %       and replaces nan's with LPF'd result.  
        %       Should really detrend the data first, though, to avoid problems from circular convolution.
        %       Should also probably be a little smarter and replace single-NAs with linear interpolated values.
        %       Revisit this, Ian!

            npts = yrlen*nyrs;
            data = reshape(data,npts,1);
            na_map = isnan(data);
            if (sum(na_map) == 0)                
                return; 
            end
            data   = reshape(data,   yrlen,nyrs);
            na_map = reshape(na_map, yrlen, nyrs);

                % first, replace NAs with daily average
            daily_avg = mean(data,2,"omitnan");
            nnandays = sum(isnan(daily_avg));
                % if we have any days w/ no data points, linear-interpolate to get a reasonable estimate.
            if (nnandays > 0)
                davg = fillmissing([daily_avg; daily_avg; daily_avg],'linear');
                daily_avg = davg(yrlen+1:2*yrlen);                
            end
            daily_avg = repmat(daily_avg,1,nyrs);
            data(na_map) = daily_avg(na_map);

            data = reshape(data, npts,1);
            na_map = reshape(na_map, npts,1);
                % now we low-pass filter the data and replace the original NAs with the low-pass value for the specific day.
            dfilt = fft_filt([data;flipud(data)], 6*nyrs);         % use first 6 terms of mirrored FFT.  mirror to avoid step function @ ends of series.
            dfilt =reshape(dfilt(1:npts), nyrs,yrlen);
            data(na_map) = dfilt(na_map);

        end
        
                    
        function [ord, valid_flags, fracs, seq_range, seq_used] = calc_order(data, valid_yr_na_thresh, yrlen)
        % function [ord, yrs, valid_flags, fracs, seq_range, seq_used] = ARRM_V2_calc_order(data_yrs, data, valid_frac)
        %
        %   Function to calculate the polynomial fit order for ARRM_V2 OBS data.
        %   Note:  ARRM_V2_run calls this function to get the valid_flags to flag which years have enough data, even if
        %           the user specified the trend_orders.  ARRM_V2_run ignores the calculated order if specified in the run
        %           Params.
        %
        %   Inputs:
        %       data_yrs    year (x-values) for each data point  (e.g., just the year column of datevecs)
        %                       can also be [yr1,yr_end] or yr1:yr_end
        %       data        y-values, usually temperature
        %       valid_yr_na_thresh  fraction of data points needed in a given year to consider the year as valid.
        %                       E.g. 0.5 (use year if at least 50% of points are valie) or 0.75 (use if at least 75% of points valid.
        %
        %   Outputs
        %       ord         order:  0, 1, 3 or nan
        %                       3   >= 50 yrs (fit3=50) valid data
        %                       1   >= 30 yrs valid data or >= 20 sequential yrs data
        %                       0   >= 10 yrs valid data.
        %                       nan <->  "insufficient data" (< 10 yrs with enough data)
        %       yrs         list of (all) years in the data  (1 for each year, not 1 for each point)
        %       valid_yrs   boolean flags;  true if year has >= valid_yr_na_thresh good points.  
        %       fracs       fraction of points valid in each year
        %       seq_range   range of years for longest sequence in the data
        %       seq_used    boolean flag.  True if > 20 sequential years, but < 30 valid years total.
        %       
            % data must be:  complete years, starting jan 1st.
            % if data not complete years, prepend/append NAs
        %
        %   To calculate the trend from this data:
        %                 d = reshape(data, yrlen, nyrs);
        %                 yr_means = nanmean(d)';
        %                 [ord, yr_flags, valid_yrs, fracs, seq_range, maxseq, use_seq_range] = ARRM_V2_DisaggregateSignal.calc_order(yrs, temps, order_thresh);
        %                 [p,S,mu] = polyfit(valid_yrs,yr_means(yr_flags)',ord);
        %                 data_yrs = yr_range(1):yr_range(2);
        %                 trend_all   = polyval(p, data_yrs,S,mu);
        % 

%            VERIFY THIS, IAN!  

            npts = numel(data);
            nyrs = npts / yrlen;
            if (mod(nyrs,1.0) ~= 0)
                fprintf(2,"error: callstack: \n");
                dbstack();
                error('error:  ARRM_V2_calc_order():  data length (%d) is not multiple of yrlen days (%d)', npts, yrlen);
            end
        %     valid_yr_na_thresh = 0.5;     % year's data valid if > 50% valid points
        %                           % valid_yr_na_thresh is now set in ARRM_V2_run_params.m, and passed in above.  
            fit3 = 50;              % do 3rd order if >= 50 yrs valid data
            fit1 = 30;              % do 1st order if >= 30 yrs valid data
            fit1_seq = 20;          % do 1st order if >= 20 yrs sequential valid data
            minvalid = 10;           % mininum # of years to consider data usable.

            data = reshape(data, yrlen, nyrs);
            fracs = (sum(~isnan(data)) / yrlen)';

            valid_flags = fracs > valid_yr_na_thresh;
            nvalid = sum(valid_flags);
            seq_used = false;
            [seq_len, seq_range] = ARRM_V2_DisaggregateSignal.longest_sequential(1:nyrs, valid_flags);
            if (nvalid >= fit3)
                ord = 3;
            elseif (nvalid > fit1)
                ord = 1;
            elseif (seq_len > fit1_seq)
                seq_used = true;
                ord = 1;
            elseif (nvalid < minvalid)
                ord = nan;
            else
                ord = 0;
            end
        end

        function [seq_len, seq_range] = longest_sequential(yrs, valid_flags)

            ix1 = find(valid_flags, 1);
            ix2 = find(valid_flags, 1, 'last');

            seq_len = 0;
            nseq = 0;
            istart=0;
            seq_range=[];

            if (isempty(ix1)), return; end

            for i=ix1:ix2
                if (~valid_flags(i))
                    nseq = 0;
                    istart=0;
                else
                    nseq = nseq + 1;
                    if (istart==0)
                        istart=i;
                    end
                    if (nseq >seq_len)
                        seq_len = nseq;
                        seq_range = [yrs(istart),yrs(i)];
                    end
                end
            end
        end      
        

        function [ysurf_out] = lpf_surf(ysurf, sig)
        % [ysurf_out] = ARRM_V2_DisaggregateSignal.lpf_surf(ysurf, sig)
        %
        % function to low-pass filter along the year dimension (horizontally) with gaussian filter.
        %   This routine uses fourier domain to circular-convolve a gaussian(sig) with each row of ysurf.
        %
        %   inputs:
        %       ysurf       s/b 2-D surface of size yrlen  X nyrs 
        %       sig         is time-domain sigma, in years, to use along the year dimension.
        %   outputs
        %       ysurf_out   filtered surface

            [ndays,nyrs] = size(ysurf);
            SIG = calc_SIGMA(sig, nyrs);
            mid=ceil((nyrs+1)/2);       % for length 256, s/b 127.
            G = gauss(nyrs, SIG,mid);   % create gaussian
            G=G/max(G);                % normalize to max of 1.
            G = ifftshift(G);           % shift so centered at (1) instead of midpoint, so fft doesn't introduct a phase shift
            GG = repmat(G,ndays,1);     % replicate it to same size as data surface.
            ysurf_out = ifft(GG .* fft(ysurf,nyrs,2), nyrs, 2); 
            ysurf_out = real(ysurf_out);              % SHOULD be purely real, but limited precision math means we have some tiny (~10**-15) imag. part, which we discard.

        end

        function y_ext = extend_data(y, yrlen, n_ext_yrs, nadd_start, nadd_end)
        % y_ext = ARRM_V2_DisaggregateSignal.extend_data(y, n_ext_yrs, nadd_start, nadd_end)
        %
        % extends the data by prepending the average of the first years 
        % and appending the average of the last years, so we have 'valid' data to convolve with.
        %   Inputs:
        %       y               data to extend
        %       n_ext_yrs     # of years to average at beginning & end to extend with
        %       nadd_start      # of years to add at start
        %       nadd_end        # of years to add at end
        %   Ouptputs:
        %       y_ext           extended data.
        %

            nyrs = length(y(:))/yrlen;
            min_ext_yrs = min(nyrs, n_ext_yrs);
            clen=yrlen*min_ext_yrs;
            ystart=reshape(y(1:clen),yrlen,min_ext_yrs);                          % reshape 1st years to 365xnstart
            start_mean = mean(ystart,2);                                        % then calculate mean  (result is 365 daily means)
            start_mean = repmat(start_mean,nadd_start,1);                           % and duplicate it
            yend=reshape(y((end-clen+1):end),yrlen,min_ext_yrs);                  % and do the same for the end few years
            end_mean = mean(yend,2);
            end_mean = repmat(end_mean,nadd_end,1);
            y_ext = [start_mean; y; end_mean];                                  % append mean of early years at front of signal, and mean of last years at end.
        end

        function pdfs_t  = normalize_pdfs_means(xq, clim_mus, pdfs_t, fwd)
        % shifts the pdfs to a common mean, or shifts back.
        % Called by create_kde_pdfs_lpf
            [~,     yrlen] = size(pdfs_t);
            mu_bar = mean(clim_mus);
            if (fwd)
                for i=1:yrlen
                    xx = xq - clim_mus(i) + mu_bar;           % subtracting the mean here will center the adjusted pdf.
                    pdfs_t(:,i) = interp1((xx), pdfs_t(:,i), xq,     'pchip',0);
                    pdfs_t(:,i) = pdfs_t(:,i)/sum(pdfs_t(:,i));     % renormalize.
                end
            else
                for i=1:yrlen
                    xx = xq - clim_mus(i) + mu_bar;           % locate it back to its original offset from zero-mean.
                    pdfs_t(:,i) = interp1(xq,    pdfs_t(:,i), xx,     'pchip',0);
                    pdfs_t(:,i) = pdfs_t(:,i)/sum(pdfs_t(:,i));     % renormalize.
                end
            end
        end
        
        function pdfs_t  = normalize_pdfs_medians(xq, pdf_medians, pdfs_t, fwd)
        % shifts the pdfs to a common mean, or shifts back.
        % Called by create_kde_pdfs_lpf
            [~,     yrlen] = size(pdfs_t);
            med_bar = mean(pdf_medians);
            if (fwd)
                for i=1:yrlen
                    try
                    xx = xq - pdf_medians(i) + med_bar;           % subtracting the mean here will center the adjusted pdf.
                    pdfs_t(:,i) = interp1(xx, pdfs_t(:,i), xq,     'pchip',0);
                    pdfs_t(:,i) = pdfs_t(:,i)/sum(pdfs_t(:,i));     % renormalize.
                    catch
                        oops();
                    end
                end
            else
                for i=1:yrlen
                    xx = xq - pdf_medians(i) + med_bar;           % locate it back to its original offset from zero-mean.
                    pdfs_t(:,i) = interp1(xq, pdfs_t(:,i), xx,     'pchip',0);
                    pdfs_t(:,i) = pdfs_t(:,i)/sum(pdfs_t(:,i));     % renormalize.
                end
            end
        end
        
%        function [anoms, sigmas] = normalize_anoms_by_sigmas(anoms, yrlen, fwd, sigmas, ttl, figno)
function [anoms, sigmas] = normalize_anoms_by_sigmas(anoms, yrlen, fwd, sigmas)
            
            % normalizes the anomalies 
            
            anoms = reshape(anoms, yrlen, []);
%           a1 = anoms;
            if (fwd)
                sigmas = std(anoms, [], 2);
%              s1 = sigmas;
                sigmas = rect_filt(sigmas, ceil(yrlen/18), 1, true);
%              s1a = sigmas;
%                sigmas = climatology(sigmas,4,1, yrlen);       % low pass filter, keeping lowest 4 terms only,  
%               sigmas = climatology(sigmas,12,8, yrlen);       % low pass filter, keeping lowest 4 terms only,  
                anoms = anoms ./ sigmas;
            else
                anoms = anoms .* sigmas;
            end

%             a2 = anoms;
%             s2a = std(anoms, [],2);
%             s2 = rect_filt(s2a, ceil(yrlen/18), 1, true);
% %           s2  = climatology(s2a,4,1, yrlen);
              anoms = reshape(anoms, [],1);
% 
%             if (fwd)
%                 figure(figno); 
%                 subplot(3,1,1);
%                 plot(1:yrlen,s1, 1:yrlen, sigmas, 1:yrlen, s2a, 1:yrlen, s2);
%                 title(ttl,"interpreter","none");
%                 subplot(3,1,2);
%                 nyrs=size(a1,2);
%                 y=reshape(repmat(1:yrlen,1,nyrs),1,length(anoms));
%                 x=reshape(repmat(1:nyrs,yrlen,1),1,21900);
%                 scatter3(x,y,a1(:),5,"filled")
%                 title("orig anoms");
%                 view([90,0]);
%                 subplot(3,1,3)
%                 scatter3(x,y,a2(:),5,"filled")
%                 title("orig anoms");
%                 view([90,0]);
%             end
        end
        
%         function pdfs_t  = normalize_pdfs_by_sigmas(x, clim_mus, clim_sigs, pdfs_t, sig0, fwd)
%         % normalize or denormalize pdfs based on sig
%         % Called by create_kde_pdfs_lpf
%         %       (and also create_kde_pdfs_norm_sigmas, which is superceded bycreate_kde_pdfs_lpf)
%  % z=0;          
% %             mu0 = mean(clim_mus);
% %             [nbins, yrlen] = size(pdfs_t);
%             [~,     yrlen] = size(pdfs_t);
% %             xx = repmat(x,1,yrlen) -  repmat(clim_mus,nbins,1);
% %             xq = repmat(x,1,yrlen) .* repmat(clim_sigs,nbins,1) / sig0;
% %             s1=zeros(yrlen,1);
% %             s2=zeros(yrlen,1);
%             if (fwd)
%                 for i=1:yrlen
%                     xx = x - clim_mus(i);           % subtracting the mean here will center the adjusted pdf.
%                     xq = (x) * clim_sigs(i)/sig0;
% %                     s1(i) = sqrt(sum(xx.^2 .* pdfs_t(:,i)) - mus(i).^2);
%                     pdfs_t(:,i) = interp1((xx), pdfs_t(:,i), xq,     'pchip',0);
%                     pdfs_t(:,i) = pdfs_t(:,i)/sum(pdfs_t(:,i));     % renormalize.
% %                     s2(i) = sqrt(sum(xx.^2 .* pdfs_t(:,i)) - mean(pdfs_t(:,i)).^2);
% %                     s1(i) = sqrt(sum(xx(:,i).^2 .* pdfs_t(:,i)) - mus(i).^2);
% %                     pdfs_t(:,i) = interp1(xx(:,i), pdfs_t(:,i), xq(:,i),'pchip',0);                     
% %                     s2(i) = sqrt(sum(xx(:,i).^2 .* pdfs_t(:,i)) - mean(pdfs_t(:,i)).^2);
% %                   fprintf(obj.DP.fidlog, '%3d: s1,s2; %8.5f %8.5f  / %8.5f  %8.5f\n', i, s1(i), s2(i), s1(i)/s2(i), clim_sigs(i)/sig0);
% 
%                 end
% %                 cs=climatology(s2,4,2,128);
% %                 figure(99);  
% %                 subplot(2,1,1);  
% %                 plot(1:yrlen,s1,'b-', 1:yrlen,s2, 'r-', [1,128],[mean(s2),mean(s2)],'r-', 1:yrlen, cs,'r-', ...
% %                      1:yrlen, clim_sigs,'g-',[1,128],[mean(clim_sigs),mean(clim_sigs)],'g-')   ;             
%             else
% %                 pdfs1=pdfs_t;
%                 for i=1:yrlen
%                     xx = x - clim_mus(i);           % locate it back to its original offset from zero-mean.
%                     xq = (x) * clim_sigs(i)/sig0;
%                     pdfs_t(:,i) = interp1(xq,    pdfs_t(:,i), xx,     'pchip',0);
%                     pdfs_t(:,i) = pdfs_t(:,i)/sum(pdfs_t(:,i));     % renormalize.
% %                     if (i<10)
% %                         figure(99);
% %                         plot(pdfs1(:,i));
% %                         hold on;
% %                         plot(pdfs_t(:,i));
% %                         hold off;
% %                     end
% %                     pdfs_t(:,i) = interp1(xq(:,i), pdfs_t(:,i), xx(:,i), 'pchip',0);
% %                     pdfs_t(:,i) = pdfs_t(:,i)/sum(pdfs_t(:,i));     % renormalize.
%                 end
%             end
%         end

%         function [tot_mem, raw_mem, histogram_mem, pdf_mem, problines_mem] = calc_mem_usage(RP, DP)
% % THIS NEEDS UPDATING, IAN!
%             nrolling_steps  = length(DP.rolling_steps);
%             nbins           = length(RP.bins);
%             yrlen           = RP.yrlen;
%             nrolling_yrs    = DP.rolling_yrs(2) - DP.rolling_yrs(1) + 1;
%             do_rolling      = nrolling_yrs > 0;
%             nall_yrs        = DP.data_yrs(2) - DP.data_yrs(1) + 1;
%         %   nbase_yrs       = DSP.base_yrs(2) - DSP.base_yrs(1) +1;
%             do_pdfs         = RP.do_pdfs; 
%             do_problines    = RP.do_pdfs & RP.do_problines;
%             nprobs          = length(RP.probs);
%             do_keep_all     = DP.keep_all;
% 
%             raw_mem             = RP.yrlen * nall_yrs * 4;  % Raw data only.  
% 
%             anoms_mem           = yrlen * nall_yrs * 8;
%             na_map_in_mem       = yrlen * nall_yrs * 1;
%             trend_mem           = yrlen * nall_yrs * 8;
% 
%             base_clim_mem       = yrlen * 8;
%             base_hist_mem       = nbins * yrlen * 8;
%             moving_clim_mem     = yrlen * nrolling_yrs * 8;
%             moving_clim_1d_mem  = yrlen * nrolling_yrs * 8;
%             rolling_clim_mem    = yrlen * 8 * nrolling_steps;
%             rolling_hist_mem    = nbins * yrlen * nrolling_steps * 8;
% 
%             pdf_base_mem            = nbins * yrlen * 8;
%             CDF_base_mem            = nbins * yrlen * 8;
%             base_okflags_mem        = nbins * yrlen * 1;
%             pdf_rolling_mem         = nbins * yrlen * nrolling_steps * 8;
%             CDF_rolling_mem         = nbins * yrlen * nrolling_steps * 8;
%             rolling_okflags_mem     = nbins * yrlen * nrolling_steps * 1;
% 
%             problines_base_mem      = yrlen * nprobs * 8;
%             pdf_xlines_base_mem     = yrlen * nprobs * 8;
%             CDF_xlines_base_mem     = yrlen * nprobs * 8;
%             probcounts_base_mem     = nprobs * 8;
% 
%             problines_rolling_mem   = yrlen * nprobs * nrolling_steps * 8;
%             pdf_xlines_rolling_mem  = yrlen * nprobs * nrolling_steps * 8;
%             CDF_xlines_rolling_mem  = yrlen * nprobs * nrolling_steps * 8;
%             probcounts_rolling_mem  = nprobs * nrolling_steps * 8;
%             totcounts_rolling_mem   = (nrolling_steps+1) * 8;
% 
%             anoms_all_mem           = yrlen * nall_yrs * 8;
%             trend_all_mem           = yrlen * nall_yrs * 8;
%             moving_clim_all_mem     = yrlen * nall_yrs * 8;
%             moving_clim_1d_all_mem  = yrlen * nall_yrs * 8;
% 
% 
%             histogram_mem = anoms_mem + na_map_mem + trend_mem + ...
%                             base_clim_mem + base_hist_mem + moving_clim_mem + moving_clim_1d_mem + ...
%                             (rolling_clim_mem + rolling_hist_mem) * do_rolling;
% 
%             pdf_mem = (pdf_base_mem + CDF_base_mem + base_okflags_mem + ...
%                       (pdf_rolling_mem + CDF_rolling_mem + rolling_okflags_mem) * do_rolling) * do_pdfs;
% 
%             problines_mem = (problines_base_mem + pdf_xlines_base_mem + CDF_xlines_base_mem + probcounts_base_mem + 1 + ...
%                              (problines_rolling_mem + pdf_xlines_rolling_mem + CDF_xlines_rolling_mem + probcounts_rolling_mem + totcounts_rolling_mem) * do_rolling) * do_problines;
% 
%             keep_mem = (anoms_all_mem + trend_all_mem + moving_clim_all_mem + moving_clim_1d_all_mem) * do_keep_all;
% 
%             tot_mem = raw_mem + histogram_mem + pdf_mem + problines_mem + keep_mem;
% 
%         end   
        
%         function  pdfs_t = smooth_along_time(yrlen, sigrange, pdfs_t, use_freq, do_variable_gaussians, mu0, sig1, sigA, sigB, sig_scale)
%             % smooths along time.
%             % If do_variable_gaussians is true, then it increases the standard deviation of the smoothing kernel by a
%             % factor of sig_scale for isolines between sigA and sigB, where sig1 is 1 sigma, in bins, and mu0 is center (in bins)
%             % Notes:    sigma   percentile inside
%             %        +/-1.645   .90
%             %        +/-1.96    .95
%             %        +/-2.575   .99
%             %           3.29    .999
%             %
%             %   sigrange of 9 days  * sqrt(12) ->  31 days rectangular filter
%             %       tripling @ .999 -> 93 days
%             %   sigrange of 14 days * sqrt(12) ->  48 days rectangular filter
%             %       doubling @ .999 -> 96 days
%             %       tripling @ .999 -> 144 days
%             
%             [nbins, pdf_yrlen, nsets] = size(pdfs_t);

%             sigrange = sigrange * pdf_yrlen/365;
%             for j=1:nsets
%                 if (use_freq)
%                     error('DisaggregateSignal:smooth_along_time: use_freq: not ready for this yet.');
%                 else
%                     if (~do_variable_gaussians)
%                         G = GAUSS_fft(pdf_yrlen, sigrange, true);
%                         pdfs_t(:,:,j) = ifft(G .* fft(pdfs_t(:,:,j), yrlen, 2),yrlen, 2);      % smooth along time (day-of-year) axis
%                     else
% 
%                         ctr = mu0;                          % center, in bin counts
%                         sigs = ((0:(nbins-1))-ctr)/sig1;        % distance from mean of each bin, relative to sig1.
% 
%                         ksigs=zeros(nbins,1);
%                         for r=1:nbins
%                             sig = abs(sigs(r));
%                             frac = max(0, min(1, (sig-sigA)/(sigB-sigA)));
%                             krnl_sig = sigrange * (1  + (sig_scale-1)*frac);
%                             ksigs(r)=krnl_sig;
%                             G = GAUSS_fft(pdf_yrlen, krnl_sig,true);     % FFT of gaussian, length pdf_yrlen & sigma.
%                             pdfs_t(r,:,j) = ifft(G .* fft(pdfs_t(r,:,j), pdf_yrlen, 2), pdf_yrlen, 2);  % check this, ian!
%                         end
% %                         figure(99);
% %                         plot(sigs, ksigs);
% %                         title(sprintf('%.4f, from %.4f to %.4f, x %.2f', sigrange, sigA, sigB, sig_scale));
% %                         grid on;
% %                         fprintf('here!\n');
%                     end
%                 end
%             end
%         end       
%         
%         function  pdfs_t = smooth_along_time(pdfs_t, yrlen, sigrange, sigrange_probs, do_linear)
%             % smooths along time.
%             % If do_variable_gaussians is true, then it increases the standard deviation of the smoothing kernel by a
%             % factor of sig_scale for isolines between sigA and sigB, where sig1 is 1 sigma, in bins, and mu0 is center (in bins)
%             % Notes:    sigma   percentile inside
%             %        +/-1.645   .90
%             %        +/-1.96    .95
%             %        +/-2.575   .99
%             %           3.29    .999
%             %
%             %   sigrange of 9 days  * sqrt(12) ->  31 days rectangular filter
%             %       tripling @ .999 -> 93 days
%             %   sigrange of 14 days * sqrt(12) ->  48 days rectangular filter
%             %       doubling @ .999 -> 96 days
%             %       tripling @ .999 -> 144 days
%             %   sigrange of 17 days * sqrt(12) -> 
%             %
%             %   probs:  [.001,.1,.9, .999]  stretch points for expanding sigrange.  
%             %               inside .1 & .9, = sigrange
%             %               outside .001 & .0999, sig_scale * sig_range
%             %               in between:  scaled linearly.
%             
%             %
%             % CURRENTLY NOT USED.  This is called by create_kde_pdfs_norm_sigmas, which does smoothing along the time
%             % dimension after normalizing by Std Dev.  Superceded by smoothing in probability space, done by create_kde_pdfs_lpf
%             
%             [nbins, pdf_yrlen, nsets] = size(pdfs_t);
%             sigrange = sigrange * pdf_yrlen/yrlen;
%             for j=1:nsets
%                 if (length(sigrange)==1)
%                     G = GAUSS_fft(pdf_yrlen, sigrange, true);
%                     pdfs_t(:,:,j) = ifft(G .* fft(pdfs_t(:,:,j), yrlen, 2),yrlen, 2);      % smooth along time (day-of-year) axis
%                 else
% 
%                     my_pdf = sum(pdfs_t,2);   
%                     my_pdf = my_pdf / sum(my_pdf);
%                     my_cdf = cumsum(my_pdf,'omitnan');
%                     ix1 = find(my_cdf < .01* sigrange_probs(1),1,'last');
%                     ix2 = find(1-my_cdf < .01 * (1-sigrange_probs(4)),1,'first'); 
%                     prob_locs = interp1(my_cdf(ix1:ix2), ix1:ix2, sigrange_probs, "pchip");
% 
%                     ksigs=zeros(nbins,1);
%                     sigrange_dif = sigrange(2) - sigrange(1);
%                     for r=1:nbins
%                         if (r <= prob_locs(1) || r >= prob_locs(4))
%                             frac = 1;
%                         elseif (r >= prob_locs(2) && r <= prob_locs(3))
%                             frac = 0;
%                         elseif (r <prob_locs(2))
%                             if (do_linear)
%                                 frac = (r-prob_locs(2))/(prob_locs(1)-prob_locs(2));
%                             else
%                                 frac = (my_cdf(r)-sigrange_probs(2))/(sigrange_probs(1)-sigrange_probs(2));
%                             end
%                         else
%                             if (do_linear)
%                                 frac = (r-prob_locs(3))/(prob_locs(4)-prob_locs(3));
%                             else
%                                 frac = (my_cdf(r)-sigrange_probs(3))/(sigrange_probs(4)-sigrange_probs(3));
%                             end
%                         end
%                         krnl_sig = sigrange  + sigrange_dif*frac;
%                         ksigs(r)=krnl_sig;
%                         G = GAUSS_fft(pdf_yrlen, krnl_sig,true);     % FFT of gaussian, length pdf_yrlen & sigma.
%                         pdfs_t(r,:,j) = ifft(G .* fft(pdfs_t(r,:,j), pdf_yrlen, 2), pdf_yrlen, 2);  % check this, ian!
%                     end
% %                         figure(99);
% %                         plot(-49.975:.25:50, ksigs);
% %                         grid on;
%                 end
%             end
%         end 
        
        function [sigmin, fourier_sigmin] = silversig_min(N)
            % minimum sigma for silverman's rule of thumb to keep the gaussian from getting too large in the fourier
            % domain.  If the silverman ROT sigma is less than this, then the resulting gaussian doesn't go 
            % in close enough to zero in the frequency domain to filter the data properly.
            % N is number of bins in the histogram
            % sigmin is the minimum sigma (in bin counts) for a gaussian of length N in the time domain
            %
            %       example:  N=401, sigmin = 2.6744 in time domain.
            %                         sigma = N/(2*pi*sigma) = 23.8637 in fourier domain
            %                         g=gauss(401, 23.8637) produces a gaussian centered in a vector length 401 
            %                                               with sigma 23.8637, g(1) = 9.34e-18
            %   output is about 2.63 for 1000, 2.67 for 400, 2.71 for 200.
            %   For temperature data, bin sizes less than 400 will generate silversigs too small.  Limiting
            %   silversig to this value will result in oversmoothing in the KDE step.
            %
            % This equation is approximate, and valid in range 100 to 2000.
            % determined empirically, to produce a value of < ~1e-17 at the ends of the Fourier domain gaussian for a given N
            % (for 2001, g(1) is 1.3e-17, slightly above 1.0e-17, but still satisfactorily low.)
            %
            sigmin = 2.97728-.050531*log(N);

            fourier_sigmin = N/2/pi/sigmin;

        end
        
        function plot_sigmas(fignum, mus, clim_mus,  sigs, clim_sigs)

            yrlen = size(mus,2);
            figure(fignum);
            clf;
            subplot(3,1,1);
            plot(1:yrlen,mus(1,:),'b:',1:yrlen, mus(2,:),'r:','linewidth',2);
            hold on;
            plot(1:yrlen,clim_mus(1,:),'b-',1:yrlen,clim_mus(2,:),'r-');
            hold off;
            title('distribution means');
            xlabel('time');
            ylabel('means');
            grid on;
            subplot(3,1,2);
            h=plot(1:yrlen,sigs(1,:),'b:',1:yrlen, sigs(2,:),'r:', 1:yrlen, sigs(3,:),'g:', 'linewidth',2);
            hold on;
            plot(1:yrlen,clim_sigs(1,:),'b-', 1:yrlen, clim_sigs(2,:),'r-', 1:yrlen, clim_sigs(3,:),'g-');
            hold off;                    
            title('distribution sigmas');
            xlabel('time');
            ylabel('standard deviation');
            grid on;
            legend(h,'raw','normalized','kde-smoothed');

            % should add figures of skewness & kurtosis.

            subplot(3,1,3)
            plot(1:yrlen,clim_sigs(1,:),'b-',1:yrlen,clim_sigs(2,:), 'r-', 1:yrlen, clim_sigs(3,:), 'g-', 1:yrlen, sigs(4,:),'m-', 1:yrlen, sigs(5,:),'k:', 1:yrlen, sigs(6,:),'k--');
            grid on;

            title('distribution sigmas');
            xlabel('time');
            ylabel('standard deviation');
        end
        
    end
    
    
    
    
end

