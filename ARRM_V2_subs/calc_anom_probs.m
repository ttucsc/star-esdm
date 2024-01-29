function [probs, mapped_nas, far_outliers, outliers] = calc_anom_probs(anoms, bins, CDF, na_map, far_outlier_thresh, outlier_thresh, isPrecipRun, prcp_min, prcp_median)
%       Calculates probability values for each anomaly in each CDF step.
%
%       returns probabilities and a map of the anomalies that fall outside min_thresh.
%       Anomalies outside the mappable range are set to nan, and flagged in mapped_nas.
%       Anomalies outside min_thresh are flagged in far_outliers (including mapped_nas).  
%       Anomalies outside outlier_thresh are flagged in outliers (including far_outliers and mapped_nas).
%
%       added 7/1/23 icsf:  precip anoms below min thresh are not mapped to NAs.  Instead, they are set to
%       probability of prcp_min, the min precip detectable.
%
%       Inputs:
%           anoms           anomalies (raw-data - mean - trend)
%           bins            bin values CDF (for interpolating)
%           CDF             Cumulative Distribution Function (size yrlen X nbins)
%           na_map          boolean map, true where raw input was originally nan or missing
%           far_outlier_thresh      probability beyond which to define point as far-outlier.  
%                               This can be used to invoke the alternate probability calculation (see below).  Set this to
%                               the point beyond which your CDF is unreliable,  usually somewhere between 1e-3 and 1e-6
%                               optional.  supply if your code is reading expecting far_outliers output array
%           outlier_thresh  probability beyond which to label a point an outlier
%                               (optional.  supply if your code uses outliers output array
%           
%
%       Outputs:
%           probs           probabilities for each anomaly, or nan if anomaly was originally nan or it cannot be
%                               calculated
%           mapped_nas      boolean map, true where probability mapping failed (interpolation produced NAs)
%           far_outliers    boolean map. true where probability is beyond min_thresh
%           outliers        boolean map.  true where probability is beyond outlier_thresh
%                               NOTE:  outliers includes far_outliers.
%                               NOTE:  outliers defined at top (1-thresh) as well as bottom of probabilities.
%
%_____________________________old version's output:
%           na_info:  
%               (:,1)   index of point outside pthresh_vals  (or outlier_thresh).
%               (:,2)   anomaly value for point
%               (:,3)   ratio of anomaly / pthresh_vals(day,1-or-2).  negative for below mean, positive for above mean.
%                           this is for gathering statistics on the outliers...shows how much further out than the day's
%                           pthresh_val the anomaly is.  Can also be used to calculate an alternate mapping point when
%                           downscaling outliers.
%               (:,4)   original interpolated probability (replaced in probs(...) array with adjusted probability
%                           approx. if actual prob. is nan, & based on (:,3) & pthresh
%               (:,5)   mapped probability (where the basic interpolation would put the point).
%
%           outlier_info:
%               same info as for na_info.  This array simply identifies "outliers", and does not substitute an alternate
%               value for the identiried outliers.
%
%       pthresh_vals is an array of values, usually for the 2.5% and 97.5% lines (~2-sigma probabilities).
%           which is used (see above) to calculate a normalized measure of how far out an outlier is.
%           It can be used to calculate a vaguely reasonable value when mapping points out in the wings.
%           pthresh_vals can  be a single column (such as std. deviations), or 2 columns, first for values < mean, 
%           2nd for values greater than the mean.  This allows for correcting differently if the data is skewed.
%       pthresh is probability point beyond which a point is considered an outlier.  (defaults of 1e-3 and 1-(1e-3) )
%

    [~, yrlen] = size(CDF);
    nyrs = length(anoms)/yrlen;
    anoms = reshape(anoms, yrlen, nyrs);
    na_map = reshape(na_map, yrlen, nyrs);
    probs    = zeros(size(anoms));
    
    if (length(outlier_thresh) == 1)
        outlier_thresh = sort([outlier_thresh, 1-outlier_thresh]);
    end
    if (length(far_outlier_thresh) == 1)
        far_outlier_thresh = sort([far_outlier_thresh, 1-far_outlier_thresh]);
    end

    okflags = ~isnan(CDF);  % valid region of CDF, not flags of anomalies.

    for i=1:yrlen           % for each day of year

                         % find probability for each model anomaly.  
        try
            flags = okflags(:,i);
            probs(i,:) = interp1(bins(flags), CDF(flags,i), anoms(i,:), 'makima',nan); % maybe need linear, Ian?
            if (isPrecipRun)
                    % fix any small precip amounts that either map to NAs
                    % or map to less than prcp_min.
                prcp_min_prob = interp1(bins(flags), CDF(flags,i), prcp_min, 'makima',nan);
                if (isnan(prcp_min_prob))
                    error("error: prcp_min_prob is NAN");
                end
                
                drizzle = (probs(i,:) <= prcp_min_prob) | (isnan(probs(i,:)) & anoms(i,:)<prcp_median(i)) & anoms(i,:) > 0; % set prob to prcp_min's prob if not dry day and mapped to NAN or very small probability
                if (any(drizzle))
                    probs(i,drizzle) = prcp_min_prob;
                end
            end

        catch
            oops();
        end
    end
    
    probs    = reshape(probs,   yrlen*nyrs,1);
    
    if (nargout > 1)
        mapped_nas = isnan(probs)& ~na_map(:);      % points outside min_thresh which weren't originally nans.
        if (nargout > 2)
            far_outliers = probs < far_outlier_thresh(1) | probs > far_outlier_thresh(2);
            if (nargout > 3)
                outliers = probs <     outlier_thresh(1) | probs >     outlier_thresh(2);
            end
        end
    end    
end


