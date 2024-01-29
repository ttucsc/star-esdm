function [ord, yrs, valid_flags, fracs, seq_range, seq_used] = ARRM_V2_calc_order(data_yrs, data, valid_yr_na_thresh, yrlen)
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
%                 [ord, yr_flags, valid_yrs, fracs, seq_range, maxseq, use_seq_range] = calc_order(yrs, temps, order_thresh);
%                 [p,S,mu] = polyfit(valid_yrs,yr_means(yr_flags)',ord);
%                 data_yrs = yr_range(1):yr_range(2);
%                 trend_all   = polyval(p, data_yrs,S,mu);
% 

    [nr,nc] = size(data);
    if (nr ~= 1 && nc ~= 1)
        yrlen = nr;
        nyrs = nc;
    else
        
        if (~exist('yrlen','var') || isempty(yrlen)), yrlen = 365; end
        npts = length(data(:));
        nyrs = npts / yrlen;
        if (mod(nyrs,1.0) ~= 0)
            throw(MException('ICSF:BAD_DATA','error:  ARRM_V2_calc_order():  data length is not multiple of yrlen days'));
        end
    end    
%     valid_yr_na_thresh = 0.5;     % year's data valid if > 50% valid points
%                           % valid_yr_na_thresh is now set in ARRM_V2_run_params.m, and passed in above.  
    fit3 = 50;              % do 3rd order if >= 50 yrs valid data
    fit1 = 30;              % do 1st order if >= 30 yrs valid data
    fit1_seq = 20;          % do 1st order if >= 20 yrs sequential valid data
    minvalid = 10;           % mininum # of years to consider data usable.
    
    
    if (length(data_yrs)==2)
        yrs = data_yrs(1):data_yrs(2);  % given 1st and last year
    elseif (length(data_yrs) == nyrs)
        yrs = data_yrs;                 % given list of years, one per year
    else
        yrs = data_yrs(1:yrlen:end);      % list of years, one for each data point.
    end
        
    if (length(yrs) ~= nyrs)
        throw(MException('ICSF:BAD_DATA','error:  calc_order():  invalid year data'));
    end
    
    data = reshape(data, yrlen, nyrs);
    fracs = (sum(~isnan(data)) / yrlen)';

    valid_flags = fracs > valid_yr_na_thresh;
    nvalid = sum(valid_flags);
    seq_used = false;
    [seq_len, seq_range] = longest_sequential(yrs, valid_flags);
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
