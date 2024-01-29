 function [clim, clim_ratio, phase, avg, clim_raw] = climatology(y, nterms, sig_terms, yrlen)        
% clim = climatology(y, nterms, sig_terms,  ylen) 
%
% calculates the average daily value (365 x 1) for the data in y  (365*nyrs x 1)  (or yrlen x 1 and yrlen*nyrs x 1)
% then low-pass filters it with a filter defined by nterms and sig_terms.
%
%   Note that the average here is a straight average for each day of the year (not gaussian-weighted), which is
%   then smoothed by circular convolution along the 365 days.
%
%   Data should be an exact multiple of yrlen.  Use 365.25 for data with leap-days included, and make sure the data is a
%   multiple of 4 years long.
%
% Returns a 365x1 vector representing the smoothed daily average of the data in y.
%   Inputs:
%       y           data to work with.  must be of multiple of 365 days long
%       nterms      equivalent # of frequency terms to retain
%       sig_term    sigma for gaussian to smoothe rectangular filter.
%                       good values to use for nterms & sig_term are 5 & 2.0
%                       this gives the equivalent filtering of an ideal
%                       filter keeping 1st 6 terms of fft (up to 6 cycles
%                       per year)
%                       for pure gaussian filtering, use nterms=0, sig_term set to desired gaussian frequency domain sigma to smooth with.  
%
%                   See "math notes" for equations to go from
%                   either time-domain gaussian sigma or idea-filter-equivalent sigma.
%
%   Outputs:
%       clim        low-pass filtered climatology
%       clim_ratio  inverse of the ratio of 1st frequency term to the rms sum of the remaining terms
%                       NOTE:  'remaining terms' are after low-pass filtering, so high-frequency terms will be reduced or eliminated.
%                       SMALLER values of clim_ratio imply stronger annual cycle.
%                       Useful for determining whether there is a strong annual cycle effect.
%       phase       phase of signal relative to pure cosine wave.
%                       NOTE: units of phase are days, not radians or degrees, and is always positive.
%                             phase == +10 means the climatology is delayed 10 days relative to a cosine wave.
%                             high-in-January-low-in-July will have positive phase near 0;
%                             high-in-December-low-in-June will have positive phase near 365;
%                             low-in-January-high-in-July will have postivie phase near 182.
%                             phase == 0 means no shift.
%
%                       Temperature climatology averaged over CONUS region:  
%                           nearly sinusoidal, with clim_ratio of about .09 and phase of 208 (~July 27) 
%
%                       precision of phase is 1/10 day.  (Accuracy may be much less if clim_ratio is not small...
%                           i.e. if climatology does not have a strong annual component. 
%
%                       To compare 2 climatologies and get their phase, use ic_corr_1d(clim1, clim2...) directly, as it
%                       will correlate the two phases against each other, whereas the phase here is compared to a pure
%                       cosine wave, so in only comparing the primary frequency.
%                           NOTE:  ic_corr_1d(...) returns phase in range -length/2 < phase <= length/2.  i.e., not
%                           always positive.
%
%       avg             mean value of y(:)
%       clim_raw        unfiltered climatology.  Simply the average value for each day-of-year.


    avg = mean(y(:), "omitnan");
    [nr,nc] = size(y);
    if (nr == 1 || nc == 1)        
        if (~exist('yrlen','var') || isempty(yrlen)), yrlen = 365; end
            % if keeping all terms, just do straight average:
        nyrs = length(y)/yrlen;    
        clim = mean(reshape(y,yrlen,nyrs),2,"omitnan");     % get average daily value for data range.
    else
        clim = mean(y,2, "omitnan");
    end
    clim_raw = clim;
    
        % Create the filter, and low-pass filter the daily averages
    FILT = calc_filter(nterms, sig_terms, 1, yrlen);
    [clim, CLIM] = lpf_FILT(clim, FILT);

    if (nargout > 1)
        clim_ratio = calc_clim_ratio(CLIM);
    else
        clim_ratio = nan;
    end
    
        
    if (nargout > 2)
        if (any(isnan(clim)))
            phase = nan;
        else
%             avg = mean(clim);
%             mycos = cos(2*pi*(0:yrlen-1)/yrlen);
%             phase = ic_corr_1d(mycos, (clim-avg)/range(clim)*2, nterms, 10);
%             phase = mod(phase, yrlen);
            phase = calc_phase(clim, nterms, 10, yrlen);

        end
    end
            
    
end

function r = calc_clim_ratio(Y)
%r = calc_clim_ratio(Y)
%
%       Returns the climatology ratio:  the ratio of the magnitude of the primary term of the FFT to the rms sum of the
%       remaining terms.  Note that if Y is already filtered, this will be really be the ratio of the first term to the
%       next few terms
%
%       NOTE IAN:  this is 1/2 the value from my original program!
%
    if (abs(Y(2))==0); r = inf; return; end
    sinmag = abs(Y(2));             % magnitude of the first term
    mid = ceil(length(Y)/2);        % end of the positive terms
    r = sqrt(sum(abs(Y(3:mid)).^2)/sinmag^2);   % ratio of 1st term to remainder (ignoring DC term)
end


