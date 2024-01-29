function clim = climatology_3d_gauss_V2(data, years, nclim_years, nterms, sig_terms, n_ext_years)
%        clim = climatology_3d_gauss_V2(data, years, nclim_years, nterms, sig_terms, n_ext_years)
%
% function to calculate rolling climatology.  Long term trend should be removed first.
%   Creates 3D climatology surface from rolling low-pass gaussian filter of specified data.
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
%           years           start- and end-years of entire dataset ( [1950,2100] is good)
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
%           n_ext_years    data will be extended to length of 256 years using average of 1st n_ext_years and last n_ext_years.
%
%       Outputs
%
%           clim            long-term climatology, as linear sequence
%                               smoothed as a 2-D surface, first along the day, then smoothed along the year, then
%                               reshaped back to a linear sequence.  To reshape into surface:  reshape(clim,365, nyrs);
%                               Does NOT have the base climatology removed.
%
% Data is extended to a n even power-of-two years, then after all smoothing, results clipped back to original years.
% FFT of even power of 2 is faster than fft of 151 or 191
%
%   7/17 icsf detrending removed, done as separate step now.
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
            
    nyears = years(2)-years(1) + 1;
    dyears = length(data)/365;
    if (nyears ~= dyears)
        throw(MException('ICSF:BAD_INPUT',sprintf('data length (%d, %f yrs) and # of years (%d) do not match.', length(data), dyears, nyears)));
    end
        
            % Extend the data either end to a power-of-2 length years so we have something to convolve with on the ends.
            % This appends the first few years' mean climatology at the beginning, and last few years' mean climatology
            % at the end.
            % Because of the efficiency of the power-of-2 FFT, it is faster to extend the series to a power of 2 years than
            % just adding 21 years either end to get 193 years total.
    
    len=64;                
    nadd = len - nyears;        
    while (nadd < 2*n_ext_years)    % find smallest power of 2 long enough for extended data.
        len = len*2;
        nadd = len - nyears;
    end
    nadd_start = floor(nadd/2);     % # of years to extend at start
    nadd_end   = nadd - nadd_start; % # of years to extend at end.
          
    y_ext = extend_data(data, n_ext_years, nadd_start, nadd_end);     % data, extended so circular convolution of FFT doesn't affect opposite ends of data.

    mpts = length(y_ext);
    myrs = mpts/365;            % this should be 256 (or 512)
    
            % now, low-pass filter entire series with gaussian along the days,
            % then reshape to a 2-D surface, and low-pass filter along the years dimension.
            
                    % get Fourier domain filter for extended sequence for smoothing along the days. 
    FILT_long = calc_filter(nterms, sig_terms, myrs);
%                     % get frequency sigma for extended sequence for smoothing along the days.
%     [~,    SIG_LONG, ~    ] = calc_sigma(myrs*nterms, myrs*365);        % ignore some of the returned values...
            % convolve original (extended) signal along the days w/ gaussian based on nterms
            % we do the convolution here in the fourier domain.  This is OK because we've extended the
            % signal at either end so there's no problem with the circular effect of fft.
            % (as long as n_ext_years*365 >> sigt_long.)

%     y2 = lpf_SIGMA(y_ext, SIG_LONG);
    y2 = lpf_FILT(y_ext, FILT_long);
 
            % now, reshape and lpf along the years.
            % See Math Notes above.
    ysurf = reshape(y2, 365, myrs);
    year_sig = sqrt((nclim_years^2-1)/12);
    clim_ext = lpf_surf(ysurf, year_sig);

%---------------------------------------
%           For now:  skip removing shifts.  We can add this in later if we find it is helpful.
%           But removing shift is computationally expensive, and probably gains very little.
    
%         % calculate and remove any phase shift in the signal.
%     shifts = calc_shifts(clim_base, clim_ext, year_sig, day_sig);   %smooths by yrsig & daysig, then calculate shifts. 
%          
%     clim_ext = remove_shift(clim_ext_raw, shifts);
%             
            % remove the data extension from start and end of shifts.
%     shifts = shifts(n_ext_years+(1:nyears));
%---------------------------------------

            % remove the data extension from start and end
    clim  = clim_ext(:,nadd_start+(1:nyears));
    
end

function [ysurf_out] = lpf_surf(ysurf, sig)
% [ysurf_out] = lpf_surf(ysurf, sig)
%
% function to low-pass filter along the year dimension (horizontally) with gaussian filter.
%   This routine uses fourier domain to circular-convolve a gaussian(sig) with each row of ysurf.
%
%   inputs:
%       ysurf       s/b 2-D surface w/ 365 rows and nyears columns.
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

function y_ext = extend_data(y, n_ext_years, nadd_start, nadd_end)
% y_ext = extend_data(y, n_ext_years, nadd_start, nadd_end)
%
% extends the data by prepending the average of the first years 
% and appending the average of the last years, so we have 'valid' data to convolve with.
%   Inputs:
%       y               data to extend
%       n_ext_years     # of years to average at beginning & end to extend with
%       nadd_start      # of years to add at start
%       nadd_end        # of years to add at end
%   Ouptputs:
%       y_ext           extended data.
%

    nyrs = length(y(:))/365;
    min_ext_years = min(nyrs, n_ext_years);
    clen=365*min_ext_years;
    ystart=reshape(y(1:clen),365,min_ext_years);                          % reshape 1st years to 365xnstart
    start_mean = mean(ystart,2);                                        % then calculate mean  (result is 365 daily means)
    start_mean = repmat(start_mean,nadd_start,1);                           % and duplicate it
    yend=reshape(y((end-clen+1):end),365,min_ext_years);                  % and do the same for the end few years
    end_mean = mean(yend,2);
    end_mean = repmat(end_mean,nadd_end,1);
    y_ext = [start_mean; y; end_mean];                                  % append mean of early years at front of signal, and mean of last years at end.
end

% 
%---------------------------------------Code for calculating and removing shifts--------------------------
%         % needs to be smarter.  This is inefficient.
% function shifts = calc_shifts(clim_base_raw, clim_ext_raw, sigday, sigyear)
% %
% %   smooths climatology surface in both directions, then does correlation to find shift
% 
%     mlen = round(7 * sigyear);  
%     myrs = size(clim_ext_raw,2);
%     if (mod(mlen,2)==0); mlen = mlen+1; end;
%     gday = ifftshift(gauss(365, sigday))';
%     gyear = gauss(mlen, sigyear);
%     
%     clim_base = cconv(clim_base_raw, gday, 365);
%     clim_ext = zeros(size(clim_ext_raw));
%     for iyr = 1:myrs
%         clim_ext(:,iyr) = cconv(clim_ext_raw(:,iyr), gday, 365);        % circular convolution along the day axis
%     end
%     for iday=1:365
%         clim_ext(iday,:) = conv(clim_ext(iday,:), gyear, 'same');       % but regular convolution along the year axis
%     end
%     
%     shifts = zeros(1,myrs);
%     for iyr=1:myrs
%         shifts(iyr) = ic_corr_1d(clim_base, clim_ext(:,iyr), [], 100);
%     end
% end
%  
% function clim_ext = remove_shift(clim_ext, shifts)
% %
% %   returns clim. surface, w/ each year shifted (back) by -shifts(i)
% %
%     myrs = size(clim_ext,2);
%     for i=1:myrs
%         clim_ext(:,i) = ph_shift(clim_ext(:,i), -shifts(i));
%     end
% end

% function [ shifted ] = ph_shift( y, shift_dist)
% % [ shifted ] = ph_shift( y, shift_dist, fignum )
% %
% %   circularly shifts signal y by shift_dist elements (rounded to nearest .01 elements)
% %
% %   y           input signal
% %   shift_dist  # of elements to shift (can be fractional)
% 
%     Y = fft(y);
%     len = length(y);
%     mid = ceil(length(Y)/2);
%     if (isrow(Y))
%         Y2=[Y(1:mid),zeros(1,99*len),Y(mid+1:end)];     % extend by factor of 100
%     else
%         Y2=[Y(1:mid);zeros(99*len,1);Y(mid+1:end)];
%     end
%     y2 = 100*ifft(Y2);
%     sh = round(shift_dist*100);
%     if (isrow(y))
%         y3=circshift(y2,sh,2);
%     else
%         y3=circshift(y2,sh);
%     end
%     shifted=y3(1:100:end);
% end    

% function [ shift ] = ic_corr_1d( base, changed, nterms, nterp, no_DC, fignum )
% %   Does fft-based correlation to get subpixel registration on 1-D data
% %
% %   Correlation done by standard fourier-domain correlation, equivalent of circular correlation is spatial domain.
% %   This version works well if signal is cyclical and changed signal is reasonably similar to base signal.
% %   Does not produce a spike at correlation point...instead provides a smooth continuum, and simply returns the point
% %   with the highest correlation value.
% %
% %   does:  ifft( fft(changed) .* conj(fft(base)) ), then finds max.
% %
% %               using only the first nterms terms of the fft.
% %
% %   If nterp > 1, expands signal by factor of nterp, so effectively interpolates to sub-cell precision.
% %
% %   Assumes changed and base are same size, and of similar signal levels.
% %
% %   To interpolate to sub-step ('sub-pixel') accuracy, set nterp to expand signal as needed.  (nterp = 10 will return result to
% %   nearest 1/10 cell.
% %   
% %       (As an alternative, see e.g. phase_corr_1d(...), based on Foroosh, "Extension of Phase Correlation to Subpixel
% %       Registration", IEEE 2002.)
% %
% % intputs:
% %   base        original signal
% %   changed     shifted signal.  response is amount that this has shifted relative to base
% %
% %   nterms      # of terms to keep (1 = DC, 2 = DC + 1 cycle, etc.);  0 or empty:  keep all terms
% %   nterp       # of times to expand signal to find max.  
% %   no_DC       exclude DC term (sets average of each signal to zero)
% %   fignum      figure to plot in.  if absent, no plot is generated.
% %
% % 
%     if (~exist('nterms','var') || isempty(nterms) || nterms < 1); nterms= 0;     end;
%     if (~exist('nterp', 'var') || isempty(nterp)  || nterp  < 1); nterp = 1;     end;
%     if (~exist('no_DC', 'var'));                                  no_DC = false; end;
%     
%     if (~isrow(base))
%         base = base';
%     end
%     if (~isrow(changed))
%         changed = changed';
%     end
%     
%     len = length(base);
%     
%     B = fft(base);
%     C = fft(changed);
%     if (no_DC)
%         B(1)=0;
%         C(1)=0;
%     end
%     
%     BC = C .* conj(B);
%     
%     if ( nterms > 0)
%         bc = ifft_resample(BC, nterp*len, nterp*(nterms-1));
%     else
%         bc = ifft_resample(BC, nterp*len);
%     end
%     
%     ix = find(bc == max(bc),1);
%     
%     if (exist('fignum','var'))
%         bf = fft_filt(base,nterms-1);       % fft_filt's nterms doesn't count DC term.
%         cf = fft_filt(changed, nterms-1);
%         figure(fignum);
%         subplot(2,1,1);
%         plot(1:length(base), base, 1:length(bf), bf, 1:length(changed), changed, 1:length(cf), cf);
%         subplot(2,1,2);
%         plot(1:length(bc), bc, [ix,ix],[min(bc),max(bc)]);
%     end
% 
%     shift = (ix-1)/nterp;
%     
%     if (shift > len/2)
%         shift = shift - len;
%     end
% end
