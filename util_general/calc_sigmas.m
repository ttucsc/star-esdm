function [sig_t, SIG_F, sig_0] = calc_sigmas(nterms, len)
% [sig_t, SIG_F, sig_0] = calc_sigmas(nterms, len)
%
% function to calculate sigmas for given # of fourier terms, for a sequence of length len.
%   inputs:
%       nterms          # of positive fourier terms to keep for sequence (not counting DC term)
%       len             length of series to calculate the sigmas for.
%   outputs
%       sig_t           sigma for time domain gaussian, # of days, 
%       SIG_F           sigma for freq domain gaussian, # of cycles
%       sig_0           sigma0 for gaussian (value where sig_t == SIG_F)
%
%   Four our purposes, we're interested in specifying the our smoothing as a function of the number of terms we want to
%   keep for a 365-day year.  A reasonable value to use is 6.  This gives us up to 6 cycles per year, which can be
%   thought of as defining our climatology as events of ~2 months duration or longer;  anything shorter than that we'll
%   consider a short-term 'weather' variation.  A gaussian filter with the equivalent frequency-domain sigma will avoid
%   the ringing in the time domain that an ideal rectangular filter will produce.  Other research shows that 6 terms
%   captures 95% of the power of a typical "seasonal" event.

%       MATH Notes:
%       The following gives an explanation of how to get a sigma for a gaussian with the equivalent standard deviation
%       of a rectangular filter of a given length.  There are 2 sigmas of interest:  the sigma for a time-domain
%       gaussian to convolve with the time-domain signal, or the sigma in the frequency domain, to multiply
%       (point-by-point) with the fourier transform of the time-domain signal.  There is a useful value, sigma-0, which
%       is the sigma at which the time domain and frequency domain sigmas are the same.
% notes re sigmas, gaussians & fft's of gaussians:
%   1.  std dev (sigma) of uniform distribution of length n: 
%           for continuous:  n/2/sqrt(3) = n/sqrt(12)
%           for discrete:  sqrt(((n^2-1))/12) , 
%           which, as m -> large, -> n/sqrt(12).
%   1.a. n, corresponding length of uniform (rectangular) distribution for given sigma  (solving discrete for n):  
%           n = sqrt(12*(sigma^2) + 1) 
%   2.  sigma-0 of gaussian s.t. sigma(fft(gaussian)), in freq = sigma in time (for n points): sqrt(n/(2*pi))
%   3.  sigma (freq. domain) for a sigma=1 (time domain) is .1592 = 1/(2*pi)
%   4.  sigma-fft for arbitrary time domain sigma:  sigma-fft = sigma0^2 / sigma(time)
%   4.a.    and same for given sigma-fft:  sigma(time domain) = sigma0^2 / sigma(freq domain)
%
%   revised icsf 3/10/2016 to fix error in math.  had n^12/12-1 instead of (n^2-1)/12
%   

        
        sig_0 = sqrt(len/2/pi);
        SIG_F = sqrt((((2*nterms+1)^2)-1)/12);    % nterms is # of positive frequencies.   need DC term and negative freq's, so 2*nterms+1
        sig_t = (sig_0^2)/SIG_F;

end

