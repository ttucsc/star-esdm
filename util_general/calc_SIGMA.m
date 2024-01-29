function [SIG, sig0] = calc_SIGMA(sig, seq_len)
% [SIG, sig0] = calc_SIGMA(sig, seq_len)
%
% Same as calc_sigmas, but starts from time-domain sequence length, rather than # of terms.
% function to calculate the frequency-domain sigma for a given time-domain sigma, for series of length seq_len
%   inputs:
%       sig         time-domain sigma
%       seq_len     length of sequence
%   outputs:
%       SIG         frequency-domain sigma
%       sig0        same as sig0's in calc_sigmas...the scaling term relating fourier-domain and time-domain sigmas.

%       MATH Notes for gaussian:
% notes re sigmas, gaussians & fft's of gaussians:
%   1.  sigma of uniform distribution of length n: 
%           for continuous:  n/2/sqrt(3) = n/sqrt(12)
%           for discrete:  sqrt((n^2-1)/12) , 
%           which, as m -> large, -> n/sqrt(12).
%   1.a. n, length of uniform dist for given sigma:  
%           n = sqrt(12*(sigma^2 + 1)) 
%   2.  sigma-0 of gaussian s.t. sigma(fft(gaussian)), in freq = sigma in time (for n points): sqrt(n/(2*pi))
%   3.  sigma (freq. domain) for a sigma=1 (time domain) is .1592 = 1/(2*pi)
%   4.  sigma-fft for arbitrary time domain sigma:  sigma-fft = sigma0^2 / sigma(time)
%   4.a.    and same for given sigma-fft:  sigma(time domain) = sigma0^2 / sigma(freq domain)
%   

    
    sig0 = sqrt(seq_len/2/pi);
    SIG  = (sig0^2)/sig;

end
