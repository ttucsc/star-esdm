function G = GAUSS_fft( len, sigma, do_shift )
% G = GAUSS_fft( len, sigma, do_shift )
%
%   Returns G, the FFT of a gaussian with length len & sigma.
%   If do_shift is true, shifts G so it is centered at the origin.
%   This results in the FFT of a gaussian with no phase delay.

%       MATH Notes:
%
%       FFT of gaussian is a gaussian, so we can actually do this more
%       quickly by calculating the corresponding sigma in the frequency
%       domain and creating a gaussian with that sigma directly.
%       This saves having to take the FFT of a series.
%   1.  sigma-0 of gaussian s.t. sigma(fft(gaussian)), in freq = sigma in time:
%           (for n points): sqrt(n/(2*pi))
%   2.  sigma (freq. domain) for a sigma=1 (time domain) is .1592 = 1/(2*pi)
%   3.  sigma-fft for arbitrary time domain sigma:  sigma-fft = sigma0^2 / sigma(time)
%           alt:  sigma_fft = N/(2*pi*sigma_time)
%   3.a.    and same for given sigma-fft:  sigma(time domain) = sigma0^2 / sigma(freq domain)
%           alt:  sigma_time = N/(2*pi*sigma_fft)
%   
%__________________________________________________________________________________________________

    if (mod(len,2)==1)
        shft = floor(len/2);
    else
        shft = len/2-1;
    end

%           The obvious way:  create a gaussian, take it's FFT.
%     if (do_shift)
%         G = fft(circshift(gauss(len,sigma),-shft, 2));
%         G = G/G(1);     % normalize to max of 1
%     else
%         G = fft(gauss(len,sigma));
%         G = G/G(shft+1);
%     end
% 

%           The faster way:  create a gaussian w/ freq-domain sigma.
%           This is because the fft of a gaussian is a gaussian.
%           See math notes above.
%     sigma_0 =  sqrt(len/2/pi);
%     SIGMA = (sigma_0^2)/sigma = N/(2*pi*sigma);                  % sigma in freq domain
    SIGMA = len/(2*pi*sigma);
    if (do_shift)
        G = circshift(gauss(len,SIGMA,shft+1),-shft,2);    
        G = G/G(1);     % normalize to max of 1;
    else
        G = gauss(len,SIGMA,shft+1);
        G = G/G(shft+1);     % normalize to max of 1;
    end


end

