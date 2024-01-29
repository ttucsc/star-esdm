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
    % determined empirically, to produce a value of < ~1e-17 and the ends of the Fourier domain gaussian for a given N
    % (for 2001, g(1) is 1.3e-17, slightly above 1.0e-17, but still satisfactorily low.)
    %
    sigmin = 2.97728-.050531*log(N);
    
    fourier_sigmin = N/2/pi/sigmin;
        
end
