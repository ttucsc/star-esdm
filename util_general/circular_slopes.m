function slopes = circular_slopes(d, dx)
%   returns slope of data.  
%   if dx is missing, scales like a sinusoid:
%       x = 0-2:pi, y ranging from -1 to 1.
%
%   NOTE:  d should be smoothed before calling this function. 
%          One approach:  d = climatology(your_signal, nterms, sig_term).  6,2 reasonable for nterms, sig_term.
%           See comments in climatology(...).  nterms is # of frequency terms to use for low-pass filter.  sig_term is
%           gaussian rolloff of ideal filter to reduce ringing.

    npts = length(d);

    if (iscolumn(d))
        d = [d(end);d;d(1)];
    else
        d = [d(end),d,d(1)];
    end
    if (nargin ==1)   % scale to sinusoid
        dx = 2*pi/npts;                                 % diff is 2 steps wide.
        d = (d-min(d))/range(d);                        % so dx s/b 4 PI, and y range of 2.
                                                        % for simplicity, we use 2 pi and range of 1
    end
    df = d(1:(end-2)) - d(3:end);
    slopes = df./dx;
end