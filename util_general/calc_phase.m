function phase = calc_phase(x, nterms, precision, yrlen, y)
%   Returns phase of x relative to y.  
%       if y not provided, calculates vs pure cosine of length (length(x))
%   nterms:  # of terms of FFT to use.
%   precision:  precision to calculate phase, in terms of 1 step ( #steps = length(x))
%       e.g., to get phase to 1/10 of a step, set precision to 10.
%
%   Phase returned is in units of x, not in radians or degrees.
%       i.e., it gives the relative offset from either y or a cosine of length yrlen.
%       
    len = length(x);
    if (~exist('y','var') || isempty(y))
        y = mean(x) + cos(2*pi*(0:len-1)/len) * range(x)/2;
    end
    phase = ic_corr_1d(y, x, nterms, precision);
    if (exist('yrlen','var') || isempty(yrlen))
        phase = mod(phase, yrlen);
    end            
end

