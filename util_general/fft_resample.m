function [yout] = fft_resample(y, outlen, terms, dim)
%function [yout] = [yout] = fft_resample(y, outlen, terms, dim)
%   Resamples y to length outlen using fft.
%   operates along 1st non-singleton dimension of y.
%   returns:  
%       yout        resampled data
%
%   inputs:
%       y           data to resample or filter
%       outlen      length of desired output
%       terms       either an array of 1's identifying frequency terms to keep, or
%                   a scalar representing the highest frequency term to keep.  
%                       note:  term 0 is DC, 1 is 1 cycle, 2 is 2 cycles, etc.
%                       note:  because of ambiguity, terms == 1 means 'DC + first term'.
%                       note:  if terms is non-scalar, it can either be a vector (applied to all fft sets)
%                              or matrix with flags unique to each set.
%                       note:  terms only needs to be as long as needed to specify all terms desired.  Will be filled
%                                   with zeros as needed.
%

    if (length(size(y))>2)
        throw(MException('FFT_RESAMPLE','input dimension too high.  Can only operate on vectors or 2-D matrices'));
    end
      
    if (~exist('terms','var'))
        terms=[];
    end
    if (~exist('dim','var') || isempty(dim))
        dim = find(size(y)>1,1);
    end
    n=size(y,dim);

    scaling=outlen/n;
    if (isempty(terms) && outlen==n)      % nothing to do.  copy input & return.
        yout=y;
        return;
    end
    
    if (mod(n,outlen)==0 && isempty(terms))     % n is multiple of outlen, and no filtering, so subsample and return
        yout = y(1:(n/outlen):end);
        return;
    end
     
    Y = fft(y, n, dim);

    yout = ifft_resample(Y, outlen, terms, dim) * scaling;

end
