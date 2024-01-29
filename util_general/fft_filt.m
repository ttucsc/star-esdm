function [yout] = fft_filt(y, terms, dim)
%[yout] = fft_filt(y, terms, dim)
%   
%   Filters y in fourier domain, keeping only terms specified
%   operates along 1st non-singleton dimension of y.
%   returns:  
%       yout        filtered data
%
%   inputs:
%       y           data to filter

%       terms       either an array of 1's identifying frequency terms to keep, or
%                   an array of frequency weights or
%                   a scalar representing the highest frequency term to keep.  
%                       -1:  keep all terms
%                       note:  term=0 keeps DC, term=1 keeps DC + 1st term cycle, term=2 keeps DC & first 2 cycles, etc.
%                       note:  because of ambiguity, terms == 1 means 'DC + first term'.
%                       note:  if terms is non-scalar, it can either be a vector (applied to all fft sets)
%                              or matrix with flags or weights unique to each set.
%                       note:  terms only needs to be as long as needed to specify all terms desired.  Will be filled
%                                   with zeros as needed.
%

    if (length(size(y))>2)
        throw(MException('ICSF:FFT_FILT','input dimension too high.  Can only operate on vectors or 2-D matrices'));
    end
      
    if (~exist('terms','var'))
        help(mfilename('fullpath'));
        yout=[];
        return;
    end
    if (~exist('dim','var') || isempty_s(dim))
        dim = find(size(y)>1,1);
    end
    n=size(y,dim);
     
    Y = fft(y, n, dim);

    yout = ifft_filt(Y, terms, dim);

end
