function [yout] = ifft_filt(Y, terms, dim)
% [yout] = ifft_filt(Y, outlen, terms, dim)

% Filters and inverse-FFT's Y in fourier domain, keeping only terms specified, starting from FFT of some data
%   See fft_filt(...) for complete task of filtering from time domain.
%   operates along 1st non-singleton dimension of y.
%   returns:  
%       yout        filtered data
%
%   inputs:
%       Y           fourier domain (FFT) of data to filter

%       terms       either an array of 1's identifying frequency terms to keep, or
%                   -1 :  keep all terms
%                   an array of frequency weights or
%                   a scalar representing the highest frequency term to keep.  
%                       note:  terms=0 keeps DC, terms=1 keeps DC & 1 cycle, terms=2 keeps DC + first 2 cycles, etc.
%                       note:  because of ambiguity, terms == 1 means 'DC + first term'.
%                       note:  if terms is non-scalar, it can either be a vector (applied to all fft sets)
%                              or matrix with flags unique to each set.
%                       note:  terms only needs to be as long as needed to specify all terms desired.  Will be filled
%                                   with zeros as needed.
%
    if (length(size(Y))>2)
        throw(MException('ICSF:IFFT_FILT','input dimension too high.  Can only operate on vectors or 2-D matrices'));
    end
      
    if (~exist('terms','var'))
        help(mfilename('fullpath'));
        yout=[];
        return;
    end

    if (~exist('dim','var') || isempty_s(dim))
        dim = find(size(Y)>1,1);
    end
    len=size(Y,dim);

            % create mask(s) to keep requested terms, zero out all other terms
            % then zero out the unwanted terms and inverse-FFT.
    
    if (isscalar(terms))        % if scalar, then keep all terms up to single value given in terms
        if (terms<0)
            terms=floor(len/2);
        end
        terms=ones((terms+1),1);
    end

    [nrows,ncols] = size(Y);
            
        % if terms-to-keep were supplied, turn them into a mask to zero out terms not wanted.
%     if (sum(terms ~= 0) ~= sum(terms == 1))
%         terms(terms ~= 0) = 1;      % make sure terms are all zeros or ones.
%     end
    if (isvector(terms))            % if given a vector, make sure it's the right type (col or row), 
                                    % then replicate it into a matrix if needed.
        if (dim==1)
            if (isrow(terms))       % could be wrong orientation.  fix it!
                terms=terms.';
            end
            if (ncols > 1 && isvector(terms))
                terms=repmat(terms, 1,ncols);
            end
        else
            if (iscolumn(terms))       % could be wrong orientation.  fix it!
                terms=terms.';
            end
            if (nrows > 1 && isvector(terms))
                terms=repmat(terms, nrows,1);
            end
        end
    end

    nt = size(terms, dim);               % terms is now matrix of positive-frequency terms to keep.
    nzeros = max(0, len - (2*nt-1));     % now clone it to the negative frequencies.
    if (dim==1)
        terms=[terms;zeros(nzeros,ncols);flipud(terms(2:end,:))];
    else
        terms=[terms,zeros(nrows,nzeros),fliplr(terms(:,2:end))];
    end
    Yout = Y .* terms;
    yout = ifft(Yout);
    
end
