function distrib = calc_day_mapping_distrib(d, yrlen, method, nterms, sig_terms)
%   returns a set of values for remapping the time dimension of histograms to a different number of days
%   Same as is used in ARRM_V2_disaggregate_signal
%
%   methods:
%       linear      % for simply changing the year length to something that runs faster for fourier transforms (yrlen of 128 is good)
%       cos         % uses cosine(slope(climatology)).
%       clim        % uses climatology.  Good for rescaling precip to get non-constant smoothing.
%

    [ndays, nsets] = size(d);

    if (strncmpi(method, 'lin',3))
        distrib = ones(nsets,yrlen);
    else    
        distrib = zeros(nsets,yrlen);
        for j=1:nsets
            clim = climatology(d(:,j), nterms, sig_terms,ndays);
            if (strncmpi(method, 'cos',3))
                slopes = circular_slopes(clim);
                distrib(j,:) = 1+sin(abs(slopes)*pi/4);
            elseif (strncmpi(method, 'clim',4))
                distrib(j,:) = clim;
            end
        end 
    end
end            


