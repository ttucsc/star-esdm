function whist = weighted_hist_mat(hhh, weights)
%   returns 2-D weighted histogram surface of hhh
%   Inputs
%       hhh         matrix of histograms, of size nbins X ndays X nsets
%       weights     weights to use for histogramming.  One weight for each column.  
%                       Weights do not need to sum to 1... they will be normalized to sum to 1.
%                       so one can use a vector of all 1's to use even weighting.
%
%   Returns:
%       whist       weighted histogram of input data, of size ( yrlen, nbins )

    
    [~,~, nsets] = size(hhh);
    
    if (isempty(weights))
        weights = ones(1,nr)./nsets;
    end

    for i=1:nsets
        if (i==1)
            whist   = squeeze(hhh(:,:,i)) .* weights(i);
        else
            whist   = whist   + squeeze(hhh(:,:,i)) .* weights(i);
        end
    end
end

