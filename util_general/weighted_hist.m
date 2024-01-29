function whist = weighted_hist(s, weights, fieldname)
%   returns yrlen weighted histograms of s or of s.(fieldname), if fieldname is specified.
%   Inputs
%       s           A cell array of DA objects, with histograms in field fieldname
%                       s.(fieldname) can be a single histogram or a matrix of histograms.
%       weights     weights to use for summing histograms, 1 for each element of s..  
%                       Weights do not need to sum to 1... they will be normalized to sum to 1.
%                       so one can use a vector of all 1's to use even weighting.
%       fieldname   field within s containing histograms to be summed.
%
%   Returns:
%       whist       weighted sum of the histograms

    n = numel(s);
    
    if (isempty(weights) || numel(weights)==1)
        weights = ones(1,n);
    end
    weights = weights/sum(weights);     % normalize weights to 1.
 
    for i=1:n
        if (i==1)
            whist = weights(1) * s{i}.(fieldname);
        else
            whist = weights(i) * s{i}.(fieldname) + whist;
        end
    end
end

