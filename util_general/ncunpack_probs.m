function probs = ncunpack_probs(packed_probs, offset, scaling, nbytes)
% packs probabilities into an integer value, by first coverting the probabilities to a z-score, 
% then packing the z-score.  This will compress the central probabilities while giving 
        % check the input
    if (nargin < 3), error('ncunpack_probs:  error:  missing arguments'); end
    if (~exist('nbytes','var') || isempty(nbytes)), nbytes=8; end
    
    if (nbytes ~= 4 && nbytes ~= 8)
        error('ncunpack_probs:  bad nbytes (%d):  must be 4 or 8');
    end
    
    zscores = ncunpack(packed_probs, nbytes, true, offset, scaling);
    
    probs = normcdf(zscores);
    
end
