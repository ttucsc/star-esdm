function [packed_probs, nanval, offset, scaling] = ncpack_probs(probs, minz, maxz, nbytes, offset, scaling)
% packs probabilities into an integer value, by first coverting the probabilities to a z-score, 
% then packing the z-score.  
%
%   minz, maxz of -7.5 and 7.5 is good for nbutes=2;

        % check the input
    if (nargin < 3), error('ncpack_probs:  error:  missing arguments'); end
    if (~exist('nbytes','var') || isempty(nbytes)), nbytes=2; end
    if (minz >= maxz), error('ncpack_probs: bad z-range'); end
    
    if (nbytes == 1)
        maxval = 254;
    elseif (nbytes == 2)
        maxval = 65534;
    elseif (nbytes == 4)
        maxval = 2^32-2;
    else
        error('ncpack_probs:  bad nbytes (%d):  must be 1, 2 or 4');
    end
    
        % set any probs outside of range to nan's. 
    probs(probs < normcdf(minz)) = nan;
    probs(probs > normcdf(maxz)) = nan;
       
    z_score = norminv(probs);
    
    if (~exist('offset','var')  || isempty(offset)),  offset  = minz;                 end
    if (~exist('scaling','var') || isempty(scaling)), scaling = (maxz - minz)/maxval; end
    
    [packed_probs, nanval, offset, scaling] = ncpack(z_score, nbytes, true, offset, scaling);
    
end

