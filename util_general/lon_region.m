function [lonix, qlons, meridian_flag, max_ix] = lon_region(lonrange, lons, closure_flags, wrap_ok)
% function to return vector of lons within lonrange (subject to closure_flags).  
% This is just a wrapper to use ll_region(...).
% See comments in ll_region(...) for explanation of inputs and outputs.
    if (~exist('wrap_ok','var')), wrap_ok = true; end
    
    [lonix, qlons, meridian_flag, max_ix] = ll_region(lonrange, lons, closure_flags, true, wrap_ok);
end
