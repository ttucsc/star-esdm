function [latix, qlats, end_flag, max_ix] = lat_region(latrange, lats, closure_flags, wrap_ok)
% function to return vector of lats within latrange (subject to closure_flags).  
% This is just a wrapper to use ll_region(...).
% See comments in ll_region(;;;) for explanation of inputs and outputs.
    if (~exist('wrap_ok','var')), wrap_ok = false; end
    
    [latix, qlats, end_flag, max_ix] = ll_region(latrange, lats, closure_flags, false, wrap_ok);
end
