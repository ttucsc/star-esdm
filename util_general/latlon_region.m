function [latix, lonix, lon_meridian_flag, qlats, qlons, lat_end_flag, lat_max_ix, lon_max_ix] = latlon_region(latrange, lonrange, lats, lons, closure_flags, wrap_ok_flags)
% function to return vectors of lats lons within latrange and lonrange (subject to closure_flags).  

% returns indexes and values of of lats& lons within the region defined by latrange and lonrange.
% Assumes
%       lats are in ascending order -90 to 90               (do not need to cover entire range)
%       lons are in ascending order 0-360 or -180 to 180    (do not need to cover entire range)
%
%   Inputs:
%       l??range        range of latitudes && longitudes to bound, as in [-45,67]; or  [210,330];
%                                   to cross meridian:       [350,20]  will return data for range [350-360, 0-20]
%                                   lonrange can be in range -180-180 even if lons are 0-360 && vice-versa.
%                                       qlons will be in same range as input lons.
%                           can also be vectors of lats or lons;  only 1st and last values of l??range are used.
%                           but should  be in sequential order (but may cross from 359 to 0)
%       l??s            lats && lons to select bounded region from.
%       closure_flags   1- or 2-element vector specifying how to manage closure.  If 1 value, applies to both top and
%                           bottom of range.  If 2 values, 1st defines closure at bottom, 2nd defines closure at top.
%                           if only 1 row, applies to both lat & lon.  If 2 rows, 1st is for lats, 2nd for lons.
%                       -2, -1, 0 or 1 or 2. absolutely inside, up to, exactly, at least, or absolutely outside  
%                           -2:  absolutely inside:  largest range completely inside lat or long range  
%                           -1:  up to:              largest range, potentially with closed boundaries (up-to-and-including)
%                            0:  exactly             lat or long must match latrange(1) and latrange(end) exactly
%                            1:  at least:           smallest range at least as large as range (potentially closed)
%                            2:  absolutely outside: smallest range completely containing range, open boundaries.
%       wrap_ok_flags   boolean, & optional.  Defaults to wrap for longitude, don't wrap for latitude).  
%                                   If true, allow values to wrap (e.g., for longitude to cross meridian).
%                                   if false, truncate at top & bottom values.
%                                   NOTE:  If only 1 element, applies to both lats & lons!!!
%

%   Outputs
%       l??ix         indexes of lats in region  (not just 1st & last, but all indexes for region )
%       lon_meridian_flag
%       lat_end_flag  true if output range crosses meridian (180 to -180 or 360 to 0) for lons, or
%                           is above top latitude or below bottom latitude  (lats have 
%       ql??s          vector of lats or lons corresponding to latix and lonix
%       l??_max_ix     index (in the output) of max ql??s.  This is the point where the output longitudes wrap.
%                           example:  if qlons are [357, 358, 359, 0, 1, 2] , then lon_max_ix will be 3
%
%   if lonrange(1) > lonrange(2), then selects the region that crosses the 360 -> 0 degrees,
%       and sets end_flag to true.
%
%   For latitudes, if latrange goes above or below the values in lls, qlats stops at top (or bottom) value, unless wrap_ok
%   is set to true.  Then it wraps, just as longitudes do.  This would be useful to select arctic && antarctic regions
%   together and exclude midlle latitudes.
%
%   If using the output to interpolate, and meridian_flag is true, you can wrap the longitudes to the range
%   [-180,180] with:
%       qlons = mod(qlons+180,180) - 180;
%
%   When using output to find the correct gridpoints to use for a lon:
%       if (lonpt >= max(qlons) || lonpt < min(qlons))
%           ix1 = find(qlons==max(qlons));
%           ix2 = find(qlons==min(qlons));
%       else
%           if (lonpt >= qlons(1))
%               ix1 = find(qlons(1:lon_max_ix) <= lonpt,1,'last');      % find west boundary in 1st part of list
%           else
%               ix1 = max_ix + find(qlons(maxix+1:end) <= lonpt,1,'last');  % fine west boundary in 2nd part of list
%           end
%           ix2 = mod(ix1,npts)+1;  % make sure we wrap back to 1st point, in case ix1 is last point.
%       end
%

%       % earlier version of this code took a single boolean for closure flags.  true -> [1, 1];  false -> [2, 2]
    if (islogical(closure_flags) && length(closure_flags)==1)   
        if (closure_flags)
            closure_flags = [1,1];
        else
            closure_flags = [2,2];
        end
    end
    if (~exist('wrap_ok_flags','var')), wrap_ok_flags = [false, true]; end      % don't wrap lats, do wrap lons.
        
    [latix, qlats, lat_end_flag,      lat_max_ix] = ll_region(latrange, lats, closure_flags(1,:),   false, wrap_ok_flags(1));
    [lonix, qlons, lon_meridian_flag, lon_max_ix] = ll_region(lonrange, lons, closure_flags(end,:), true, wrap_ok_flags(end));
end
