function [lix, qlls, end_flag, max_ix] = ll_region(llrange, lls, closure_flags, is_lons, wrap_ok)

% returns indexes and values of of lats or lons within the region defined by llrange.
% Assumes
%       lats are in ascending order -90 to 90               (do not need to cover entire range)
%       lons are in ascending order 0-360 or -180 to 180    (do not need to cover entire range)
%
%   Inputs:
%       llrange        range of latitudes or longitudes to bound, as in [-45,67]; or  [210,330];
%                                     to cross meridian:       [350,20]  will return data for range [350-360, 0-20]
%                           can also be vectors of lats or lons;  only 1st and last values of llrange are used.
%                           but should  be in sequential order (but may cross from 359 to 0)
%       lls             lats or lons to select bounded region from.
%       closure_flags   1- or 2-element vector specifying how to manage closure.  If 1 value, applies to both top and
%                           bottom of range.  If 2 values, 1st defines closure at bottom, 2nd defines closure at top.
%                       -2, -1, 0 or 1 or 2. absolutely inside, up to, exactly, at least, or absolutely outside  
%                           -2:  absolutely inside:  largest range completely inside lat or long range  
%                           -1:  up to:              largest range, potentially with closed boundaries (up-to-and-including)
%                            0:  exactly             lat or long must match latrange(1) and latrange(end) exactly
%                            1:  at least:           smallest range at least as large as range (potentially closed)
%                            2:  absolutely outside: smallest range completely containing range, open boundaries.
%       is_lons         boolean.  If true, lls are longitude values (0-360 or -180 - 180);  false:  lls are latitude values (-90 - 90)
%       wrap_ok         boolean, & optional.  Defaults to is_lons (wrap for longitude, don't wrap for latitude).  
%                                   If true, allow values to wrap (e.g., for longitude to cross meridian).
%                                   if false, truncate at top & bottom values.
%
%           NOTE:  comparisons for closure flags are "inexact", using a fudge factor.
%                   if lats/lons and range are both doubles, the fudge is about 1.1 cm
%                   if either lats/lons or range are float, the fudge is about 1.1 m
%

%   Outputs
%       lix           indexes of lats or lons in region  (not just 1st & last, but all indexes for region )
%       end_flag      true if output range crosses meridian (180 to -180 or 360 to 0) for lons, or
%                           is above top latitude or below bottom latitude  
%       qlls           vector of lats or lons in the region
%       maxix          index (in the output) of max qlls.  This is the point where the output longitudes wrap.
%                           example:  if qlls are [357, 358, 359, 0, 1, 2] , then maxix will be 3
%
%   if llrange(1) > llrange(2), then selects the region that crosses the 360 -> 0 degrees,
%       and sets end_flag to true.
%
%   For latitudes, if llrange goes above or below the values in lls, data stops at top (or bottom) value, unless wrap_ok
%   is set to true.  Then it wraps, just as longitudes do.  This would be useful to select arctic regions only.
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
%   modified 7/6/22      icsf to make sure returned lats,lons and ranges are doubles.
%   modified 2022-09-19  icsf matches now have tolerances so disconnects between single & double precision works properly.
%                               fixed bug where single/double could miss exact comparisons and return wrong range.


    close1 = closure_flags(1);      % -2 absolutely inside, -1 inside or equal, 0 exactly equal, 1 outside or equal, 1 absolutely outside
    close2 = closure_flags(end);
    
        % make sure lats & lons and ranges are doubles.
        
            % problems with exact comparisons when comparing floats and doubles.
            % so we set a fudge factor so our comparisons don't make an error when numbers aren't exact.
            %  For lat & lon, 1e-5 is about 1.1m, and 1e-7 is 1.1 cm.   
    if (all(strcmp("double", [class(llrange), class(lls)])))
        fudge = 1e-7;
    else
        fudge = 1e-5;
    end
    llrange = double(llrange);
    lls = double(lls);
    
    if (~exist('wrap_ok','var') || isempty(wrap_ok)), wrap_ok = is_lons; end
    
    latlons = to_row(lls);
    npts = length(lls);
    if (is_lons)
%       llrange = mod(llrange,360);
        latlons = [latlons-360, latlons, latlons+360];
        if (llrange(1) > llrange(end))
            llrange(end) = llrange(end) + 360;
        end
    else
%       llrange = mod(llrange+90, 180) - 90;
        latlons = [latlons-180, latlons, latlons+180];
        if (llrange(1) > llrange(end))
            llrange(end) = llrange(end) + 180;
        end
    end
    
    llrange1 = llrange(1);
    llrange2 = llrange(end);
    
    if (close1 == -2)       
        ix1 = find(latlons > llrange1+fudge, 1);  % absolutely inside
    elseif (close1 == -1)
        ix1 = find(latlons >= llrange1-fudge, 1);  % up-to and including
    elseif (close1 == 0)
        ix1 = find(abs(latlons-llrange1) <= fudge,1); % exactly equal
    elseif (close1 == 1)
        ix1 = find( latlons <= llrange1+fudge, 1, 'last'); % at least
    elseif (close1 == 2)
        ix1 = find( latlons < llrange1-fudge, 1, 'last'); % absolutely outside
    else
        error('error:  bad closure flag:  %f', close1);
    end
    
    if (close2 == -2)      
        ix2 = find(latlons < llrange2 - fudge, 1,'last');  % absolutely inside
    elseif (close2 == -1)
        ix2 = find(latlons <= llrange2 + fudge, 1, 'last'); % up-to and including
    elseif (close2 == 0)
        ix2 = find(abs(latlons - llrange2) <= fudge,1); % exactly equal
    elseif (close2 == 1)
        ix2 = find( latlons >= llrange2 - fudge, 1); % at least
    elseif (close2 == 2) 
        ix2 = find( latlons > llrange2 + fudge, 1); % absolutely outside
    else
        error('error:  bad closure flag:  %f', close1);
    end
    
    if (isempty(ix1) || isempty(ix2))
        lix = [];
        qlls = []; 
        end_flag = true;
    else
        if (ix1 > ix2)      % will happen if all points fall within same step inside the 
            tix=ix2;
            ix2=ix1;
            ix1=tix;
        end
        lix = (ix1:ix2) - npts; 
        if (wrap_ok)
            lix = mod((lix-1), npts) + 1;
        else
            lix = unique(max(1,min(npts, lix)));
        end

        qlls = lls(lix);
        
        end_flag = lix(1) > lix(end);
    end    
    
    if (~isempty(qlls)), max_ix = find(qlls == max(qlls)); end
end
    