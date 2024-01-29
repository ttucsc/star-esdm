function [latix, lonix, latout, lonout] = closest(lat_pts, lon_pts, lats, lons, nlat, nlon)
%   Returns indexes of closest points in matrices lons & lats to the points (lat_pts(i), lon_pts(i))
%     Inputs:
%       lon_pt, lat_pt  lon, lat of point to locate 
%       lons, lats      matrices of lons & lats values
%       nlon, nlat      # of points (indexes) to return [1, 1]
%                           for bilinear interpolation, use [2,2].
%     Outputs:
%       lonix, latix    indexes of points closest to lons & lats
%                       indexes are 1-based.  Subtract 1 to use for netcdf files, which use 0-based arrays.
%       lonout, latout  values of lon & lat arrays corresponding to lonix and latix   
%
%   Note:  always returns nlon & nlat values.
%           for latitude, if lat1 is below lats(1) or above lats(end), returns 1st or last nlat indexes and values
%                           For points below 1st lat, returns 1st lats; or above last lat, returns last lats.
%                           So bilinear interpolation will be extrapolating from beginning or end.
%                              Perhaps should return points 180 degrees away in longitude instead?
%           for longitude, it wraps around at the end, which is presumed to be just below 360 degrees.
%                           if lons don't cover full range of 0-360, results will not be useful.
%               Thus, if trying to interpolate points near prime meridian, make sure lons start near 0 and end near 359. 
%                           if lons cover prime meridian properly, then it returns correct values for locations between
%                           last lon. value and first lon. value.  e.g., returns 359 & 0 for lon_pt 359.5. 
%
%   If all lat_pts or all lon_pts fall exactly on a lat or lon and nlat or nlon are even, then the output range is
%   semi-open at the top.  For example, if lon_pt is 5 and grid step is 1, then the returned lats and lons are [5,6),
%   not (4,5];
    npts = length(lat_pts);
    
    if (~exist('nlat','var') || isempty(nlat)), nlat = 1; end
    if (~exist('nlon','var') || isempty(nlon)), nlon = 1; end
    
    if (npts == 1)
        [latix, lonix, latout, lonout] = closest_sub(lat_pts, lon_pts, lats, lons, nlat, nlon);
    else        
        latix = zeros(npts, nlat);
        lonix = zeros(npts, nlon);
        latout = zeros(npts, nlat);
        lonout = zeros(npts, nlon);

        for i=1:npts
            [latix(i,:), lonix(i,:), latout(i,:), lonout(i,:)] = closest_sub(lat_pts(i), lon_pts(i), lats, lons, nlat, nlon);
        end
    end
end
    

    % does the work for a single (lat,lon) point.
function [latix, lonix, latout, lonout] = closest_sub(lat_pt, lon_pt, lats, lons, nlat, nlon)
%
%   Returns indexes of closest points in matrices lons & lats to the points lon1 and lat1.
%     Inputs:
%       lon_pt, lat_pt  lon, lat of point to locate 
%       lons, lats      matrices of lons & lats values
%       nlon, nlat      # of points (indexes) to return [2, 2]
%     Outputs:
%       lonix, latix    indexes of points closest to lons & lats
%                       indexes are 1-based.  Subtract 1 to use for netcdf files, which use 0-based arrays.
%       lonout, latout  values of lon & lat arrays corresponding to lonix and latix
%
%   Note:  always returns nlon & nlat values.
%           for latitude, if lat1 is below lats(1) or above lats(end), returns 1st or last nlat indexes and values
%           for longitude, it wraps around at the end, which is presumed to be 360 degrees.
%                           if lons don't cover full range of 0-360, results will not be useful.
%           Thus, if longitudes are not a full circle, make sure that all points are inside the longitude range! 

    if (~exist('nlon','var') || isempty(nlon)), nlon=1; end
    if (~exist('nlat','var') || isempty(nlat)), nlat=1; end

    % find indexes of closest lons & lats points.
    londist=abs(lons-lon_pt);
        % in case closest point is on far side of 360-degree wraparound for longitude:
    wrappers= londist>180.0;
    londist(wrappers) = abs(londist(wrappers)-360);      % wraps the distances
    lonix=find(londist==min(londist), 1);    

    latdist=abs(lats-lat_pt);
    latix=find(latdist==min(latdist), 1);
        % get closest lons,lats location.
    lonout=lons(lonix);
    latout=lats(latix);

        % if more than 1 lons or lats requested, get full range

    if (nlat > 1)
        nnlat = length(lats);
        dlat = floor(nlat/2);  % # above & below if odd # requested.
            % deal with odd, even # of points requested
        if (mod(nlat,2)==1)     % odd
            l1=latix-dlat;
            l2=latix+dlat;
        elseif (lat_pt>=latout)    % even, and above closest
            l1=latix-dlat+1;
            l2=latix+dlat;
        else                    % even, and below closest
            l1=latix-dlat;
            l2=latix+dlat-1;
        end
        latix = l1:l2;
            % make sure indexes are in range.  shift up or down if indexes out of range
        if (latix(1) < 1)
            latix = 1:length(latix);
        elseif (latix(end) > nnlat)
            latix = (nnlat-length(latix)+1):nnlat;
        end
        latout = lats(latix);
    end

    if (nlon > 1)
        nnlon = length(lons);
            % figure out whether we're above or below.  wrap if dist > 180.
        mylondist = lon_pt - lons(lonix);
        if (mylondist > 180)
            mylondist = mylondist-360;
        elseif (mylondist < -180)
            mylondist = mylondist+360;
        end
        dlon = floor(nlon/2);
        if (mod(nlon,2)==1)     % odd
            l1 = lonix-dlon;
            l2 = lonix+dlon;
        elseif (mylondist >= 0)  % even, and above closest
            l1 = lonix-dlon+1;
            l2 = lonix+dlon;
        else                    % even, and below closest
            l1 = lonix-dlon;
            l2 = lonix+dlon-1;
        end
            % make sure indexes are in range.  lons (longitude) wraps around, so we use modulus here
        lonix = mod((l1:l2)-1, nnlon) + 1;
        lonout = lons(lonix);
    end
end
