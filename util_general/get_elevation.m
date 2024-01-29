function elev = get_elevation(lat,lon,do_google,do_all,google_apikey)
%query for online elevation info
%
%   lat, lon    vectors of locations to query for.
%   do_google   use googlemaps for query.  Limited to 2500/day.
%                   but faster for multiple queries than nationalmap
%                   Use googlemaps if locations are outside US.
%   do_all      query googlemaps for all elevations at once.  Faster, but less accurate.
%
%   By default, queries national map first, then goes to google for any points 
%   that national map returned NA for.
%
%
    if (~exist('lon','var') || isempty(lon))
        z=lat;
        if (isrow(z) || iscolumn(z))
            lat=z(1:2:end);
            lon=z(2:2:end);
        else
            lat=z(:,1);
            lon=z(:,2);
        end
    end
    
    
    if (nargin <3), do_google = false; end
    if (nargin < 4), do_all = false; end
    if (nargin < 5), google_apikey=[]; end
    fixers=lon>180.0;
    lon(fixers)=lon(fixers)-360;    
    if (do_google)
        elev = get_googlemaps_elevation(lat,lon,do_all,google_apikey);
    else
        elev = get_nationalmap_elevation(lat,lon);
        missing=find(isnan(elev));
        if (~isempty(missing))
            e2 = get_googlemaps_elevation(lat(missing),lon(missing),do_all,google_apikey);
            elev(missing) = e2;
        end
    end
end

