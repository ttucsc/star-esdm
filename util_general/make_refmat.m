function R = make_refmat(lats, lons, nlats, nlons, dlat, dlon)
% creates a matlab GeoReference object for drawing maps from gridded data.
% note that latrange and lonrange should actually start and end at latbnds and lonbnds, so we generate that here by
% adding/subtracting 1/2 the lat or lon step.
%
%   Inputs:
%       lats    vector of latitudes such as read from netcdf file
%       lons    vector or longitudes such as read from netcdf file
%     (optional)
%       nlats   # of latitudes to draw.  [length(lats)], or, if lats is simply [lat1, lat2], calculated using dlat
%       nlons   # of longitudes to draw  [length(lons)], or, if lons is simply [lon1, lon2], calculated using dlon
%       dlat    grid latitude spacing [calculated from lats, nlats]
%       dlon    grid longitude spacing [calculated from lons, nlons]

    if (~exist('nlats','var') || isempty(nlats))
        if (~exist('dlat','var') || isempty(dlat))
            nlats = length(lats);
        else
            nlats = round((lats(end)-lats(1))/dlat) + 1;
        end
    end
    
    if (~exist('nlons','var') || isempty(nlons))
        if (~exist('dlat','var') || isempty(dlat))
            nlons = length(lons);
        else
            nlons = round((lons(end)-lons(1))/dlon) + 1;
        end
    end
    
    if (~exist('dlat','var') || isempty(dlat))
        dlat=0;
        dlon=0;
    end
    
    latrange = [lats(1)-dlat/2, lats(end)+dlat/2];
    lonrange = [lons(1)-dlon/2, lons(end)+dlon/2];
    latrange = max(-90,min(90,latrange));
    if (lonrange(1) < 0)
        lonrange = max(-180,min(180,lonrange));
    else
        lonrange = max(0,min(360,lonrange));
    end
    R = georefcells(latrange, lonrange, [nlats, nlons]);
end
    