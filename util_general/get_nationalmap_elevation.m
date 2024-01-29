function elev = get_nationalmap_elevation(lat,lon)
%   queries us national map (gov. svc) for elevation of location lat,lon
%   returns elevation in meters.
%   if unsuccessful (lat/lon not on map) returns nan.
%   otherwise returns elevation
%
%   For info on querying nationalmap.gov for elevation, see https://nationalmap.gov/epqs/


    url="https://nationalmap.gov/epqs/pqs.php";
    npts = length(lat);
    if (~exist('lon','var'))
            z=lat;
        if (isrow(z) || iscolumn(z))
            lat=z(1:2:end);
            lon=z(2:2:end);
        else
            lat=z(:,1);
            lon=z(:,2);
        end
        npts=length(lat);
    end
    elev = nan(npts,1);
    fixers=lon>180.0;
    lon(fixers)=lon(fixers)-360;
    for i=1:npts
        query=sprintf("?x=%f&y=%f&units=Meters&output=json",lon(i),lat(i));
        uri=sprintf('%s%s',url,query);
        resp=urlread(uri);
        e=jsondecode(resp);
        try
            el = e.USGS_Elevation_Point_Query_Service.Elevation_Query.Elevation;
            if (isnumeric(el)), elev(i)=el; end
        catch
        end
    end
end

