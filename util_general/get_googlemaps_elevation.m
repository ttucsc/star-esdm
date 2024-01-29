function [elev, resolution, key] = get_googlemaps_elevation(lat,lon, do_all, key)
%returns elevation from lat/lon
%   Data is obtained via http query to Google Maps API.
%
%   For google query to work, you need an API key from google elevation API services
%   and put it in a file called api_keys.com
%   File should have 3 columns:
%           key         key (provided from API source)
%           email       email registered to key
%           API         name of API for ARRM_V2's web-based query code
%
%   for info on keys and queries, see:
%       key                 data source
%       ===                 ===========
%       elevation           google maps elevation web service https://developers.google.com/maps/documentation/elevation/intro
%                               elevation also available from https://nationalmap.gov/epqs/ (no key required)

    if (~exist('lon','var') && size(lat,2)==2)
        lon=lat(:,2);
        lat=lat(:,1);
    end
    if (~exist('do_all','var') || isempty(do_all) || ~do_all)
        do_all=false;
    else
        do_all=true;
    end
    
    lon = mod(lon+180,360)-180;     % lons must be in range [-180,180);
    
    if (~exist('key','var') || isempty(key))
        key=get_api_key('elevation', true);
        if (isempty(key))
            error('get_googlemaps_elevation:  No API key found.  You need a google maps API key for this.  Type "help get_googlemaps_elevation" for more info');
        end
    end
    
    n=numel(lat);
    if (do_all || n <=10)
        [elev,resolution] = getElevations(lat,lon,'key',key);
    else
        elev=nan(size(lat));
        resolution=nan(size(lat));
        for i=1:10:n
            ix=(0:9)+i;
            ix(ix>n)=[];
            [elev(ix),resolution(ix)] = getElevations(lat(ix),lon(ix),'key',key);
        end
    end
            
    
%     locations=sprintf("%.4f,%.4f",lat(1),lon(1));
%     npts=size(lat,1);
%     for i=2:npts
%         ss=sprintf("|%.4f,%.4f",lat(i),lon(i));
%         if (strlength(ss)+140 > 8192), break; end
%         locations=locations+ss;
%     end
%     
%     uri=sprintf("https://maps.googleapis.com/maps/api/elevation/xml?locations=%s&key=%s",locations,key);
%     
%     
%     str=webread(uri);
%     
%     status = regexp(str, '<status>([^<]*)<\/status>', 'tokens');
%     switch status{1}{1}
%         case 'OK'
%             res = regexp(str, '<elevation>([^<]*)<\/elevation>', 'tokens');
%             elev = cellfun(@str2double,res);
%             if (nargout>1)
%                 res = regexp(str, '<resolution>([^<]*)<\/resolution>', 'tokens');
%                 resolution = cellfun(@str2double,res);
%             end
%         case 'INVALID_REQUEST'
%             error('Google Maps API request was malformed');
%         case 'OVER_QUERY_LIMIT'
%             error('Google Maps API requestor has exceeded quota');
%         case 'REQUEST_DENIED'
%             error('Google Maps API did not complete the request (invalid sensor parameter?)');
%         case 'UNKNOWN_ERROR'
%             error('Google Maps API: an unknown error.');
%     end       
end

% https://maps.googleapis.com/maps/api/elevation/json?locations=39.7391536,-104.9847034&key=

% google developers account:  iancsf@gmail.com

% api key:  google_latlon_api_key.csv