function [tz, dst, tzName, dst_start, dst_end, apikey] = getTZ(lat, lon, tstamp, apikey)
%   Function to return the timezone info for a set of lat/lon locations
%
%   For this to work, you need an API key from timezoneDB  https://timezonedb.com/
%   and put it in a file called api_keys.com
%   File should have 3 columns:
%           key         key (provided from API source)
%           email       email registered to key
%           API         name of API for ARRM_V2's web-based query code
%
%   for info on keys and queries, see:
%       key                 data source
%       ===                 ===========
%       timezoneDB          timezoneDB.com
%       time_zone           google maps timezone web service  https://developers.google.com/maps/documentation/timezone/intro


    if (~exist('tstamp','var') || isempty(tstamp)), tstamp=now(); end
    
    lon_fixers = lon > 180;
    lon(lon_fixers) = mod(lon(lon_fixers)+180,360) - 180;
    
    if (~exist('apikey','var') || isempty(apikey))
        apikey = get_api_key('timezoneDB', true);
        if (isempty(apikey))
            error('getTZ:  this code needs a timezoneDB API key.  type "help getTZ" for more information');
        end
    end
    
    if (~exist('tstamp','var') || isempty(tstamp))
        tstamp = now;
    end
    
    [tz,dst, tzName, dst_start, dst_end] = get_timezone(lat, lon, tstamp, apikey);
    
end
