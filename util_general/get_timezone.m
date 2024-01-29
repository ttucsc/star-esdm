function [tz, dst, tzName, dst_start, dst_end] = get_timezone(lat,lon, tstamp, key)
%returns timezone info from lat/lon
%   Data is obtained via http query to timzoneDB.com 
%   
%   Inputs:
%       lat,lon         column vectors of lat/lon locations
%       tstamp          either:
%                           1 datenum (applies to all locations)
%                           1 datevec ([yyyy,mm,dd]) (applies to all locations)
%                           col. vector of datenums, 1 per location
%                           matrix of size n x 3 (or n x 6) w/ datevec for each location.
%       key             timzoneDB user key.
%
%       note:  timezone info is also available via Google Maps API.  https://developers.google.com/maps/documentation/timezone/intro
%   this code queries timezoneDB for timezone info.


    if (~exist('tstamp','var') || isempty(tstamp))
        tstamp=posixtime(datetime(datevec(now)));
    elseif (~isdatetime(tstamp))
        if (length(tstamp)==1)
            tstamp=posixtime(datetime(datevec(tstamp)));
        else
            tstamp=posixtime(datetime(tstamp));
        end
    end
        
% icsf modified from:
%GETELEVATIONS queries Google Maps API webservice for ground elevations.
%
%   elevation = GETELEVATIONS(latitude, longitude) returns ground 
%   elevation for latitude and longitude arrays.
%
%   [elevation, resolution] = getElevations(latitude, longitude, ...
%     'key', 'AIzaSyCuN8tjAVEaXorgNjS1tDiVC-oc0QBJoYc' );
%   is an example of a call passing additional attributes to Google Maps 
%   API webservice, and capturing also the resolution of the data. 
%   See https://developers.google.com/maps/documentation/elevation/
%   for details.
%
% Author: Jarek Tuszynski (jaroslaw.w.tuszynski@leidos.com)
% License: BSD (Berkeley Software Distribution)
% Documentation: https://developers.google.com/maps/documentation/elevation/



% Check inputs
    nPos = numel(lat);
    assert(nPos>0, 'lat and longitude inputs can not be empty')
    assert(nPos==numel(lon), 'lat and lon inputs are not of the same length')
    assert(min(lat(:)) >= -90 && max(lat(:)) <= 90, 'lats has to be between -90 and 90')
    assert(min(lon(:))>=-180 && max(lon(:))<=180, 'lon has to be between -180 and 180')

% Query timezoneDB
    tz  = nan(size(lat));
    dst = nan(size(lat));
    dst_start = nan(size(lat));
    dst_end   = nan(size(lat));
    tzName = strings(size(lat));
    
    for i=1:length(lat)

        website = 'http://api.timezonedb.com/v2/get-time-zone?';
        keyStr = sprintf('key=%s', key);
        fmt='&format=xml';
        by='&by=position';
        coord=sprintf('&lat=%f&lng=%f',lat(i),lon(i));
        if (nargin > 2)
            ts=sprintf('&timestamp=%d',floor(tstamp));
        else
            ts=[];
        end
        fields='&fields=dst,gmtOffset';
        if (nargout > 2)
            fields=sprintf('%s,%s', fields,'zoneName');
            if (nargout>3)
                fields=sprintf('%s,%s', fields,'dstStart,dstEnd');
            end
        end

      % create query string and run a query
        url = [website, keyStr, fmt, by, coord, ts, fields];
        str = urlread(url);

      % Parse results
        status = regexp(str, '<status>([^<]*)<\/status>', 'tokens');
        switch status{1}{1}
        case 'OK'
            res = regexp(str, '<gmtOffset>([^<]*)<\/gmtOffset>', 'tokens');
            tz(i) = str2double(res{1})/3600;
            if (nargout>1)
                res = regexp(str, '<dst>([^<]*)<\/dst>', 'tokens');
                dst(i) = str2double(res{1});
            end
            if (nargout>2)
                res = regexp(str, '<zoneName>([^<]*)<\/zoneName>', 'tokens');
                tzName(i) = string(res{1});
            end
            if (nargout>3)
                res = regexp(str, '<dstStart>([^<]*)<\/dstStart>', 'tokens');
                dst_start(i) = datenum([1970,1,1,0,0,0]) + str2double(res{1})/3600/24;
            end
            if (nargout>4)
                res = regexp(str, '<dstEnd>([^<]*)<\/dstEnd>', 'tokens');
                dst_end(i) = datenum([1970,1,1,0,0,0]) + str2double(res{1})/3600/24;
            end
        case 'INVALID_REQUEST'
            error('Google Maps API request was malformed');
        case 'OVER_QUERY_LIMIT'
            error('Google Maps API requestor has exceeded quota');
        case 'REQUEST_DENIED'
            error('Google Maps API did not complete the request (invalid sensor parameter?)');
        case 'UNKNOWN_ERROR'
            error('Google Maps API: an unknown error.');
        end
        if (length(lat)>1), pause(1); end
    end
end


% https://maps.googleapis.com/maps/api/timezone/json?locations=39.7391536,-104.9847034&timestamp=1331766000&key=

% google developers account:  kanesf@gmail.com

% api key:  google_latlon_api_key.csv