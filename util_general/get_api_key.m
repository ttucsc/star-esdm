function apikey = get_api_key(api, noerror)
%   reads file api_keys.csv for a specific key to make a web query for information such as timezone, elevation, etc.
%
%   Possible Keys:
%       key                 data source
%       ===                 ===========
%       timezoneDB          timezoneDB.com
%       time_zone           google maps timezone web service  https://developers.google.com/maps/documentation/timezone/intro
%       elevation           google maps elevation web service https://developers.google.com/maps/documentation/elevation/intro
%                               elevation also available from https://nationalmap.gov/epqs/ (no key required, but slower than google)
%
%   returns a struct containing 3 fields:
%       key         API key from web-query API source
%       email       email registered with key
%       API         name of API.

    keytbl = readtable('api_keys.csv');
    ix = find(strcmpi(keytbl.API, api), 1);
    if (isempty(ix))
        if(exist('noerror','var') || noerror)
            apikey=[];
        else
            error('get_api_key:  key %s does not exist.  Type help get_api_key for more info', api);
        end
    else
        apikey=keytbl.key{ix};
    end
end

