function seed = latlon_seed(lat, lon)
    % generates a random seed from lat & lon
    % returns an integer between 0 & 2^32 - 1
    %   This should return a unique integer number for locations at least 100 cm apart.

    lat = mod(lat+90,180);
    lon = mod(lon, 360);
    lat_seed = floor((lat)/180 * 65534);
    lon_seed = floor(lon/360 * 65534);
    
    seed = mod(floor(lat_seed + 64000*lon_seed + lat*10000 + lon*10000), 2^32); 
end

