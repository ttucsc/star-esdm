function lon = longitude_360(lon)
% function lon = lon_180_wrap(lon)
%   returns longitude adjusted to range -180 - +179.999999....
    lon = mod(lon,360);
end

