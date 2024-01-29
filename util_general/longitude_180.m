function lon = longitude_180(lon)
% function lon = lon_180_wrap(lon)
%   returns longitude adjusted to range -180 - +179.999999....
    lon = mod(lon+180,360)-180;
end

