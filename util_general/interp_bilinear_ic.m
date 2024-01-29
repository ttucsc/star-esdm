function [ d_out, wts ] = interp_bilinear_ic( latpt, lonpt, lats, lons, d11, d21, d12, d22 )
%       Linear interpolation of point within rectangular grid.
%       interpolates for all timesteps at same lat & lon, 
%       Correctly handles the case where the longitude wraps (such from 359 to 0).
%           (via fix_range(...) function below).
%
%   Note:  inputs are:  y-values first, x-values second -- (lat, lon):
%           y, x, y1, x1, y2, x2, vals(y1,x1), vals(y1,x2), vals(y2,x1), vals(y2,x2);
%   Note2:  equations usually written in terms of (x,y), and when done so,
%               Tlat1,lon2) is the value at (x2, y1)
%               Tlat2,lon1) is the value of (x1, y2)
%                   watch indexes carefully when comparing to, e.g., wikipedia page
%                   https://en.wikipedia.org/wiki/Bilinear_interpolation
%
%   Inputs:
%       location:
%           latpt, lonpt        lat, lon of location
%           lat1, lon1          coords of one corner of 4 gridcells  (usually bottom left)
%           lat2, lon2          coords of opposite corner (usually upper right)
%       data tp omter[p;ate:
%         either:
%           data                npts x 2 x 2 matrix  where 2x2 grid is (2 lats, 2lons)
%         or
%           data                2x2 set of single points, (2 lats, 2lons)
%         or
%           Tlat1lon1           value at lat1, lon1
%           Tlat2lon1           value at lat2, lon1
%           Tlat1lon2           value at lat1, lon2
%           Tlat2lon2           value at lat2, lon2
%
%   Outputs:
%       Tout                output value or output series.
%       wts                 weights used for interpolation.

    wts = bilinear_weights_ic(latpt, lonpt, lats, lons);

%     fprintf('interp_bilinear weights: %8.6f %8.6f %8.6f %8.6f\n', w11, w12, w21, w22); 

    nd = ndims(d11);
    [nd1, nd2, nd3] = size(d11);
    if (ndims(d11)==3 && nd2==2 && nd3 == 2)
        d_out = w(1,1)*d11(:,1,1) + w(2,1)*d11(:,2,1) + w(1,2)*d11(:,1,2) + w11(2,2)*d11(:,2,2);
    elseif (nd==2 && nd1==2 && nd2 == 2)
        d_out = w(1,1)*d11(1,1) + w(2,1)*d11(2,1) + w(1,2)*d11(1,2) + w(2,2)*d11(2,2);
    else
        d_out = wts(1,1)*d11 + wts(2,1)*d21 + wts(1,2)*d12 + wts(2,2)*d22; 
    end

end


