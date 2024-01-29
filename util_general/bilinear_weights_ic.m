function wts  = bilinear_weights_ic( latpt, lonpt, lats, lons )
%       Returns weights for bilinear interpolation of  point within rectangular grid.
%       Correctly handles the case where the longitude wraps (such from 359 to 0, or from 180 to -180.)
%
%   This solution is not quite correct for spherical coordinates;  A quadrangular solution would be better, especially
%   near the poles.  See for example, https://www.researchgate.net/publication/331602694_An_Alternative_Bilinear_Interpolation_Method_Between_Spherical_Grids
%
%   Note:  inputs are:  y-values first, x-values second -- (lat, lon):
%           y, x, [y1, y2], [x1, x2];
%   Note2:  equations are usually written in terms of (x,y), which is reverse of this.
%                   watch indexes carefully when comparing to, e.g., wikipedia page
%                   https://en.wikipedia.org/wiki/Bilinear_interpolation
%                   Comments below show how to match the inputs & outputs with Wikipedia's equations.
%   NOTE3:  longitudes *MUST* be western 1st, eastern 2nd.  
%           If 1st > 2nd, then code assumes longitudes crosses the 180->-180 wrapping point 
%                                                           or the 360->0 wrapping point.
%
%   Inputs:
%       latpt, lonpt        lat, lon of location
%       lats                2-element vector of latitudes  (should be [south-er, north-er]   (y1, y2)
%       lons                2-element vector of longitudes (should be [west-er, east-er]     (x1, x2)
%                               These coords define the four gridcell coordinates
%                                   [lats(1), lons(1)]  lower left  (takes w_y1_x1)
%                                   [lats(2), lons(1)]  upper left  (takes w_y2_x1)
%                                   [lats(1), lons(2)]  lower right (takes w_y1_x2)
%                                   [lats(2), lons(2)]  upper right (takes w_y2_x2)
% 
%   Outputs:
%       wts                 2x2 matrix [w_y1_x1, w_y1_x2;  w_y2_x1, w_y2_x2]
%                               w_y1_x1:  lower left weight
%                               w_y2_x1:  upper left weight
%                               w_y1_x2:  lower right weight
%                               w_y2_x2:  upper right weight
%                               note that viewed as a matlab 2x2 matrix, weights could be considered flipped up-down,
%                               if you're expecting the 1st row to be north, and the second south.
%                               1st row is bottom lat (southwest wt, southeast weight)
%                               2nd row is top    lat (northwest wt, northeast weight)
%
%                               This corresponds to the way the data is read out of our rotated netcdf files
%                               where the first matrix comes out (time, lat, lon), with time varying most quickly,
%                               followed by latitude (y), then longitude (x)

        % special cases:  if lats(1) & lats(2) are the same, or lons(1) and lons(2) are the same.
        % in that case, we just do standard linear interpolation instead of bilinear.
        % but split the weights between the 2 identical points.  This is to handle case
        % where the points are identical, but the data is not.  This weights each co-located
        % set evenly.
    wts=nan(2,2);
        % put lons in range -180 to 180
    lons=longitude_180(lons);
    lonpt = longitude_180(lonpt);
    if (lons(2) < lons(1))          % but if points cross the 180 -> -180 wrapping point,
        lons=longitude_360(lons);   % then switch longitudes to 0-360 range so code works for all cases.
        lonpt = longitude_360(lonpt);
    end
    if (lons(1) == lons(2) && lats(1) == lats(2))       % all 4 grid points are identical.
        wts = ones(2,2)/4;                              % simply weight each data set identically.

    elseif (lons(1) == lons(2))                   % lon pts are identical.  do simple linear interpolation
        denom = (lats(2) - lats(1));              % on latitude.
        wts(1,1) = (lats(2)-latpt)/denom/2;  
        wts(1,2) = wts(1,1);
        wts(2,1) = (latpt-lats(1))/denom/2;
        wts(2,2) = wts(2,1);
    elseif (lats(1) == lats(2))                   % lat pts are identical. do simple linear interpolation        
        denom = (lons(2) - lons(1));              % on longitude.
        wts(1,1) = (lons(2)-lonpt)/denom/2;
        wts(2,1) = wts(1,1);
        wts(1,2) = (lonpt-lons(1))/denom/2;
        wts(2,2) = wts(1,2);
    else                                    % four distinct points.
        denom = (lons(2)-lons(1))*(lats(2)-lats(1));
        wts(1,1) = (lons(2)-lonpt)*(lats(2)-latpt)/denom;   % Wikipedia's Q11  bottom left
        wts(2,1) = (lons(2)-lonpt)*(latpt-lats(1))/denom;   % Wikipedia's Q12  upper left
        wts(1,2) = (lonpt-lons(1))*(lats(2)-latpt)/denom;   % Wikipedia's Q21  bottom right
        wts(2,2) = (lonpt-lons(1))*(latpt-lats(1))/denom;   % Wikipedia's Q22  upper right
    end
end
