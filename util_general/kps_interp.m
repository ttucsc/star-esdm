function [q_out, B] = kps_interp(lat1, lat2, lon1, lon2, q11, q12, q21, q22, lats,lons)
% generalized interpolation routine.  equivalent to bilinear interpolation when the four points are in a standard
% rectangle.  Seems to have artifacts when the four points are not rectangular.
%   Interpolates for all lats & lons inside the quadrilateral defined by lat1,lat2, lon1, lon2.  This code can only work
%   on a horizontal parallelogram, but the algorithm should work with any four points, as long as they are not
%   co-linear.
%
%   See test program kps_interps(...) which calls this to interpolate over a wide area.
%
%   based on:   Kim, K-H;  Shim, P; Shin, S.  An Alternative Bilinear Interpolation Method Between Spherical Grids, 
%               Atmosphere 2019, 10(3), 123; https://doi.org/10.3390/atmos10030123
%   also at:    https://www.mdpi.com/2073-4433/10/3/123/htm
%
%       A method for interpolating a point from the four nearest neighbors when the neighbors are not on a true
%       rectangular grid.  This is useful for interpolating on a lat/lon grid.
%    dx1 = distance(lat1, lon1, lat1, lon2);     % horizontal distance in degrees
%    dx2 = distance(lat2, lon1, lat2, lon2);     % horizontal distance in degrees
    dx1 = (lon2 - lon1) * cos(lat1*pi/180);
    dx2 = (lon2 - lon1) * cos(lat2*pi/180);
    dd  = (dx1 - dx2)/2;
    
    x1 = 0;
    x2 = dx1;
    x3 = dd;
    x4 = dx1-dd;
    
    
    y1 = 0;
    y2 = 0;
%   y3 = distance(lat1, lon1, lat2, lon1);                         % distance in degrees.
    y3 = lat2 - lat1;
    y4 = y3;
    
      A = [ 1, x1, y1, x1*y1; ...
            1, x2, y2, x2*y2; ...
            1, x3, y3, x3*y3; ...
            1, x4, y4, x4*y4;
        ];
    
    C = [q11; q12; q21; q22];
    
    B = linsolve(A,C);

    if (~exist('lats','var') || isempty(lats))
        q_out = [];
    else
        nx = length(lons);
        ny = length(lats);
        q_out = nan(ny,nx);
        lat1s = repmat(lat1, ny,1);
        lon1s = repmat(lon1, 1,ny);
        lon2s = repmat(lon2, 1,ny);
%       dx = distance(lats,lon1s, lats,lon2s);
        dx = cos(lats*pi/180)*(lon2 - lon1); 
        ddx = (dx1-dx)/2;
        xfrac = (lons-lon1)/(lon2-lon1);
%       dy = distance(lats,lon1s, lat1, lon1); 
        dy = lats - lat1s';
        for i=1:ny
            for j=1:nx
                x =  xfrac(j) * dx(i) + ddx(i);
                y = dy(i);
                q_out(i,j) = sum( B .* [1, x, y, x*y]');
            end
        end
    end
end

