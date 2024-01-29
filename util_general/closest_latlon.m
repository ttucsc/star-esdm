function [latlon_ix, keepers, out_latlons, dist] = closest_latlon(latlon_pts, lats, lons, tolerance, make_grid)
%   Returns indices of closest location in [lats,lons] to each [latpt, lonpt]
%       lats,lons can be a set of unique lats & lons, or gridpoints (such as from model netcdf data).
%
%       distance is distance in degrees from nearest point. 
%   
%   [latlon_ix, keepers, out_latlons, dist] = closest_latlon(latlon_pts, lats, lons, tolerance, only_matches)
%
%           Note:  lats, lons & tolerance are assumed to be in degrees.
%   Inputs:
%       latlon_pts      nx2 array of lats & lons to look for in lats & lons
%       lats, lons      nx1 arrays of lats & lons to search
%                           if ~make_grid, #lats must equal #lons, & code returns closest point in [lats,lons].
%                           if  make_grid, then lats & lons represent a grid of points, and will be expanded into a set
%                                          of [lat,lon] points using meshgrid(...) 
%       tolerance       accept closest point within this tolerance (default:  1e-4 in degrees, or approx 11 m)
%       make_grid       t/f [false]. if true, lats & lons, generate a grid of all possible locations from lats & lons
%                           use make_grid=true if lats & lons are from a gridded dataset instead of a set of [lat,lon]
%                           values to search for the minimum distance.
%
%   Outputs:
%       latlon_ix       array of indices of point in lats,lons closest to each latlon_pt.
%                         NOTE:  each point  in latlon_pts will have a closest point, which may not  be within the 
%                                specified tolerance.  See keepers below.
%                           if ~make_grid, latlon_ix is n x 1 .  
%                           if  make_grid, latlon_ix is n x 2, giving [latix, lonix] of closest point. 
%       keepers         boolean array flagging which points are within tolerance.
%       out_latlons     nx2 array of lat & lon of closest points.
%       dist            distance of each [lat_pt, lon_pt] to nearest [lat,lon] location.
%
    if (~exist("tolerance","var") || isempty(tolerance)), tolerance = 1e-4; end
    if (~exist("make_grid","var") || isempty(make_grid)), make_grid = false; end

    nlatlons = size(latlon_pts, 1+make_grid);
    latlon_ix = nan(nlatlons,1);
    keepers = false(nlatlons,1);
    out_latlons = nan(nlatlons,2);
    dist = nan(nlatlons,1);
    
    if (make_grid)
        nlats_in = length(lats);
        nlons_in = length(lons);
        [lats,lons] = meshgrid(lats,lons);
        lats = lats(:);
        lons = lons(:);
    end
    
    for i=1:nlatlons
        d = distance(latlon_pts(i,1),latlon_pts(i,2), lats,lons);
        ix = find(d == min(d),1);
        keepers(i) = d(ix) <= tolerance;
        out_latlons(i,:) = [lats(ix), lons(ix)];
        dist(i) = d(ix);
        if (~make_grid)
            latlon_ix(i) = ix;
        else
            lat_ix = ceil(ix/nlons_in);
            lon_ix = mod(ix-1, nlats_in)+1;
            latlon_ix(i,:) = [lat_ix, lon_ix];
        end            
    end
end

