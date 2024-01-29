function [outmask, outlats, outlons]  = is_land(lats, lons, threshold, no_antarctica, ls_mask, mask_lats, mask_lons)
%   returns boolean land/sea status of each point in lats & lons.  
%       lats & lons can specify a lat/lon rectangle, or can specify individual lat/lon locations.
%
%   Interpolates from either the ERA5 .1 degree land/sea mask or a provided landsea mask,
%       if threshold is not empty,  thresholds to decide if each point is a land or sea (water) point.
%       If threshold is empty, then outmask is values between 0 & 1, representing fraction of land for each lat/lon.
%           (0 = no land;  1 = all land.)  Otherwise, mask is thresholded and outmask is logical:   true=land; false=water
%       default threshold, if not present, is 0.5.

%   If size(lats) ~= size(lons) or [size(lats) == size(lons) and lats is a column and lons is a row,
%       then  each (lat,lon) pair is treated a location to check.  
%       Otherwise, it generates a meshgrid from lats & lons, and returns a boolean (logical) mask for lats,lons
%           If lats & lons are same size but a meshgrid is desired, make lats a column vector and lons a row vector.
%
%   If no_antarctica is true [def:  false], then points south of -60 lat are returned as ocean.
%
%   If ls_mask is missing or empty, uses ECMWF's .1-degree land/sea mask thresholded at threshold (threshold default:  0.5)
%   Otherwise, uses ls_mask as a binary mask, rounds each lat/lon to nearest mask_lat,mask_lon
%

    if (~exist("no_antarctica","var") || isempty(no_antarctica)),   no_antarctica = false; end
    if (~exist("threshold","var")), threshold = 0.5;  end

    if (~exist("lsmask","var"))
        [outmask, outlats,outlons] = landsea_mask(lats,lons,threshold);
    else
        [outmask, outlats,outlons] = landsea_mask(lats,lons,threshold,ls_mask,mask_lats,mask_lons);
    end
    
    if (no_antarctica)
        south_pts = outlats < -60;
        outmask(south_pts) = false;     % works whether outmask is logical or double.
    end
        
end
        
        
    
        

    