function [lsmask, outlats, outlons] = landsea_mask(lats,lons, threshold, bigmask_name_or_bigmask, blats_or_varname,blons)
% lsmask = landsea_mask(lats,lons, threshold, bigmask_name)
%   
%   Returns land/sea mask at resolution given by lats,lons, interpolated
%       from ECMWF's ERA5 .1 degree land-sea mask available at 
%       https://confluence.ecmwf.int/pages/viewpage.action?pageId=140385202#ERA5Land:datadocumentation-parameterlistingParameterlistings
%
%       lats,lons s/b lats & lons of model  (e.g. from something linke ncdump -v lat mdlname...)
%       If lats & lons are column & row vectors (respectively) or different lengths, then a meshgrid is generated from
%       them.  Otherwise they are treated as a set of lat,lon points to be interpolated to.
%   
%       threshold   empty/missing or value between 0 & 1.
%                       if (empty/missing, mask is fraction of land within the .1x.1 degree gridcell centered @ lat,lon
%                       if present, is used to threshold map to binary (logical) mask table
%
%       bigmask_name_ot_bigmask    name of netcdf mask file, or mask matrix.  If emtpy or missing, uses ECMWF's ERA5 land/sea mask.
%                                       if mask matrix, must also provide lats & lons -- blats & blons -- for mask matrix.
%       blats_or_varname           if bigmask is alpha, then blats must be varname of variable to load for mask
%
%   Returns:
%       lsmask      either (approx) fraction of gridbox (area inside latbnds,lonbnds) which is land OR
%                          logical array, true where lat/lon location is land, false where it is water
%                               if threshold is present and not empty, fractional values are thresholded to arrive at tru/false status.
%       outlats     lat/lon pairs for each point (regardless of whether point is land or water) 
%                               depending on size & orientation of input lats & lons, can either be the input lats
%                               and lons, or can be a meshgrid of the input lats & lons.
%
%   NOTE:  Code circularly interpolates between the last lon & first longitude, so locations between 359 & 0 or 179 and -180 
%          are valid.   If bigmask does not cover the full globe, then user must make sure lats and lons fall completely
%          inside the bounds of blats and blons.
%
%
    if (~exist("lats","var")), lats = []; lons = []; end % will simply return bigmask's mask &  lat/lon range. 
    if (~exist("threshold","var")), threshold = []; end
    if (~exist("bigmask_name_or_bigmask","var") || isempty(bigmask_name_or_bigmask))
        bigmask_name_or_bigmask = "lsm_1279l4_0.1x0.1.grb_v4_unpack.nc";
        varname = "lsm";
    end
    
    if (ischar(bigmask_name_or_bigmask) || isstring(bigmask_name_or_bigmask))
        if (exist("blats_or_varname","var"))
            varname = blats_or_varname;
        end
        nc = ncdf(bigmask_name_or_bigmask);
        nc.loadvars();
        [latname,lonname] = ncdf_get_llt_dimnames(nc);
        [nclats,nclons] = ncdf_get_latlons(nc);
        if (isrow(nclats))
            nclats = nclats';
        end
        if (nclats(end)<nclats(1))
            nclats = flipud(nclats);
            flip_matrix = true;
        else
            flip_matrix = false;
        end
            
        [blons,blats] = meshgrid(nclons,nclats);
        varix = find(strcmp(varname,nc.varlist));
        latix = find(strcmp(latname, nc.Variables(varix).dimlist));
        lonix = find(strcmp(lonname, nc.Variables(varix).dimlist));
        if (lonix < latix)
            bmask=nc.Variables(varix).vdata';
        else
            bmask=nc.Variables(varix).vdata;
        end
        if (flip_matrix)
            bmask = flipud(bmask);
        end
    else
        bmask = bigmask_name_or_bigmask;
            % user must also provide blats and blons!
        blats = blats_or_varname;
        if ((iscolumn(blats) && isrow(blons)) || (numel(blats) ~= numel(blons)))
            [blons,blats] = meshgrid(blons,blats);
        end
    end
    
    if (~isempty(lats))
        nlats = numel(lats);
        nlons = numel(lons);
        if ((nlats ~= nlons) || (iscolumn(lats) && isrow(lons)))
            [outlons,outlats] = meshgrid(lons,lats);
        else
            outlons = lons;
            outlats = lats;
        end
        nblats = numel(blats);
        nblons = numel(blons);
        if ((nblats ~= nblons) || (iscolumn(blats) && isrow(blons)))
            [blons,blats] = meshgrid(nclons,nclats);
        end
        
                % duplicate mask on both ends so we don't have problems with -180 to 180 vs 0-360 longitudes
                
        bmask = [bmask,bmask,bmask];
        blats = [blats,blats,blats];
        blons = [blons-360,blons,blons+360];
            
        lsmask = interp2(blons,blats,bmask, outlons,outlats);
    else
        lsmask = bmask;
        outlats = nclats;
        outlons = nclons;
    end
    
    if (~isempty(threshold))
        lsmask(lsmask>=threshold) = 1;
        lsmask(lsmask<threshold) = 0;
        lsmask = logical(lsmask);
    end
    
end
    
    