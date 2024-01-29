function [mdl_lats, mdl_lons] = model_latlons(model, fname, dirname)

    if (~exist('dirname','var') || isempty(dirname)), dirname = "/Volumes/lacie_1/data/gcm_cmip5/rotated_netcdf"; end
    if (~exist('fname','var') || isempty(fname))
        fname = sprintf("%s.r1i1p1.tasmax.hist.day.1900.2005.llt.nc", model);
    end
    if (~strcmp(extractBefore(fname,2),"/"))
        fname = fullfile(dirname, fname);
    end
    
    if (~isfile(fname)), error("error:  file not found: %s\n", fname); end
    
    [mdl_lats,mdl_lons] = ncdf_get_latlons(fname);

end

