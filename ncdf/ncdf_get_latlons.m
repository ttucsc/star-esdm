function [all_lats, all_lons, file_lats, file_lons, latbnds, lonbnds] = ncdf_get_latlons(fnames)
%returns lats & lons from one or more filenames
% outputs:
%   all_lats        complete list of lats in all files
%   all_lons        complete list of lons in all files
%   file_lats     	cell array of lats in each file
%   file_lons       cell array of lons in each file
%   latbnds         cell array of latbnds in each file (if present)
%   lonbnds         cell array of lonbnds in each file (if present);
%
%       All outputs returned as doubles.  (CMIP6's IPSL-CM6A-LR has lats & lons stored as single-precision.)
%
%   modified 7/6/22 icsf to make sure returned lats,lons and ranges are doubles.  ncdf_getvar does that job for us.


    
    if (isa(fnames,"ncObj"))
        fnames={fnames}; 
    end
    if (ischar(fnames))
        fnames={fnames}; 
    end

    nfiles = numel(fnames);
    file_lats = cell(nfiles,1);
    file_lons = cell(nfiles,1);
    latbnds   = cell(nfiles,1);
    lonbnds   = cell(nfiles,1);

    for i=1:nfiles
        [file_lats{i},~,~,nc] = ncdf_getvar(fnames{i},"lat", true);     % will find latitude, Lat, LAT, etc.
        file_lons{i}          = ncdf_getvar(nc,       "lon", true);     % same for lon...
        
%         file_lons{i} = mod(file_lons{i},360);
        
        try
            latbnds{i} = double(ncdf_getvar(fnames{i}, "lat_bnds"));
        catch
        end
        try
            lonbnds{i} = double(ncdf_getvar(fnames{i}, "lon_bnds"));
        catch
        end
        
        if (i==1)
            all_lats = file_lats{i};
            all_lons = file_lons{i};
        else
            all_lats = union(all_lats, file_lats{i});
            all_lons = union(all_lons, file_lons{i});
        end
    end
    
    if (nfiles == 1)
        file_lats = file_lats{1};
        file_lons = file_lons{1};
        latbnds   = latbnds{1};
        lonbnds   = lonbnds{1}; 
    end
end

