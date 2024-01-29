function clim_vnames = climate_varnames()
%   list of expected variable names and their variations for TTU CSC climate data.
%   note that you may want to test with
%       any(strcmpi(myvname, std_varnames())
%   to avoid case sensitivity.
%       12/4/21 icsf:  added prcp.  (earlier added precipitation and packed_zvals)

    clim_vnames   = [  "tmax","tasmax",  "tmin","tasmin","tavg","temp_F","temp_C","Prec","prcp","precipitation",  "pr","rh","rh_F","rhsmax","rhsmin","relhum","packed_zvals","zvals"];

end