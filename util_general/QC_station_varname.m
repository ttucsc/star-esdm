function vname = obs_varname(varname)
% returns the observation variable name used in TTU CSC observation netcdf files for the equivalent model netcdf
% variable name (or any obvious equivalent)
    if (isempty(varname))
        vname = varname;
    elseif (contains(lower(varname),["tasmax","tmax","maxtemp"]))
        vname = "Tmax";
    elseif (contains(lower(varname),["tasmin","tmin","mintemp"]))
        vname = "Tmin";
    elseif (contains(lower(varname),["pr","precip","prcp","precipitation"]))
        vname = "Prcp";
    elseif (contains(lower(varname),["rh","relhum","hum"]))
        vname = "rh";
    else
        error("unknown_climate_variable: %s", varname);
    end
end