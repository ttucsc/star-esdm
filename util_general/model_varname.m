function varname = model_varname(vname)
% returns the variable name used in TTU CSC model netcdf files for the equivalent observation netcdf variable name

    if (contains(lower(vname),["tasmax","tmax","maxtemp"]))
        varname = "tasmax";
    elseif (contains(lower(vname),["tasmin","tmin","mintemp"]))
        varname = "tasmin";
    elseif (contains(lower(vname),["pr","precip","prcp","precipitation"]))
        varname = "pr";
    elseif (contains(lower(vname),["rh","relhum","hum"]))
        varname = "rh";
    else
        error("unknown_climate_variable: %s", vname);
    end
end