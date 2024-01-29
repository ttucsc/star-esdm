function std_clim = std_clim_name(vname)
% std_clim = std_clim_name(vname)
%   returns a standard name for vname to simplify testing elsewhere in the code.
%       std_clim    "tasmin", "tasmax", "prcp" etc. if it is a variation on one of the climate variables.
%
%   See also functions is_climate_variable(vname), clim_category(varname), is_temp_variable(varname), is_prcp_variable(varname) 

    [~, std_clim] = is_climate_variable(vname);
end