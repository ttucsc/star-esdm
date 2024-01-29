function clim_varname = find_climate_varname(vnames, maxnames)
% clim_varname = find_climate_varname(vnames, maxnames)
%
%       Returns first match in vnames to any of the standard climate varnames as given in is_climate_varname(...):
%   std_varnames = ["pr", "prc","precip", "prcp", "rhs", "rhsmax", "rhsmin", "rsds", "tas", "tasmax", "tasmin", "tos","tmax","tmin"];
%       returns an empty string array (strings(0)) if no match.
%
%   If there are multiple climate variables in the file, setting maxnames to a number or to "all"
%       will return a string array of the names in varname
%
%   purpose:  to extract the climate variable's name from the list of variables in a netcdf.
%             example:
%               nc =ncdf("daily_max_temp_file.nc");
%               vname = find_climate_varname(nc.varnames());
%                   this will return "tasmax" from the list of variables in the file
%                   assuming, of course, that there is a variable named tasmax in the file.
%               
%
    if (~isstring(vnames)), vnames = string(vnames); end
    if (~exist('maxnames','var') || isempty(maxnames))
        maxnames = 1;
    elseif (strcmpi(maxnames,"all"))
        maxnames = length(vnames);
    end
    clim_varname = strings(0,0);
    for i=1:length(vnames)
        if (is_climate_variable(vnames(i)))
            clim_varname = cat(2,clim_varname, vnames(i));
            if (length(clim_varname) >= maxnames), return; end
        end
    end
end

