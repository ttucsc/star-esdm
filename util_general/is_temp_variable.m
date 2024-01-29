function tf = is_temp_variable(vname)
% tf = is_temp_variable(vname)
%   returns true if vname is a temperature variable name (tmax, tmin, tasmax, tasmin, etc.), false otherwise
%
    [~,~,clim_cat] = is_climate_variable(vname); 
    tf = strcmp(clim_cat,"temp");
end

