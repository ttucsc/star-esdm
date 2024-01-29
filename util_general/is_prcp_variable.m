function tf = is_prcp_variable(vname)
% tf = is_prcp_variable(vname)
%   returns true if vname is a precipitation variable name (prcp, pr, precip, precipitation)(, false otherwise
%
    [~,~,clim_cat] = is_climate_variable(vname); 
    tf = strcmp(clim_cat,"prcp");
end

