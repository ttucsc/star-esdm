function tf = is_zvals_variable(vname)
% tf = is_zvals_variable(vname)
%   returns true if vname is a zvals variable name (zvals, packed_zvals), false otherwise
%
%   zvals variables contain the z-value-equivalent (# of std deviations from mean in a normal distribution)
%   of the probability of the corresponding datapoint  occurring.
%   ARRM_V2 (STAR_ESDM) code stores the zval for each downscaled daily value 
%   NOTE:  to get the cumulative probability of the zvalue  (integral of normal pdf from 0 to zval):
%       prob = normcdf(zval);
%
    [~,~,clim_cat] = is_climate_variable(vname); 
    tf = strcmp(clim_cat,"zvals");
end

