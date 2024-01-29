function clim_cat = clim_category(vname)
% clim_cat = clim_category(vname)
%   returns general category of climate variable vname
%       category is one of:
%           temp        temperature             (tmax, tmin, tasmax, tasmin, etc.)
%           prcp        precipitation           (pr, prcp, precip, etc.)
%           rh          relative humidity       (rh, rhs, relhum, etc.)
%           zvals       probability z-scores    (zvals, packed_zvals)
%
    [~,~,clim_cat] = is_climate_variable(vname); 
end

