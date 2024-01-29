function [yn, std_clim, clim_cat] = is_climate_variable(vnms, skip_zvals)
%   returns true if vnms matches (case-insensitive) to any of the standard climate variable names we encounter...
%
%   Also returns a 'standard name" std_clim and climate category for simpler testing elsewhere in our code.
%       std_clim    returns tasmin, tasmax, prcp  etc. if it is a variation on one of the climate variables.
%       clim_cat    "temp" for temperature variable, "prcp" for precip variable, "rh" for humidity variable, "zvals" for & packed_zvals
%                       (helpful to simplyl identify the category for setting default parameters
%
%   if yn is false, std_clim and clim_cat are zero-length strings ""
%   if vnms is an array of variable names, yn, std_clim and clim_cat values are returned as arrays.
%
%   See also functions std_clim_name(varname) and clim_category(varname), is_temp_variable(varname), is_prcp_variable(varname), is_rh_variable(varname), is_zvals_variable(vname) 
%
%     clim_vnames   = [  "tmax","tasmax",  "tmin","tasmin","tavg","temp_F","temp_C","Prec","prcp","precipitation",  "pr","precip","rh","rh_F","rhsmax","rhsmin","relhum","packed_zvals","zvals"];
%     %                      1        2        3        4      5        6        7      8      9              10     11       12   13     14       15       16        17            18      19  
%     std_clim_name = ["tasmax","tasmax","tasmin","tasmin","tavg","temp_f","temp_c","prcp","prcp",         "prcp","prcp",  "prcp","rh",  "rh","rhsmax","rhsmin",    "rh",       "zvals","zvals"];
%     clim_category = [  "temp",  "temp",  "temp",  "temp","temp", "temp",   "temp","prcp","prcp",         "prcp","prcp",  "prcp","rh",  "rh",    "rh",    "rh",    "rh"        "zvals","zvals"];
%
%   7/28/2022 icsf added skip_zvals option.

    if (~exist("skip_zvals","var") || isempty(skip_zvals)), skip_zvals = false; end
   
    clim_vnames   = [  "tmax","tasmax",  "tmin","tasmin","tavg","temp_F","temp_C","Prec","prcp","precipitation",  "pr","precip","rh","rh_F","rhsmax","rhsmin","relhum"];
    %                      1        2        3        4      5        6        7      8      9              10     11       12   13     14       15       16        17            18      19  
    std_clim_name = ["tasmax","tasmax","tasmin","tasmin","tavg","temp_f","temp_c","prcp","prcp",         "prcp","prcp",  "prcp","rh",  "rh","rhsmax","rhsmin",    "rh"];
    clim_category = [  "temp",  "temp",  "temp",  "temp","temp", "temp",   "temp","prcp","prcp",         "prcp","prcp",  "prcp","rh",  "rh",    "rh",    "rh",    "rh"];
   
    vnms = string(vnms);
    yn       = false(size(vnms));
    std_clim = strings(size(vnms));
    clim_cat = strings(size(vnms));
    
    for i=1:length(vnms)
        vnm = vnms(i);
        if (skip_zvals && ~isempty(regexp(vnm,"zvals$","ONCE")))
            yn(i) = false;
        else
            ix = find(strcmpi(vnm, clim_vnames),1,"first");
            if (isempty(ix))
                yn(i) = ~isempty(regexp(vnm,"zvals$","ONCE"));
                if (yn(i))
                    std_clim(i) = "zvals";
                    clim_cat(i) = "zvals";
                end
            else
                yn(i) = true;
                std_clim(i) = std_clim_name(ix);
                clim_cat(i) = clim_category(ix);
            end
        end
    end
        
end