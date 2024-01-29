function model_set = find_model_set(mdl_scen_ens_grgn, cmip5_tbl, cmip6_tbl)
%
%
%   returns "cmip5" or "cmip6" if the models/scenarios/ensembles/grgns define a valid set in
%   either cmip5 or cmip6 data.
%   This relies on the presence of the files valid_models_cmip5.csv and valid_models_cmip6.csv in source code folder.
%       these files can be created with the shell script "generate_valid_model_csv_files.sh
%   Also assumes that the models/scenarios/ensembles/grgns are mutually
%   exclusive between cmip5 and cmip6.  (Currently that is the case.)
%   
%
%   Inputs:
%       mdl_scen_ens_grgn   either a single value or an array of strings
%                               values can be either a model, scenario, ensemble or grgn (gridding key)
%       cmip5_tbl           (optional) read in table of all valid cmip5 combos of variable/model/ensemble/scenario/grgn's
%       cmip6_tbl           (optional) read in table of all valid cmip6 combos of variable/model/ensemble/scenario/grgn's
%                               Code will read  these tables from .csv files if not provided
%
%   returns:  

    mdl_scen_ens_grgn = string(mdl_scen_ens_grgn);
    model_set = strings(0,0);
    
    if (~exist("cmip5_tbl", "var") || isempty(cmip5_tbl))
        cmip5_tbl=readtable("valid_models_cmip5.csv");
    end
    if (~exist("cmip6_tbl", "var") || isempty(cmip6_tbl))
        cmip6_tbl=readtable("valid_models_cmip6.csv");
    end
    
    if (length(mdl_scen_ens_grgn) > 1)
        for i=1:length(mdl_scen_ens_grgn)
            if (strlength(mdl_scen_ens_grgn(i))>0)
                model_set = find_model_set(mdl_scen_ens_grgn(i), cmip5_tbl, cmip6_tbl);
            end
        end
        if (strlength(model_set)>0), return; end
%         if (~isempty(model_set))
%             model_set = unique(model_set);
%             return; 
%         end
        error("error:  unknown cmip5/cmip6 values:  %s\n", join(mdl_scen_ens_grgn," "));
    else
    
                % see if it matches any of the models
        m5=unique(cmip5_tbl.model);
        m6=unique(cmip6_tbl.model);

        if (any(strcmp(mdl_scen_ens_grgn, m5)))
            model_set = "cmip5";
            return
        end
        if (any(strcmp(mdl_scen_ens_grgn, m6)))
            model_set(end+1) = "cmip6";
            return
        end

                % see if it matches any of the scenarios
        s5 = unique(cmip5_tbl.scenario);
        s6 = unique(cmip6_tbl.scenario);

        if (any(strcmp(mdl_scen_ens_grgn, s5)))
            model_set = "cmip5";
            return
        end
        if (any(strcmp(mdl_scen_ens_grgn, s6)))
            model_set(end+1) = "cmip6";
            return
        end
    
                % see if it matches any of the ensembles
        e5 = unique(cmip5_tbl.ensemble);
        e6 = unique(cmip6_tbl.ensemble);

        if (any(strcmp(mdl_scen_ens_grgn, e5)))
            model_set = "cmip5";
            return
        end
        if (any(strcmp(mdl_scen_ens_grgn, e6)))
            model_set(end+1) = "cmip6";
            return
        end
        
                % see if it matches any of the cmip6 grgn's
        g6 = unique(cmip6_tbl.grgn);

        if (any(strcmp(mdl_scen_ens_grgn, g6)))
            model_set(end+1) = "cmip6";
            return
        end
        
        error("error:  unknown cmip5/cmip6 value:  %s\n", mdl_scen_ens_grgn);
        
    end
    
end
    
