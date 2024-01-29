function models = all_models_cmip6
%       this is a quick-and-dirty.  This table should have 5 columns:
%   model, ix, calendar, nlats, nlons

    tbl = readtable("valid_models_cmip6.csv");
    models = table(unique(tbl.model), 'VariableNames',["model"]);
end

