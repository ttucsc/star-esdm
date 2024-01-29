function models = all_ensembles_cmip6
%       this is a quick-and-dirty.  This table should have 5 columns:
%   model, ix, calendar, nlats, nlons

    tbl = readtable("valid_models_cmip6.csv");
    models = table(unique(tbl.ensemble), 'VariableNames',["ensemble"]);
end

