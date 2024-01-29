function [run_tbl, csvname] = valid_models(varargin)
%      function run_tbl =valid_models("model_set", "cmip5"|"cmip6", "varname",[var1,var2...],"models",[...],"scenario",[...],"ensemble",[...], "grgn", [...])
%                   where each kwd,value pair is optional.
%           Use a minus sign to exclude a specific value.  For example:
%
%               tbl = valid_models("ensembles","-r1i1p1");
%
%           would exclude r1i1p1 and include all other ensembles.  Note: ["all","-r1i1p1"] is identical to "-r1i1p1".
%
%       returns table with list of valid models matching selection.  
%       table columns are:  model, ensemble, varname, scenario [ & grgn for cmip6]
% 
%   NOTE:  "models" must be plural to differentiate it from "model_set"

    [model_set, varnames, models, scenarios, ensembles, grgns, excludes, valid_sets, csvname] = init_params(varargin{:});
    
    run_tbl = find_valid_sets(valid_sets, model_set, varnames, models, scenarios, ensembles, grgns, excludes);
    
    nsets = size(run_tbl,1);
    run_tbl.filename = strings(nsets,0);
    for i=1:nsets
        rt=run_tbl(i,:);
        if (strcmp(model_set,"cmip6"))
            run_tbl.filename(i) = sprintf("%s_day_%s_%s_%s_%s_%d-%d.llt.nc",rt.varname, rt.model, rt.scenario, rt.ensemble, rt.grgn, rt.start_year, rt.end_year);
        elseif (strcmp(model_set,"cmip5"))
            run_tbl.filename(i) = sprintf("%s.%s.%s.%s.day.%4d.%4d.llt.nc", rt.model, rt.ensemble, rt.varname, rt.scenario, rt.start_year, rt.end_year);
        end
    end
            
    
end

function  run_tbl = find_valid_sets(run_tbl, model_set, varnames, models, scenarios, ensembles, grgns, excludes)
    
    obsvnames=["Tmax","Tmin","Prcp"];
    vnames=["tasmax","tasmin","precip"];
    
        % in case user specified observation variable names, translate them to equivalent names in the models.
    for i=1:length(varnames)
        ix = find(strcmpi(varnames(i), obsvnames),1);
        if (~isempty(ix))
            varnames(i) = vnames(ix);
        end
    end

    if (~isempty(varnames))
        keepers = ismember(lower(run_tbl.varname),lower(varnames));
        run_tbl = run_tbl(keepers,:);
    end
    
    if (~isempty(models))
        keepers = ismember(lower(run_tbl.model),lower(models));
        run_tbl = run_tbl(keepers,:);
    end

    if (~isempty(scenarios))
        keepers = ismember(lower(run_tbl.scenario),lower(scenarios));
        run_tbl = run_tbl(keepers,:);
    end

    if (~isempty(ensembles))
        keepers = ismember(lower(run_tbl.ensemble),lower(ensembles));
        run_tbl = run_tbl(keepers,:);
    end
    
    if (strcmp(model_set,"cmip6") && ~isempty(grgns))
        keepers = ismember(lower(run_tbl.grgn),lower(grgns));
        run_tbl = run_tbl(keepers,:);
    end
    
    run_tbl = remove_excludes(run_tbl, excludes);    
end

function [model_set, varnames, models, scenarios, ensembles, grgns, excludes, valid_sets, csvname] = init_params(varargin)

    p = inputParser;
    p.StructExpand = true;
    p.KeepUnmatched = true;
    
    addParameter(p,'varnames',["tasmax","tasmin"]);
    addParameter(p,'models',strings(0,0),@(s) ischars(s));
    addParameter(p,'scenarios',strings(0,0),@(s) ischars(s));
    addParameter(p,'ensembles',strings(0,0),@(s) ischars(s));
    addParameter(p,'grgns', strings(0,0),@(s) ischars(s));
    addParameter(p,"model_set", strings(0,0),@(s) ischars(s));
    
    parse(p, varargin{:});
    varnames  = string(p.Results.varnames);
    models    = string(p.Results.models);
    scenarios = string(p.Results.scenarios);
    ensembles = string(p.Results.ensembles);
    grgns     = string(p.Results.grgns);
    model_set = string(p.Results.model_set);
    
    excludes = cell(1,5);
    [varnames,  excludes{1}] = check_for_excludes(varnames);
    [models,    excludes{2}] = check_for_excludes(models);
    [scenarios, excludes{3}] = check_for_excludes(scenarios);
    [ensembles, excludes{4}] = check_for_excludes(ensembles);    
    [grgns,     excludes{5}] = check_for_excludes(grgns);    
        
    if (~isempty(fieldnames(p.Unmatched))) 
        fprintf("error: unexpected input parameters: ");
        disp(p.Unmatched);
        error("unexpected input parameters");
    end
    
    obsvnames = ["Tmax",  "Tmin",  "Prec","Pr","Prcp","precip","precipitation"];
    vnames    = ["tasmax","tasmin","pr",  "pr","pr",  "pr"];
    if (isempty_s(varnames) || (length(varnames)==1 && strcmpi(varnames,'all')))
        varnames = ["tasmax","tasmin"];
    else
        for i=1:length(varnames)
            ix = find(strcmpi(varnames(i),obsvnames),1);
            if (~isempty(ix)), varnames(i) = vnames(ix); end
            if (~ismember(lower(varnames(i)), vnames)), error("unknown varname: %s", varnames(i)); end 
        end
    end
    
    if (isempty(model_set) || strlength(model_set) == 0)
        if (~isempty(grgns))
            model_set = find_model_set(grgns);
        elseif (~isempty(models))
            model_set = find_model_set(models);
        elseif (~isempty(scenarios))
            model_set = find_model_set(scenarios);
        elseif (~isempty(ensembles))
            model_set = find_model_set(ensembles);
        else
            error("cannot determine model set");
        end
    elseif (~any(strcmp(model_set,["cmip5","cmip6"])))
        error("bad model set:  %s.  Must be ""cmip5"" or ""cmip6"" "); 
    end
    
    if (strcmp(model_set,"cmip6"))
        csvname = "valid_models_cmip6.csv";
        valid_sets = readtable(csvname);
            % make text column strings...
        valid_sets.varname  = string(valid_sets.varname);
        valid_sets.model    = string(valid_sets.model);
        valid_sets.scenario = string(valid_sets.scenario);
        valid_sets.ensemble = string(valid_sets.ensemble);
        valid_sets.grgn     = string(valid_sets.grgn);
        valid_sets.calendar = string(valid_sets.calendar);
      
        mdls  = unique(valid_sets.model);
        scens = unique(valid_sets.scenario);
        ens   = unique(valid_sets.ensemble);
        gggs  = unique(valid_sets.grgn);        
    else
        csvname = "valid_models_cmip5.csv";
        valid_sets = readtable(csvname);
        mdls  = unique(valid_sets.model);
        scens = unique(valid_sets.scenario);
        ens   = unique(valid_sets.ensemble);
            % add empty grgn column to simplify 
        nsets = size(valid_sets,1);
        valid_sets.grgn = strings(nsets,1);
    end
    
    if (isempty_s(models) || (length(models)==1 && strcmpi(models,"all")))
        models = mdls;
    else
        for i=1:length(models)
            if (~any(strcmpi(models(i),mdls))), error("unknown model: %s", models(i)); end
        end
    end    
    
    if (isempty_s(scenarios) || (length(scenarios)==1 && strcmpi(scenarios,'all')))
        scenarios = scens;
    else
        for i=1:length(scenarios)
            if (~any(strcmp(scenarios(i),scens))), error("unknown scenario: %s", scenarios(i)); end
        end
    end
    
    if (isempty_s(ensembles) || (length(ensembles)==1 && strcmpi(ensembles,'all')))
        ensembles=ens;
    else
        for i=1:length(ensembles)
            if (~any(strcmp(ensembles(i),ens))), error("unknown ensemble: %s", ensembles(i)); end
        end
    end   
    
    if (strcmp(model_set, "cmip6"))
        if (isempty_s(grgns) || (length(grgns)==1 && strcmpi(grgns,'all')))
            grgns=gggs;
        else
            for i=1:length(grgns)
                if (~any(strcmp(grgns(i),gggs))), error("unknown grgn: %s", grgns(i)); end
            end
        end  
    end
    
end

function [mylist, myexcludes] = check_for_excludes(mylist)
%   removes name from list if preceded by negative sign.
%
    keepers = true(size(mylist));
    myexcludes = strings(0,0);
    if (isempty_s(mylist)), return; end
    for i=1:length(mylist)
        if (strcmp(extractBefore(mylist(i),2),"-"))
            keepers(i) = false;
            mylist(i) = extractAfter(mylist(i),1);
            myexcludes(end+1) = mylist(i); %#ok<AGROW>
        end
    end
    mylist = mylist(keepers);
end

function  run_tbl = remove_excludes(run_tbl, excludes)

    cols = ["varname","model","scenario","ensemble"];   % order of excludes doesn't match order of columns in table...

    keepers = true(size(run_tbl,1),1);
    for j=1:4
        colname = cols(j);
        for i=1:length(excludes{j})
            ex = excludes{j}(i);
            keeps = ~strcmp(run_tbl.(colname), ex);
            keepers = keepers & keeps;
        end
    end
    run_tbl = run_tbl(keepers,:);
end


% function run_tbl = sort_run_tbl(run_tbl)
% 
%     run_tbl.mdl = upper(run_tbl.model);
%     run_tbl.ens = upper(run_tbl.ensemble);
%     run_tbl.vnm = upper(run_tbl.varname);
%     run_tbl.scn = upper(run_tbl.scenario);
%     
%     run_tbl.ens = strrep(run_tbl.ens,"R10","R9z");
%     
%     run_tbl = sortrows(run_tbl,["scn","vnm","mdl","ens"]);
%     run_tbl.mdl=[];
%     run_tbl.ens=[];
%     run_tbl.vnm=[];
%     run_tbl.scn=[];
% end
% 
% function [base_dir, hostname, on] = data_base_dir(typ,varname,scenario, model_set)
% 
%     [hostname,on] = get_hostname();
%     if (strcmp(model_set,"cmip5"))
%         if (strcmp(typ,"obs"))    % Obs directories need fixing, Ian!
%             error("oops.  Obs folder stuff needs fixing!");
% %             if (on.neys)
% %                 base_dir="/Volumes/lacie_1/data/gcm_cmip5/rotated_netcdf";
% %             elseif (on.laptop)
% %                 base_dir="/Volumes/2018_1/data";
% %             elseif (on.dev_system)
% %                     base_dir="/Volumes/jcsf_data/data/gcm_cmip5/daily.1900.2100/rotated"; 
% %             elseif (on.KMQ)
% %                     base_dir="/data/obs/stations_netcdf";
% %             elseif (on.hpcc_system)
% %                     error("don't know the right folder on HPCC (hosrname: %s\n",hostname); 
% %             else
% %                 error("error:  unknown system:  %s", hostname);
% %             end
%         elseif (strcmp(typ,"model"))    
%             if (on.neys)
%                 base_dir="/Volumes/lacie_1/data/gcm_cmip5/rotated_netcdf";
%             elseif (on.laptop)
%                 base_dir="/Volumes/2018_1/data";
%             elseif (on.dev_system)
%                 base_dir="/Volumes/jcsf_data/data/gcm_cmip5/daily.1900.2100/rotated"; 
%             elseif (on.KMQ)
%                 base_dir=sprintf("/data/gcm_cmip5/daily.1900.2100/%s/%s/", varname, scenario);
%             elseif (on.hpcc_system)
%                 base_dir = "/lustre/scratch/iscottfl/cmip5_rotated";
%             else
%                 error("error:  unknown system:  %s", hostname);
%             end
%         else
%             error("error:  unknown typ: %s",typ);
%         end
%     else
%         if (strcmp(typ,"obs"))    % Obs directories need fixing, Ian!
%             error("oops.  Obs folder stuff needs fixing!");
% %             if (on.neys)
% %                 base_dir="/Volumes/lacie_1/data/gcm_cmip5/rotated_netcdf";
% %             elseif (on.laptop)
% %                 base_dir="/Volumes/2018_1/data";
% %             elseif (on.dev_system)
% %                     base_dir="/Volumes/jcsf_data/data/gcm_cmip5/daily.1900.2100/rotated"; 
% %             elseif (on.KMQ)
% %                     base_dir="/data/obs/stations_netcdf";
% %             elseif (on.hpcc_system)
% %                     error("don't know the right folder on HPCC (hosrname: %s\n",hostname); 
% %             else
% %                 error("error:  unknown system:  %s", hostname);
% %             end
%         elseif (strcmp(typ,"model"))    
%             if (strcmp(scenario,"hist")), scenario = "historical"; end
%             if (on.neys)
%                 base_dir=sprintf("/Volumes/lacie_1/data/cmip6_rotated/%s", scenario);
%             elseif (on.laptop)
%                 base_dir="/Volumes/2018_1/data";
%             elseif (on.dev_system)
%                 error("oops.  cmip6 files not on icsf-jmac");
% %               base_dir="/Volumes/jcsf_data/data/gcm_cmip5/daily.1900.2100/rotated"; 
%             elseif (on.KMQ)
%                 error("oops.  cmip6 files not on KMQ");
% %               base_dir=sprintf("/data/gcm_cmip5/daily.1900.2100/%s/%s/", varname, scenario);
%             elseif (on.hpcc_system)
%                 base_dir = sprintf("/lustre/research/hayhoe/cmip6_rotated/%s", scenario);
%             else
%                 error("error:  unknown system:  %s", hostname);
%             end
%         else
%             error("error:  unknown typ: %s",typ);
%         end
%     end
%         
% end

