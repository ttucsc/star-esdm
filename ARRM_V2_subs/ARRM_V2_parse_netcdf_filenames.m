function run_info = ARRM_V2_parse_netcdf_filenames(fnames,terms, run_info, no_stations)
%
%   Parses one or more netcdf filenames for model, variable, scenario and ensemble
%   Returns a run_info struct with fields model, variable, scenario and ensemble
%       fields are set to the first matching value found in the filenames.
%
%   If an array of terms is provided, you can limit the terms it searches for.  Otherwise searches for all four terms.
%   If a run_info struct is passed in, then only the missing fields are filled.
%
%   fnames and terms can be either cell arrays of chars or string arrays.
%   fnames should either be all cmip5 or all cmip6 files.  The first filename is checked to see if it is cmip5 or cmip6.
%
%   Inputs:
%
%       fnames              array of filenames (strings or cell array of chars)
%       terms               (optional)  array of fields desired from filename (model, varname, etc.)
%       run_info            (optional)  run_info struct or class object with fields (or properties) to store results in.
%                                           NOTE:  if parsing multiple filenames, runinfo must be empty or missing.
%                                           NOTE:  if run_info is a class, then it must have a property for each term
%       no_stations         (optional)  true/false.  If true, does not consider "stations" as a valid model type.
%
%   Outputs:
%       run_info            a struct with fields (or properties) for each term requested
%                               if terms empty or not present, will parse for varname, model, scenario, ensemble and grgn
%                               if multiple filenames, each field of run_info will be an array.
%


    if (~exist('terms','var') || isempty(terms))
        terms= ["varname","model","scenario","ensemble","grgn"];
    end
    if (~exist('run_info','var') || isempty(run_info))
        run_info = struct;
    end
    if (~exist('no_stations','var') || isempty(no_stations))
        no_stations = false; 
    end
    
    
    if (isempty_s(fnames) && ~exist("run_info","var")), run_info = struct; return; end
    
    if (ischar(fnames) && ~iscell(fnames)), fnames={fnames}; end  
    
    if (any(contains(fnames,"_day_")) && any(contains(fnames,".day.")))
        error("error:  can't parse both cmip5 and cmip6 filenames\n");    
    elseif (contains(fnames{1},"_day_"))
        model_set = "cmip6";
    elseif (contains(fnames{1},".day."))
        model_set = "cmip5";
    else
        error("error: can't determine model set");
    end
    
    
        % now read the valid_model table (or tables)
%   for i=1:length(model_set)
        i=1;    % quick-fix, since we only handle cmip5 or cmip6 now.
        tblname = sprintf("valid_models_%s.csv", model_set(i));
%       if (any(strcmp(model_set,"cmip6")))
            try
                t = readtable(tblname);
            catch
                t = valid_models("model_set",model_set(i));
            end
%            if (~exist("tbl","var"))
                tbl = t;
%            else
%                tbl = [tbl;t]; %#ok<AGROW>
%            end
%       end
%    end
    
    parts.varname  = unique(tbl.varname);
    parts.model    = unique(tbl.model);
    parts.scenario = unique(tbl.scenario);
    parts.ensemble = unique(tbl.ensemble);
    try
        parts.grgn = unique(tbl.grgn);
    catch
        parts.grgn = strings(0);
    end
%     parts.varname = ["pr", "precip", "prcp", "rhs", "rhsmax", "rhsmin", "rsds", "tas", "tasmax", "tasmin", "tos","tmax","tmin",];
%     parts.model   = ["ACCESS1-0", "ACCESS1-3", "BNU-ESM", "CCSM4", "CESM1-BGC", "CESM1-FASTCHEM", "CMCC-CESM", "CMCC-CM", ...
%                       "CMCC-CMS", "CNRM-CM5", "CSIRO-Mk3-6-0", "CanCM4", "CanESM2", "EC-EARTH", "GFDL-CM3", "GFDL-ESM2G",  ...
%                       "GFDL-ESM2M", "HadCM3", "HadGEM2-CC", "HadGEM2-ES", "IPSL-CM5A-LR", "IPSL-CM5A-MR", "IPSL-CM5B-LR",  ...
%                       "MIROC-ESM", "MIROC-ESM-CHEM", "MIROC4h", "MIROC5", "MPI-ESM-LR", "MPI-ESM-MR", "MPI-ESM-P",  ...
%                       "MRI-CGCM3", "NorESM1-M", "bcc-csm1-1", "bcc-csm1-1-m", "inmcm4","stations"];
%     parts.scenario = ["hist","esmrcp85", "rcp26", "rcp45", "rcp60", "rcp85", "sst2090rcp45"];
%     parts.ensemble = ["r10i1p1", "r11i1p1", "r12i1p1", "r13i1p1", "r14i1p1", "r1i1p1", "r1i1p2", "r2i1p1", "r2i1p2", ...
%                       "r3i1p1", "r3i1p2", "r4i1p1", "r4i1p2", "r5i1p1", "r5i1p2", "r6i1p1", "r6i1p2", "r6i1p3", "r7i1p1", ...
%                       "r8i1p1", "r9i1p1"];
     
    if (no_stations)
        parts.model = parts.model(~strcmp("stations",parts.model));
    end
    for t=1:length(terms)
        term=terms{t};
        for f=1:length(fnames)
            fname=fnames{f};
%           if (~any(strcmp(model_set,"cmip6")) && strcmp(term,"grgn")), continue; end
            [~,filebase,~] = fileparts(fname);
            if (strcmp(model_set, "cmip6"))
                fparts = split(filebase,'_');
            else
                fparts = split(filebase,".");
            end
            for i=1:length(fparts)
%               if (~isempty(run_info.(term))), continue; end
                fpart=fparts{i};
                ix=find(strcmpi(parts.(term), fpart), 1);
                if (~isempty(ix))
                    if (isfield(run_info,term) || isprop(run_info, term))
                        if (isempty(run_info.(term)))
                            run_info.(term) = strings(0,0);
                        end
                        run_info.(term)(f) = string(fpart);
                    else
                        run_info.(term)    = string(fpart);
                    end
                end
            end
        end
    end
end