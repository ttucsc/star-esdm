function [fnames, base_dirs] = ARRM_V2_make_netcdf_filenames(varname, model, ensemble, scenario, data_yrs, region)     % creates netcdf filenames from model, varname, etc.
% function [fnames, base_dir] = ARRM_V2_make_netcdf_filenames(varname, model, ensemble, scenario, data_yrs, region)     % creates netcdf filenames from model, varname, etc.
%
%   returns a cell array of strings & base_dirs of the files needed to read the data covering data_yrs for varname, model, ensemble & scenario.
%   if model is "Stations", then only 1 filename and dir is in the cell.
%   else 1 or 2 filenames and dirs are returned, depending on whether all the data is in hist, or in future, or if both
%   are needed.

    fnames    = {};
    base_dirs = {};
    nncs = 0;
    [~,hostname]=system('hostname');
    hostname = strtrim(hostname);
    if (~exist("model","var"))
        model=strings();
    else
        model = string(model);
    end
    if (isempty(model) || any(strncmpi(model,"station",7)))  % create a stations filename
        if (~exist('region','var') || isempty(region)), region = "ncamerica"; end
        region = string(region);
        varname = string(varname);

        for v=1:length(varname)
            for r=1:length(region)
                nncs=nncs+1;
                vname = stations_varname(varname(v));
                fnames{nncs}(1) = sprintf("stations.%s.%s.1850.2017.nc", vname, region(r)); %#ok<AGROW>
                base_dirs{nncs}(1) = get_basedir(true, varname, hostname); %#ok<AGROW>
            end
        end
        model(strncmpi(model,"station",7)) = [];
    end
    if (~isempty(model))
        model = string(model);
        ensemble = string(ensemble);
        varname = string(varname);
        scenario = string(scenario);
        if (data_yrs(1) <= 2005 && ~contains(scenario,"hist"))
            scenario(end+1) = "hist";
        end
        for m=1:length(model)
            for e=1:length(ensemble)
                for v=1:length(varname)
                    nncs = nncs + 1;
                    mncs = 0;
                    for s=1:length(scenario)
                        if (scenario(s) == "hist")
                            yrs = [1900,2005];
                        else
                            yrs = [2006,2100];
                        end
                        mncs=mncs+1;
                        fnames{nncs}(mncs,1) = sprintf("%s.%s.%s.%s.day.%d.%d.llt.nc", model(m),ensemble(e),varname(v),scenario(s),yrs(1),yrs(2)); %#ok<AGROW>
                        base_dirs{nncs}(mncs,1) = get_basedir(false, varname, hostname); %#ok<AGROW>
                    end
                end
            end
        end
    end
    
end

function base_dir = get_basedir(is_stations, varname, hostname)

    if (~exist('hostname','var') || isempty(hostname))
        [~,hostname]=system('hostname');
        hostname = strtrim(hostname);
    end
        
    on_laptop      = contains(hostname,"jcsf") || contains(hostname,"jpro");
    on_dev_system  = contains(hostname,"icsf");
    on_hpcc_system = contains(hostname,"compute");
    
    if (on_laptop)
        base_dir  = "/Volumes/2018_1/data";
    elseif (on_dev_system)
        if (is_stations)
            base_dir = "/Volumes/jcsf_data/data/obs";
        else
            base_dir = "/Volumes/jcsf_data/data/gcm_cmip5/daily.1900.2100/rotated";
        end
    elseif (on_hpcc_system)
        base_dir = "/lustre/scratch/iscottfl/cmip5_rotated";
    else
        base_dir = sprintf("/data/gcm_cmip5/daily.1900.2100/%s", varname);
    end
end
        

function vname = stations_varname(varname)

    vnames=["Tmax","Tmin", "Prec"];
    varnames = ["tmax","tasmax", "", ""; ...
                "tmin","tasmin", "", ""; ...
                "prec","pr","prcp","precip"; ];
            
    for i=1:size(varnames,1)
        if (any(strcmpi(varname, varnames(i,:))))
            vname = vnames(i);
            return;
        end
    end
    error("error:  unknown varname %s\n", varname);
end



