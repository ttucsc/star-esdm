function model_join_minifiles(models, scenarios, varnames, ensembles, region, do_overwrite, runix, dry_run)
%   setups up and runs ARRM_V2 for specified models, scenarios, varnames & ensembles.
%       By default will run tasmin & tasmax, for rcp85 and rcp45 and r1i1p1.
%
%   Currently set up to use livneh data, and downscale for CONUS.
%       
%   Inputs:
%       models          string array of models to run
%       scenarios       optional.  string array of scenarios to run.  ["rcp45, "rcp85"].
%       varnames        optional.  string array of varnames to run.   ["tasmax","tasmin"]
%       ensembles       optional.  string array of ensembles to run.  ["r1i1p1"]
%       do_overwrite    optional.  set to true to allow overwrite of existing file.
%       runix           optional.  specify a list of indexes of run_table to be run.
%       dry_run         optional.  if true, prints out run_table so you can verify the correct runix indexes.
%
%   Edit as needed!  Generally will need to change template and directories.
%
    if (~exist("scenarios","var")),                             scenarios = ["rcp45","rcp85"];      end
    if (~exist("ensembles","var")),                             ensembles = "r1i1p1";               end
    if (~exist("varnames","var")),                              varnames  = ["tasmin","tasmax"];    end
    if (~exist("do_overwrite","var") || isempty(do_overwrite)), do_overwrite = false;               end
    if (~exist("runix",   "var")),                              runix = [];                         end
    if (~exist("dry_run","var")   || isempty(dry_run)),     dry_run = false;                        end
    
    
    run_table=valid_models("models",models, "scenarios", scenarios, "varnames", varnames,"ensembles",ensembles);
        
    nruns = size(run_table,1);
    if (isempty(runix))
        runix = 1:nruns;
    else
        runix = runix(runix >= 1 & runix <= nruns);
    end    
    
    if (dry_run)
        fprintf("\n%s: dry-run:  full run would run the following:\n", mfilename);
        disp(run_table(runix,:));
        return;
    end
        
    for r = runix
        
        model    = run_table.model{r};
        ensemble = run_table.ensemble{r};
        varname  = run_table.varname{r};
        scenario = run_table.scenario{r};
            
        basedir = "/Volumes/lacie_1/data/downscaled/arrm_v2";
        subdir  = sprintf("test/%s/%s/%s/%s", model,ensemble,varname,scenario);
        srcdir  = fullfile(basedir, subdir);
        outdir  = fullfile(basedir,"livneh");
        
        if (~isfolder(srcdir))
            error("error:  source dir %s does not exist", srcdir);
        end
        if (~isfolder(outdir))
            fprintf("NOTE:  output folder %s does not exist.  Creating.\n", outdir);
            setup_dirs(basedir,outdir)
        end
    
        templ1=sprintf("downscaled.%s.%s.%s.%s.ncamerica.1950.2100.grid_*.nc", model,ensemble, varname,scenario);
        template=fullfile(srcdir,templ1);
        
        outbase = sprintf("downscaled.livneh_test.%s.%s.%s.%s.conus.1950.2100.nc", model, ensemble, varname, scenario);
        outname=fullfile(outdir,outbase);
        
        if (~do_overwrite && isfile(outname))
            fprintf("warning:  file %s already exists;  set do_overwrite to true to overwrite it\n", outname);
        else
            finfo=dir(template);
            if (isempty(finfo))
                fprintf("warning:  cannot find any files matching input template %s.  Skipping this set\n", template);
                continue
            else
                fprintf("\njoining %4d files from : %s\n", length(finfo), srcdir);
                fprintf( "writing to              : %s\n\n", outname);
            end
        
            join_minifiles(template,varname,outname,"verbose",true, "overwrite", do_overwrite,'format','netcdf4','rotate',false,"clean_metadata",true);
        end
    end
end