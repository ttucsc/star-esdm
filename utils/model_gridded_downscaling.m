function retval = model_gridded_downscaling(obs_src, model, ensemble, scenario, varname, model_set, varargin)
%    model_gridded_downscaling(obs_src, model, ensemble, scenario, varname,
%    model_set, + optional keyword/value pairs for:
%               grgn, dry_run, region, latrange, lonrange, gridsize, obsdir, mdldir, do_parfor, do_clean_files, do_continue, do_extended)
%
%   Runs ARRM_V2 for specified models, scenarios, varnames & ensembles.
%       By default will run ARRM_V2_wrapper with inputs supplied.
%       Alternately, a 'dry-run' will print out the settings for the run(s), but not do the actual run.
%       (dry-run can be used to make sure you're getting the run you want before actually doing the run.
%
%       This function calls ARRM_V2_wrapper(...) to get the downscaling done, passing in the info needed to find
%       the desired data.
%
%   The best way to run this for lots of models, ensembles, scenarios, etc., is by running 
%   
%       setup_slurm_run(...)     which creates HPCC qsub scripts to run each model/ensemble/scenario/variable and join
%                               the output files into a single netcdf file.
%   and then sbatch the first script of each set of scripts.  Those scripts will email you as they start and complete.
%
%   Inputs:
%     required:
%
%       obs_src         string          observation source:  currently only "nrcan" or "livneh" or "nclimgrid"
%                                           obs_src identifies which gridded observation set to use
%       model           string    model to run.  
%                                           defines  model to be downscaled against obs_src
%       scenario        string array    scenario to run
%       varname         string array    variable to run
%       ensemble        string array    ensemble to run
%       model_set       string          "cmip5" or "cmip6".  Default:  "cmip5"
%                                           will specify the directories of where to find the model data
%     Optional keyword/value pairs
%       "grgn"          string          gridding key (one of "gr","gn", "gr1" etc);  for cmip6 runs.
%                                               if not present, will run all possible grgn's for the sets.
%                                               However, on HPCC, may run beyond 48-hour limit if more than 1 grgn is
%                                               being run.
%       onhost          string          "neys","nocona","quanah" or "xlquanah","icsf".  Name of system where run will be done.
%                                               defaults to "nocona".  (buty cannot check for files present on HPCC if
%                                               run on Neys.)
%       dry_run         true/false      [false]  if true, prints all params that would be run, but doesn't do actual run
%                                                   NOTE:  can't check existence of folders and files if not on same host as onhost.
%       region          string          "conus", "concan","canada","livneh","nrcan", "nclimgrid", "sheffield", "central", "test_25", or "test". see below
%                                           region defaults to the region defined by obs_src
%                                           region is used to label the output, and defines the lat/lon area used
%                                           (unless a latrange and lonrange is provided).  
%       latrange        [lat1, lat2]    range of latitudes to run.  Default based on obs_src or region
%       lonrange        [lon1, lon2]        will override the default lats & lons normally  specified by obs_src or region.
%       gridsize        [deglat,deglon] gridbox size, in degrees.  defaults to [2,2]
%                                           each child process will downscale all obs gridpoints within the gridbox.
%                                           Size of gridbox defines how much model data to read, and hence how much
%                                           memory each child process will need.
%       obsdir          string          directory of where to find obs data.
%       mdldir          string          directory of where to find model data
%       baseout         string          base directory of where to write output. 
%                                           no need to set these on neys or HPCC if using the standard data and locations
%                                           see code in init_params(...) for default values
%                                           or do a dry-run to see where it will look or write data                                           
%       do_parfor       true/false      run in parallel or not.  default:  true for neys & hpcc, false otherwise
%                                           if you've got a multi-core system, set do_parfor to true to use all cores.
%       cluster         'local'         if do_parfor is true, which cluster to run on.
%                    or 'redraider R2020b'            default:  'redraider R2020b' for HPCC (starts 'redraider R2020b' cluster) or 'local' for all other.
%                    or 'threads'       to run in threaded mode instead of child processes. (experimental)
%                                           NOTE:  'threads' won't work on current ARRM_V2 code.  icsf 1/22
%       nworkers        number          # of worker threads to start. Default is max possible on target system
%                                           18 on Neys, 120 on HPCC with 'local', or 128*nnodes for nocona 
%                                           generally no reason to specify this unless you want to run with fewer nodes for some reason.
%       nnodes          1-4             # of nodes to run on.  Default:  1
%                                           generally best to run on a single node, since all models will currently run in under 48 hour limit on a single node on nocona.
%       mb_per_core     2000 - 515625   Amount of memory, in MBytes, to use per node.  Applies only for runs using MPS at HPCC.
%                                           
%       do_clean_files  true/false      If true, erases all files in destination (output folder) before starting.  default:  true (false if do_continue is true).
%       do_continue     true/false      If true, restarts incomplete run;  will check for existing completed grid boxes
%                                           and only run those gridboxes which don't exist or aren't complete.
%                                           default:  false
%       do_extended     true/false      If true, adds a z-score variable for each data point to the output netcdf file.
%       use_scratch     true/false      If true:  use model files from /lustre/scratch/iscottfl.        default:  false.
%                                          false: use model files from /lustre/research/khayhoe
%       stnnums    v     [1,2,...]       vector of stations numbers to use (works with test_stations only)
%       hist_only       true/false      if true, runs entire hist period as 1 rolling step;  works only with
%                                           test_stations "region"
%       is_tll          true/false      defines input file name.  true:  "tllc" (time/lat/lon ordering, chunked).
%                                           false: "llt"  (lon/lat/time ordering).
%                                           default:  true for cmip6 files, false for cmip5 files.
%        exit_on_error  true/false      true:  aborts run if error occurs calculating results for a gridpoint;
%                                       false: assignes NA's for gridpoint if error encountered.
%                                           use false for 'production runs', but check for 'error' in the log files.
%
%       Region:  default set by obs_src as follows:
%                       "conus":        lat [20, 50],   lon [-130, -55],  Continental US, which gets parts of mexico & canada 
%                                                                               as well, but not Alaska or Hawaii, Marianas, Puerto Rico, US Virgins, etc..
%                       "livneh":       lat [24, 54],   lon [-126,-66],   which bounds livneh data
%                       "nclimgrid":    lat [24, 50],   lon [-125, -67]   which bounds nclimgrid data 
%                       "sheffield":    lat [-60,90],   lon [   0, 360]   which bounds sheffield data 
%                       "nrcan"         lat [41, 83.5], lon [-141, -52]
%                       "canada"        (same as nrcan)    
%                       "concan"        lat [41, 71],   lon [-141, -52]   continental Canada, excluding north of NorthWest Territories.
%                       "central":      lat [28, 48],   lon [-120,-90]    west/central US.  Nevada through Dakotas, LA through mid-ND.  Useful for testing
%                       "test":         lat [32, 40],   lon [-102,-94]    Small region in central US.  Approx. Midland TX to Kansas City, MO
%                       "test_stations":  (various)                         test at locations given by 25 test cities.  Output written to mat files, one per city. 
%                                                                       
%
%   Edit as needed!  Generally will need to change obsdir and obsnames below for other gridded obs data.
%   
%   9/12/2021:  added support for nclimgrid gridded obs data.  icsf
%   1/20/22:    updated code to support local cluster on HPCC's nocona group.  
%   1/20/22:    removed all references to pdf_map_method,  Proper defaultrs are now set for that later in the chain.
%                   NOTE:  can still pass pdf_map_method in;  it will be passed through in unMatched list.
%   2022-09-19: icsf    added support for Sheffield gridded obs data.
%
    if (nargin < 6)
        retval = 2;
        return;
    end
    
    retval = 1;
    try
        [obs_src, model, ensemble, scenario, varname, grgn, model_set, prcp_distrib, region, ...
         latrange, lonrange, gridsize, obsdir, histdir, mdldir, baseout, ...
         testsites, lats, lons, obs_yrs, hist_yrs, mdl_yrs, hist_only, is_tll, ...
         par_info, do_clean_files, do_continue, do_extended, this_hostname, run_on_hostname, dry_run, hist_styr, exit_on_error, isPrecipRun, Unmatched]  = init_params(obs_src, model, ensemble, scenario, varname, model_set, varargin{:});

    %------------ main program begins here------------

        run_table=valid_models("model_set", model_set, "models",model,"ensembles",ensemble, "varnames", varname, "scenarios", scenario,"grgns",grgn);
        
        if (size(run_table,1)==0)
            error("no matching data for %s %s %s %s %s %s", model_set, model, ensemble, varname, scenario, grgn);
        elseif (size(run_table,1)>1)
            disp(run_table);
            error("error:  multiple matching sets for %s %s %s %s %s %s: \n", model_set, model, ensemble, varname, scenario, grgn);
        end

        interp_method = "bilinear";

        nruns = size(run_table,1);
        if (nruns == 0)
            error("model_gridded_downscaling:  bailing out.");
        elseif (nruns ~= 1)
            fprintf(2, "error:  multiple matches to inputs:\n");
            disp(run_table);
            error("model_gridded_downscaling:  bailing out.");
        end
            

        if (dry_run)
            fprintf("\n\ndry_run:  full run would run the following:\n");
        end

        if (do_clean_files)
            fprintf("\nwill delete .nc and .log files in output folders before starting run\n\n");
        else
            fprintf("\nnot deleting .nc or .log files in output folders\n\n");
        end


        model    = run_table.model{1};
        ensemble = run_table.ensemble{1};
        varname  = run_table.varname{1};
        scenario = run_table.scenario{1};
        if (strcmp(model_set,"cmip6"))
            grgn     = run_table.grgn{1};
        else
            grgn     = "";
        end

        if (strcmp(model_set,"cmip5"))
            if (is_tll)
                mdlnames(1)=sprintf("%s.%s.%s.hist.day.1900.2005.tllc.nc",model,ensemble,varname);
                mdlnames(2)=sprintf("%s.%s.%s.%s.day.2006.2100.tllc.nc",model,ensemble,varname,scenario);
            else
                mdlnames(1)=sprintf("%s.%s.%s.hist.day.1900.2005.llt.nc",model,ensemble,varname);
                mdlnames(2)=sprintf("%s.%s.%s.%s.day.2006.2100.llt.nc",model,ensemble,varname,scenario);
            end
        else
            if (is_tll)
                mdlnames(1)=fullfile("historical", sprintf("%s_day_%s_historical_%s_%s_%d-2014.tllc.nc",varname, model,           ensemble, grgn, hist_styr));
                mdlnames(2)=fullfile(scenario,     sprintf("%s_day_%s_%s_%s_%s_2015-2100.tllc.nc",        varname, model, scenario, ensemble, grgn));
            else
                mdlnames(1)=fullfile("historical", sprintf("%s_day_%s_historical_%s_%s_%d-2014.llt.nc",varname, model,           ensemble, grgn, hist_styr));
                mdlnames(2)=fullfile(scenario,     sprintf("%s_day_%s_%s_%s_%s_2015-2100.llt.nc",        varname, model, scenario, ensemble, grgn));
            end
        end
        histnames  = mdlnames(1);

        if (strcmp(obs_src,"nrcan"))
            obsnames = sprintf("nrcan.%s.llt.nc", varname);
            if (strlength(grgn)==0)
                outname=sprintf("downscaled.%s.%s.%s.%s.nrcan_1_12th.%s.1950.2100.LLGRID.nc", model, ensemble, varname, scenario, region);
            else
                outname=sprintf("downscaled.%s.%s.%s.%s.%s.nrcan_1_12th.%s.1950.2100.LLGRID.nc", model, ensemble, varname, scenario, grgn, region);
            end
            obsvname = varname;
        elseif (strcmp(obs_src,"livneh"))
            if (strcmp(varname,"tasmax"))
                obsvname = "Tmax";
            elseif (strcmp(varname,"tasmin"))
                obsvname = "Tmin";
            elseif (isPrecipRun)
                obsvname = "Prec";
            end
            obsnames=sprintf("livneh_1_16th_%s.1950.2015.llt.nc",obsvname);
            if (strcmp(model_set, "cmip5"))
                outname = sprintf("downscaled.%s.%s.%s.%s.livneh_1_16th.%s.1950.2100.LLGRID.nc",    model, ensemble, varname, scenario,       region);
            else                    
                outname = sprintf("downscaled.%s.%s.%s.%s.%s.livneh_1_16th.%s.1950.2100.LLGRID.nc", model, ensemble, varname, scenario, grgn, region);
            end
        elseif (strcmp(obs_src,"nclimgrid"))
            if (strcmp(varname,"tasmax"))
                obsvname = "tmax";
            elseif (strcmp(varname,"tasmin"))
                obsvname = "tmin";
            elseif (isPrecipRun)
                obsvname = "prcp";
            end
            obsnames=sprintf("nclimgrid.day.%s.1951.2019.llt.nc",obsvname);
            if (strcmp(model_set, "cmip5"))
                outname = sprintf("downscaled.%s.%s.%s.%s.nclimgrid.%s.1950.2100.LLGRID.nc",    model, ensemble, varname, scenario,       region);
            else                    
                outname = sprintf("downscaled.%s.%s.%s.%s.%s.nclimgrid.%s.1950.2100.LLGRID.nc", model, ensemble, varname, scenario, grgn, region);
            end
        elseif (strcmp(obs_src,"sheffield"))
            if (strcmp(varname,"tasmax"))
                obsvname = "tmax";
            elseif (strcmp(varname,"tasmin"))
                obsvname = "tmin";
            elseif (isPrecipRun)
                obsvname = "prcp";
            end
            obsnames=sprintf("sheffield_global_%s_day_1948-2116.tllc_time.nc",obsvname);
            if (strcmp(model_set, "cmip5"))
                outname = sprintf("downscaled.%s.%s.%s.%s.sheffield.%s.1950.2100.LLGRID.nc",    model, ensemble, varname, scenario,       region);
            else                    
                outname = sprintf("downscaled.%s.%s.%s.%s.%s.sheffield.%s.1950.2100.LLGRID.nc", model, ensemble, varname, scenario, grgn, region);
            end
        elseif (strcmp(obs_src,"gfdl"))
            error("oops.  this program isn't set up for GFDL testing!");
%                 if (strcmp(varname,"tasmax") && obs_src == "livneh")
%                     obsvname = "Tmax";
%                 elseif (strcmp(varname,"tasmin"))
%                     obsvname = "Tmin";
%                 elseif (isPrecipRun)
%                     obsvname = "Prec";
%                 end
%                 obsnames=sprintf("livneh_1_16th_%s.1950.2015.llt.nc",obsvname);
%                 if (strcmp(model_set, "cmip5"))
%                     outname = sprintf("downscaled.%s.%s.%s.%s.livneh_1_16th.%s.1950.2100.LLGRID.nc",    model, ensemble, varname, scenario,       region);
%                 else                    
%                     outname = sprintf("downscaled.%s.%s.%s.%s.%s.livneh_1_16th.%s.1950.2100.LLGRID.nc", model, ensemble, varname, scenario, grgn, region);
%                 end
        else
            error("error:  obs_src must be ""ncrcan"", ""nclimgrid"", ""sheffield"" or ""livneh""");
        end

        if (strcmp(region,"test_stations"))
            outname = strrep(outname, "LLGRID.nc","%s.mat");
            histnames = mdlnames(1);
            if (hist_only)
                mdl_yrs = hist_yrs;
                pdf_yrs = range(hist_yrs)+1;
                pdf_yrstep = pdf_yrs;
                mdlnames = mdlnames(1);
            else
                mdlnames = mdlnames(2);
                mdl_yrs = [2016,2100];
            end
        end

        if (strcmp(model_set,"cmip5"))
            final_baseout = baseout;
            outsubdir= sprintf("%s_%s_%s_%s_%s", model,ensemble,varname,scenario, obs_src);
            outdir   = fullfile(final_baseout, outsubdir);
        else
            final_baseout = fullfile(baseout, scenario);
            outsubdir= sprintf("%s_%s_%s_%s_%s_%s", model,ensemble,varname,scenario, grgn, obs_src);
            outdir   = fullfile(final_baseout, outsubdir);
        end

        if (dry_run)
            fprintf("dry_run:  run %d parameters:\n", 1);
            fprintf("model_set:     %s\n", model_set);
            fprintf("model:         %s\n", model);
            fprintf("ensemble:      %s\n", ensemble);
            fprintf("varname:       %s\n", varname);
            fprintf("obsvname:      %s\n", obsvname);
            fprintf("scenario:      %s\n", scenario);
            fprintf("region:        %s\n", region);
            if (strlength(grgn)>0)
                fprintf("grgn:          %s\n", grgn);
            end
            if (is_prcp_variable(varname))
                fprintf("prcp_distrib:  %s\n", prcp_distrib)
            end
            fprintf("obsnames:      %s\n", vec2string(obsnames,  "format","%s ","brackets",'[]'));
            fprintf("histnames:     %s\n", vec2string(histnames, "format","%s ","brackets",'[]'));
            fprintf("mdlnames:      %s\n", vec2string(mdlnames,  "format","%s ","brackets",'[]'));
            fprintf("outname:       %s\n", outname);
            fprintf("obsdir:        %s\n", obsdir);
            fprintf("histdir:       %s\n", histdir);
            fprintf("mdldir:        %s\n", mdldir);
            fprintf("outdir:        %s\n", outdir);
            fprintf("\n");
            fprintf("obs_yrs:       %s\n", vec2string(obs_yrs,  "format","%4d ","brackets",'[]'));
            fprintf("hist_yrs:      %s\n", vec2string(hist_yrs, "format","%4d ","brackets",'[]'));
            fprintf("mdl_yrs:       %s\n", vec2string(mdl_yrs,  "format","%4d ","brackets",'[]'));
            if (~isempty(testsites))
                fprintf("lats:          %s\n", vec2string(lats, "format","%9.4f ", "brackets",'[]'));
                fprintf("lons:          %s\n", vec2string(lons, "format","%9.4f ", "brackets",'[]'));
                fprintf("sites:         %s\n", vec2string([testsites(1:min(5,length(testsites)));"..."]));
                fprintf("pdf_yrstep:    %s\n", vec2string(pdf_yrstep, "format","%9.4f "));
                fprintf("pdf_yrs:       %s\n", vec2string(pdf_yrs, "format","%9.4f "));
            else
                fprintf("latrange:      %s\n", vec2string(latrange, "format","%9.4f ", "brackets",'[]'));
                fprintf("lonrange:      %s\n", vec2string(lonrange, "format","%9.4f ", "brackets",'[]'));
                fprintf("gridsize:      %s\n", vec2string(gridsize, "format","%2d ",   "brackets",'[]'));
                fprintf("do_continue:   %d\n", do_continue);
                fprintf("do_extended:   %d\n", do_extended);
            end
            fprintf("interp_method: %s\n", interp_method);
            fprintf("do_parfor:     %d\n", par_info.do_parfor);
            if (par_info.do_parfor)
                fprintf("cluster:       %s\n", par_info.cluster);
                fprintf("nnodes:        %d\n", par_info.nnodes);
                fprintf("nworkers:      %d\n", par_info.nworkers);
                if (~isempty(par_info.mb_per_core))
                    tot_mem = par_info.mb_per_core * par_info.nworkers / par_info.nnodes / 1024.0;
                    fprintf("mb_per_core:   %d (%.2fGB total memory)\n", par_info.mb_per_core, tot_mem);
                end
            end
            fprintf("do_clean_files %d\n", do_clean_files);
            fprintf("exit_on_error: %d\n", exit_on_error);
            fprintf("is_tll:        %d\n", is_tll);
            fprintf("\n");
            flds = fieldnames(Unmatched);
            nflds = length(flds);
            if (nflds > 0)
                fprintf("Unmatched:\n");
                for i=1:nflds
                    if (length(Unmatched.(flds{i}))==1)
                        fprintf("\t%-12s: %s\n", flds{i}, string(Unmatched.(flds{i})));
                    else
                        fprintf("\t%-12s: %s\n", flds{i}, vec2string(Unmatched.(flds{i}),"Brackets",'[]'));
                    end

                end
                fprintf("\n");
            end


                % check for existence of data files and folders
            if (run_on_this_machine(run_on_hostname, this_hostname))
                for i=1:length(obsnames)
                    if (~isfile(fullfile(obsdir,  obsnames(i))))  
                        fprintf(2, "missing  obs file:     %s\n", fullfile(obsdir, obsnames(i))); 
                    else
                        fprintf(1, "obs file exists:       %s\n", fullfile(obsdir, obsnames(i))); 
                    end
                end
                for i=1:length(histnames)
                    if (~isfile(fullfile(histdir, histnames(i))))
                        fprintf(2, "missing hist file:     %s\n", fullfile(histdir, histnames(i))); 
                    else
                        fprintf(1, "hist file exists:      %s\n", fullfile(histdir, histnames(i))); 
                    end
                end
                for i=1:length(mdlnames)
                    if (~isfile(fullfile(mdldir, mdlnames(i))))
                        fprintf(2, "missing  mdl file:     %s\n", fullfile(mdldir, mdlnames(i))); 
                    else
                        fprintf(1, "mdl file exists:       %s\n", fullfile(mdldir, mdlnames(i))); 
                    end
                end

                if (~isfolder(final_baseout))
                    fprintf(2, "missing output folder: %s\n", final_baseout); 
                else
                    fprintf(1, "output folder exists:  %s\n", final_baseout); 
                end

                retval = 0;
            else
                fnames = strings(0,0);
                for i=1:length(obsnames)
                        fprintf(1, "obs file:       %s\n", fullfile(obsdir, obsnames(i))); 
                        fnames(end+1) = fullfile(obsdir, obsnames(i)); %#ok<AGROW>
                end
                for i=1:length(histnames)
                        fprintf(1, "hist file:      %s\n", fullfile(histdir, histnames(i))); 
                        fnames(end+1) = fullfile(histdir, histnames(i)); %#ok<AGROW>
                end
                for i=1:length(mdlnames)
                        fprintf(1, "mdl file:       %s\n", fullfile(mdldir, mdlnames(i))); 
                        fnames(end+1) = fullfile(mdldir, mdlnames(i)); %#ok<AGROW>
                end
                fprintf(1, "output folder:  %s\n", final_baseout);
                fprintf("\nto check for files and output folder on target system: \n");
                fprintf("\t ls -l ");
                fprintf("%s ", fnames);
                fprintf("\n\t ls -ld %s\n", final_baseout); 
                fprintf("\n");

            end
        else 
            setup_dirs(final_baseout, outsubdir)

                % check for existence of data files and folders
            ok=true;
            for i=1:length(obsnames)
                if (~isfile(fullfile(obsdir,  obsnames(i))))  
                    fprintf(2, "missing  obs file:     %s\n", fullfile(obsdir, obsnames(i))); 
                    ok=false;
                end
            end
                for i=1:length(histnames)
                    if (~isfile(fullfile(histdir, histnames(i))))
                        fprintf(2, "missing hist file:     %s\n", fullfile(histdir, histnames(i))); 
                        ok=false;
                    end
                end
                for i=1:length(mdlnames)
                    if (~isfile(fullfile(mdldir, mdlnames(i))))
                        fprintf(2, "missing  mdl file:     %s\n", fullfile(mdldir, mdlnames(i))); 
                        ok=false;
                    end
                end

                if (~isfolder(final_baseout))
                    fprintf(2, "missing output folder: %s\n", final_baseout); 
                    ok=false;
                end
                if (~ok)
                    error("Error:  missing files or folders");
                end

                if (do_clean_files), clean_dirs(final_baseout, outsubdir); end

                if (strcmp(region, "test_stations"))
                    [~,retval] = ARRM_V2_wrapper(model,ensemble,varname,scenario,obs_yrs,hist_yrs,mdl_yrs,...  
                                     "obsvname",obsvname, ...
                                     "obsnames", obsnames, ...
                                     "histnames",histnames, ...
                                     "mdlnames", mdlnames, ...
                                     "region", region, ...
                                     "testsites",   testsites, ...
                                     "outname", outname, ...
                                     "obsdir", obsdir,...
                                     "histdir", histdir,...
                                     "mdldir", mdldir,...
                                     "outdir", outdir, ...
                                     "prcp_distrib", prcp_distrib, ...
                                     "lats",    lats, "lons",lons, ...
                                     "llgrid_size",gridsize, ...
                                     "pdf_yrstep",  pdf_yrstep, ...     % test_stations-specific!
                                     "pdf_yrs",     pdf_yrs, ...        % test_stations-specific.  either doing all hist as 1 period, or all future as 1 period.
                                     "do_continue",do_continue, ...
                                     "extended",do_extended, ...
                                     "interp_method",interp_method, ...
                                     "do_parfor",par_info.do_parfor, "cluster",par_info.cluster, "nnodes", par_info.nnodes, "mb_per_core", par_info.mb_per_core, "nworkers", par_info.nworkers, ...
                                     "exit_on_error",exit_on_error, ...
                                     Unmatched...
                                     );
                else
                    [~,retval] = ARRM_V2_wrapper(model,ensemble,varname,scenario,obs_yrs,hist_yrs,mdl_yrs,...  
                                     "obsvname",obsvname, ...
                                     "obsnames", obsnames, ...
                                     "histnames",histnames, ...
                                     "mdlnames", mdlnames, ...
                                     "region", region, ...
                                     "outname", outname, ...
                                     "obsdir", obsdir,...
                                     "histdir", histdir,...
                                     "mdldir", mdldir,...
                                     "prcp_distrib", prcp_distrib, ...
                                     "outdir", outdir, ...
                                     "latrange", latrange, "lonrange",lonrange, ...
                                     "llgrid_size",gridsize, ...
                                     "do_continue",do_continue, ...
                                     "extended",do_extended, ...
                                     "interp_method",interp_method, ...
                                     "do_parfor",par_info.do_parfor, "cluster",par_info.cluster, "nnodes", par_info.nnodes, "mb_per_core", par_info.mb_per_core, "nworkers", par_info.nworkers, ...
                                     "exit_on_error",exit_on_error,...
                                     Unmatched...
                                     );
                end
                fprintf("ending model_gridded_downscaling().  ARRM_V2_wrapper returned %d\n", retval);
        end
    catch me
        report_me_error(me, mfilename(), 1, false);
        retval = 3;
    end
end
    
function [obs_src, model, ensemble, scenario, varname, grgn, model_set, prcp_distrib, region, ...
          latrange, lonrange, gridsize, obsdir, histdir, mdldir, baseout, ...
          testsites, lats, lons, obs_yrs, hist_yrs, mdl_yrs, hist_only, is_tll, ...
          par_info, do_clean, do_continue, do_extended, hostname, onhost, dry_run, hist_styr, exit_on_error, isPrecipRun, Unmatched]  = init_params(obs_src, model, ensemble, scenario, varname, model_set, varargin)

%      fprintf("in model_gridded_downscaling: init_params(...)\n");
%      disp(varargin);
%      
    p = inputParser;
    p.StructExpand = true;
    p.KeepUnmatched = true;
    
    model   = to_row(string(model));
    
    lats=[];
    lons=[];
    
    addRequired (p,"obs_src",               @(s) any(strcmp(s,["livneh","sheffield","nrcan","nclimgrid"])));
    addRequired (p,"model",                 @(s) isstring(s) || ischar(s));
    addRequired (p,"ensemble",              @(s) isstring(s) || ischar(s));
    addRequired (p,"scenario",              @(s) isstring(s) || ischar(s));
    addRequired (p,"varname",               @(s) isstring(s) || ischar(s));
    addRequired (p,"model_set",             @(s) any(strcmp(s,["cmip5","cmip6"])));
    addParameter(p,'prcp_distrib',strings(0),@(s) isempty(s) || strlength(s)==0 || any(strcmp(s,["generalizedpareto", "inversegaussian", "loglogistic", "lognormal","log","pwr"])));
    addParameter(p,'grgn',     strings(0),  @(s) (isempty(s) || isstring(s) || ischar(s)));
    addParameter(p,"region",   strings(0),  @(s) isempty(s)  || ((isstring(s) || ischar(s)) && any(strcmp(s,["conus", "concan","canada","livneh","nrcan", "nclimgrid", "global", "sheffield", "central","test","test_stations"]))));
    addParameter(p,'latrange', strings(0),  @(s) isnumeric(s) && (isempty(s) || length(s)==2));
    addParameter(p,'lonrange', strings(0),  @(s) isnumeric(s) && (isempty(s) || length(s)==2));
    addParameter(p,'gridsize', [2,2],       @(s) isnumeric(s) && length(s)==2);
    addParameter(p,"obsdir",   strings(0),  @(s) isstring(s) || ischar(s) || isempty(s));
    addParameter(p,"mdldir",   strings(0),  @(s) isstring(s) || ischar(s) || isempty(s));
    addParameter(p,"baseout",  strings(0),  @(s) isstring(s) || ischar(s) || isempty(s));
    addParameter(p,"dry_run",  false,       @(s) islogical(s) || (isnumeric(s) && any(s==[0,1])));
    addParameter(p,"do_parfor",true,        @(s) islogical(s) || any(s==[0,1]) || isempty(s));
    addParameter(p,"cluster", strings(0),   @(s) isstring(s) && any(strcmp(s,["local","nocona","redraider R2020b","threads"])));
    addParameter(p, "nworkers", [],         @(s) isempty(s) || (s>=1 && s<=512)); 
    addParameter(p, "nnodes", 1,            @(s) isnumeric(s) && mod(s,1==0) && s>=1 && s<=4); 
    addParameter(p, "mb_per_core",[],       @(s) isempty(s) || (isnumeric(s) && s>=0 && s<=515565));    % applies only to HPCC runs.  default will allow 48 notes in 250GB or 96 workers in 500G
    
    addParameter(p,"do_clean_files",true,   @(s) islogical(s) || any(s==[0,1]));
    addParameter(p,"do_continue",false,     @(s) islogical(s) || any(s==[0,1]));
    addParameter(p,"do_extended",false,     @(s) islogical(s) || any(s==[0,1]));
    addParameter(p,'mdl_yrs',[1950,2100],   @(s) isnumeric(s) && length(s)==2);
    addParameter(p,'hist_yrs',[1950,2000],  @(s) isnumeric(s) && length(s)==2);
    addParameter(p,'hist_styr',1950,        @(s) isnumeric(s) && length(s)==1 && s>=1850 && s<=1975);  % there are one or two files with different start years...
    addParameter(p,"use_scratch",false,     @(s) islogical(s) || any(s==[0,1]));
    addParameter(p,"hist_only",  false,     @(s) islogical(s) || any(s==[0,1]));
    addParameter(p,"stnnums",    [],        @(s) isnumeric(s));
    addParameter(p,"onhost",    strings(0), @(s) isempty(s) || (isstring(s) && (any(strcmpi(s,["neys","nocona","quanah","xlquanah","icsf-lmac","icsf-kmac"])) || strncmp(s,"cpu-",4))));
    addParameter(p,"is_tll",  [],           @(s) isempty(s) || islogical(s) || any(s==[0,1]));
    addParameter(p,"chunked", [],           @(s) isempty(s) || islogical(s) || any(s==[0,1]));  % chunked renamed to is_tll, but some scripts still use chunked, so added this back in.
    addParameter(p,"exit_on_error",false,   @(s) islogical(s) || any(s==[0,1]));
    
    if (length(varargin)<4), varargin{end+1:4}=[]; end
    parse(p,obs_src, model, ensemble, scenario, varname, model_set,varargin{:});
    obs_src     = string(p.Results.obs_src);
    model       = string(p.Results.model );
    ensemble    = string(p.Results.ensemble);
    scenario    = string(p.Results.scenario);
    varname     = string(p.Results.varname);
    grgn        = string(p.Results.grgn);
    model_set   = string(p.Results.model_set);
    prcp_distrib= string(p.Results.prcp_distrib);
    region      = string(p.Results.region);
    latrange    = p.Results.latrange;
    lonrange    = p.Results.lonrange;
    mdl_yrs     = p.Results.mdl_yrs;
    hist_yrs    = p.Results.hist_yrs;
    hist_styr   = p.Results.hist_styr;
    gridsize    = p.Results.gridsize;
    obsdir      = p.Results.obsdir;
    histdir     = p.Results.mdldir;
    mdldir      = p.Results.mdldir;
    baseout     = p.Results.baseout;
    dry_run     = p.Results.dry_run;
    do_parfor   = p.Results.do_parfor;
    nnodes      = p.Results.nnodes;
    mb_per_core = p.Results.mb_per_core;
    nworkers    = p.Results.nworkers;
    pcluster    = p.Results.cluster;
    do_clean    = p.Results.do_clean_files;
    do_continue = p.Results.do_continue;
    do_extended = p.Results.do_extended;
    use_scratch = p.Results.use_scratch;
    hist_only   = p.Results.hist_only;
    stnnums     = p.Results.stnnums;
    onhost      = p.Results.onhost;
    is_tll      = p.Results.is_tll;
    chunked     = p.Results.chunked;
    exit_on_error= p.Results.exit_on_error;

    if (do_continue) 		% don't clean files if do_continue is true.
        do_clean = false;
    end
    
    if (length(model)>1 || length(ensemble)>1 || length(scenario)>1 || length(varname)>1 || length(grgn)>1)
        error("error:  too many values for model/ensemble/scenario/varname/grgn");
    end
    
    obs_yrs     = hist_yrs;
    
    Unmatched = p.Unmatched;
    
    testsites = strings(0,0);
    
    [hostname, on] = get_hostname();
    
    hname = split(string(hostname),".");
    if (strcmp(hname(end),"local"))
        hostname = join(hname(1:end-1),".");
    end
    hostname = lower(hostname);
    
    if (isempty(onhost) || strlength(onhost)==0)
        if (on.hpcc_system)
            onhost = "nocona";
        else
            onhost = hostname;
        end
    end
    
    if (do_parfor)
        if (strcmp(onhost,"nocona") || strncmp(onhost,"cpu-",4))
            if (isempty(pcluster)), pcluster = "redraider R2020b"; end
            if (strcmp(pcluster,"redraider R2020b"))
                maxworkers = 120*nnodes;      
            else
                maxworkers = 120;
            end
            if (isempty(nworkers))
                nworkers = maxworkers;
            else
                nworkers = min(nworkers, maxworkers);
            end
        elseif (strcmp(onhost,"neys"))
            if (~isempty(pcluster) && ~any(strcmp(pcluster,["local","threads"]))), warning("cluster %s not available on neys; using local cluster", pcluster); end 
            maxworkers = 18;
        elseif (any(strcmp(onhost,["icsf-kmac","icsf-lmac"])))
            if (~isempty(pcluster) && ~any(strcmp(pcluster,["local","threads"]))), warning("cluster %s not available on kmac or lmac; using local cluster", pcluster); end
            maxworkers = 8;
        else
            error("error...don't know how many workers are available");
        end
                 
        if (isempty(nworkers))
            nworkers = maxworkers;                
        else
            nworkers = min(nworkers, maxworkers);
        end
    end
            
    par_info = struct("do_parfor", do_parfor, "cluster", pcluster, "nworkers", nworkers, "nnodes", nnodes, "mb_per_core", mb_per_core);
    
    isPrecipRun = is_prcp_variable(varname);
    if (isPrecipRun)
        if (isempty(prcp_distrib) || strlength(prcp_distrib) == 0)
            prcp_distrib = "pwr";
        end
        if (strcmp(prcp_distrib,"pwr"))
            fprintf("using root-power scaling for precip\n");
%           prcp_distrib = "loglogistic";
        elseif (strcmp(prcp_distrib,"log"))
            fprintf("using log scaling\n")
        elseif (any(strcmp(prcp_distrib,["generalizedpareto", "inversegaussian", "loglogistic"]))) 
            fprintf("using %s scaling\n", prcp_distrib);
        else
            error("unknown value for pcrp_distrib: %s", prcp_distrib);      % shouldn't happen...
        end
    else
        if (~isempty(prcp_distrib))
            if (strlength(prcp_distrib)==0)
                prcp_distrib = strings(0);
            else
                error("prcp_distrib specified for a non-precip run");
            end
        end
    end
    if (isempty(region) && isempty(latrange) && isempty(lonrange))
        region = obs_src; 
    end
    
    if (isempty(region) || strlength(region)==0)
        region = "custom";
    end
    if ((isempty(latrange) || isempty(lonrange)) && isempty_s(region))
        error("error:  custom region, but latrange or lonrange not specified");
    elseif (strcmpi(region,"conus"))
        latrange = [20,50];
        lonrange = [-130,-55];
        region = "conus";
        obs_src  = "livneh";
    elseif (any(strcmpi(region,["livneh","conus_livneh"])))
        latrange = [24,54];        % these covers the livneh data, which goes [25.15625,  52.84375] & [-124.59375,-67.03125] in steps of  .0625 (1/16 deg)
        lonrange = [-126,-66]; 
        region = "conus_livneh";            
        obs_src = "livneh";
    elseif (strcmpi(region,"nclimgrid"))%
        latrange = [24,50];        % these covers the livneh data, which goes [25.15625,  52.84375] & [-124.59375,-67.03125] in steps of  .0625 (1/16 deg)
        lonrange = [-125,-67]; 
        region = "nclimgrid";            
        obs_src = "nclimgrid";
    elseif (any(strcmpi(region,["nrcan", "canada"])))
        latrange = [41,83.5];         % these covers the nrcan data, which goes 41,83.41667] & [-141,-52.0833333] in steps of  .0833333 (1/11 deg)
        if (isempty(lonrange))
            lonrange = [-141,-52];
        end
        gridsize = min([2,2],gridsize); 
        region = "nrcan";
        obs_src = "nrcan";
    elseif (strcmpi(region,"central"))
        latrange = [28,48];
        lonrange = [-120,-90];
        region = "central";
%       obs_src  = "livneh";
    elseif (strcmpi(region,"test_stations"))
            S = load("test_stations_25.mat","siteTbl");
            siteTbl = S.siteTbl;
            if (~isempty(stnnums))
                siteTbl = siteTbl(stnnums,:);
            end
            latrange = [];
            lonrange = [];
            lats = siteTbl.lat;
            lons = siteTbl.lon;
            testsites=string(siteTbl.stnName);
            par_info.do_parfor = false;
    elseif (strcmpi(region,"test"))
        if (isempty(latrange))
            latrange = [32,40];
        end
        if (isempty(lonrange))
            lonrange = [-102,-94];
        end
        region = "test_region";
%       obs_src  = "livneh";
    elseif (any(strcmpi(region,["global","sheffield"])))
        latrange = [-90,90];
        lonrange = [0,359.99999];
        region = "global";
    end
    
    if (~strcmp(region,"test_stations") && (all(lonrange==[0,360]) || all(lonrange==[-180,180])))
        lonrange(end) = lonrange(end) - .00001;   % so modulus won't screw it up...
        region = "global";
    elseif (isempty(region))
        region = "custom";
    end

    if (isempty(is_tll)), is_tll = chunked; end         
    if (isempty(is_tll))
        if (strcmp(model_set,"cmip5"))
            is_tll = false;
        else
            is_tll = true;
        end
    end
    
    if (strcmp(onhost,"neys") || strcmp(onhost,"icsf-kmac"))
        if (isempty(mdldir) || strlength(mdldir)==0)
            if (strcmp(model_set,"cmip5"))
                mdldir  = "/Volumes/lacie_1/data/gcm_cmip5/rotated_netcdf";
            else
                if (is_tll)
                    mdldir  = "/Volumes/lacie_1/data/cmip6_tllc";
                else
                    mdldir  = "/Volumes/lacie_1/data/gcm_cmip6";
                end
            end
        end
        if (isempty(histdir) || strlength(histdir)==0)
            histdir = mdldir;
        end
        if (strcmp(obs_src,"sheffield") && (isempty(obsdir) || strlength(obsdir)==0))
            obsdir  = "/Volumes/lacie_1/data/obs/sheffield";
        elseif (strcmp(obs_src,"livneh") && (isempty(obsdir) || strlength(obsdir)==0))
            obsdir  = "/Volumes/lacie_1/data/obs/livneh";
        elseif (strcmp(obs_src,"nclimgrid") && (isempty(obsdir) || strlength(obsdir)==0))
            obsdir  = "/Volumes/lacie_1/data/obs/nclimgrid";
        elseif (strcmp(obs_src,nrcan) && (isempty(obsdir) || strlength(obsdir)==0))
            obsdir  = "/Volumes/lacie_1/data/obs/nrcan_rotated";
        end
        if (isempty(baseout) || strlength(baseout)==0)
            baseout = fullfile("/Volumes/lacie_1/data/downscaled/arrm_v2/", model_set, region, obs_src);
        end
        
    elseif (any(strcmp(onhost,["nocona","quanah","xlquanah"])) || strncmp(onhost,"cpu-",4))
        if (isempty(mdldir) || strlength(mdldir)==0)
            if (use_scratch)
                if (strcmp(model_set,"cmip6"))
                    mdldir  = sprintf("/lustre/scratch/iscottfl/cmip6_tllc"); 
                else
                    mdldir  = sprintf("/lustre/scratch/iscottfl/cmip5_rotated"); 
                end
            else
                if (strcmp(model_set,"cmip6"))
                    if (is_tll)
                        mdldir  = sprintf("/lustre/research/hayhoe/cmip6_tllc"); 
                    else
                        mdldir  = sprintf("/lustre/research/hayhoe/cmip6_rotated");
                    end
                else
                    mdldir  = sprintf("/lustre/research/hayhoe/cmip5_rotated"); 
                end
            end
        end
        if (isempty(histdir) || strlength(histdir)==0)
            histdir = mdldir;
        end
        if (isempty(obsdir) || strlength(obsdir)==0)
            obsdir  = sprintf("/lustre/research/hayhoe/obs/%s", obs_src);
            if (any(strcmp(obs_src, ["livneh","nclimgrid"]))), obsdir=strcat(obsdir,"_llt"); end
        end
        if (isempty(baseout) || strlength(baseout)==0)
            baseout = fullfile("/lustre/scratch/iscottfl/downscaled", model_set, region, obs_src);
        end
    elseif (strcmpi(onhost,"icsf-lmac"))
        username=getusername();
        if (strcmp(model_set,"cmip6"))
            mdldir  = sprintf("/Users/%s/data/cmip6_tllc",username); 
        else
            mdldir  = sprintf("/Users/%s/data/cmip5_rotated",username); 
        end        
        if (isempty(histdir) || strlength(histdir)==0)
            histdir = mdldir;
        end
        if (strcmp(obs_src,"livneh") && (isempty(obsdir) || strlength(obsdir)==0))
            obsdir  = sprintf("/Users/%s/data/obs/livneh",username);
        else
            obsdir  = sprintf("/Users/%s/data/obs/nrcan_rotated",username);
        end
        baseout = fullfile(sprintf("/Users/%s/data/downscaled/arrm_v2/",username), model_set, region);
    else
        error("can't tell what machine you're on.  Please set mdldir, obsdir and baseout appropriately");
    end 
        
end
    
