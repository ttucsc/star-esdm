function setup_slurm_gridded_run(obs_src, model, ensemble, varname, scenario, varargin)
%        setup_slurm_gridded_run(obs_src, model, ensemble, varname, scenario, keyword/value pairs for :  model_set, do_clean, gridsize, mdl_yrs, split_jobs, email_addr, use_scratch, do_continue, do_extended, is_tll...)
%                                                                        ^
%                                                                        |
%                                                                   These last few parameters are optional keyword/value pairs.
%
%
%   sets up sheffield, livneh, nclimgrid or NRCAN downscaling run for CMIP5 or CMIP6 model.
%
%   creates a set of shell scripts to be run on HPCC via sbatch, or on Neys as a shell script which invokes matlab.  
%
%   NOTE 1: This needs to be run on the system where the scripts will be submitted from.
%   NOTE 2: These scripts need to be sbatch'ed from the folder they are written to, because the scripts are in pairs.
%           The first script sbatch's the second script, and has the path to the second script hard coded.
%           (This could be fixed, at least for hpcc, using readlink.  Contact Ian if this is a problem.)
%   NOTE 3:  Currenly only sets up runs for nocona (redraider) system.  Quanah currently requires some different parameters. 
%
%   For each downscaling run created, you only need to sbatch the first script.
%
%   For livneh, the first script downscales the entire CONUS area, and automatically sbatch's the second script, which
%   merges the output minifiles into a single netcdf file.
%
%   Same for nclimgrid.
%
%   For NRCAN:
%       if running on nocona, it does the same as for livneh:  1 script to downscale, which sbatch's the merge script.
% %       if running on quanah, it creates 2 downscaling scripts.  The first script downscales half of the longitude 
% %                             range, then automatically sbatch's a 2nd script to downscale the other half of the  
% %                             longitude range, which then sbatchs the script to merge all the minifiles back into a
% %                             single netcdf file.  This is because it takes too long to do the entire area with just
% %                             36 cores.       
%
%   Inputs:
%       obs_src     "sheffield", "livneh", "nclimgrid" or "nrcan" or "stations" or "stations_25", stations_100 or stations_1000
%       model       model name, in double quotes.  example:  "CCSM4"
%       ensemble    ensemble name.  example:  "r1i1p1"
%       varname     variable name.  example:  "tasmin"
%       scenario    scenario name.  example:  "rcp45"
%
%           if any of the above are empty strings (""), then a script will be created for each possible combination for
%           that parameter.  For example, if scenario is left blank on a cmip6 run, then 4 sets of scripts will be
%           created, one each for ssp126, ssp245, ssp370 and ssp585 (at present), assuming the input data is present for 
%           each of those.  The list of valid combinations is found in valid_models_cmip5.csv and valid_models_cmip6.csv 
%
%   Optional keyword/value input pairs:
%       "grgn",       "gr"/"gn"     gr or gn
%       "model_set",  "cmip5"/"cmip6"   default:  code attempts to determine from model name, ensemble or scenario.
%       "region",     region        string defining region (used to select lat/lon range).  Should be "conus_livneh" or
%                                       "nrcan" or "nclimgrid", or "global".  ["sheffield" == ""global""] 
%                                       If empty, defaults to obs_src value.
%       "gridsize",   gridsize      size of minifile's gridbox, in lat, lon degrees.  example:  [2,2] 
%       "mdl_yrs",    mdl_yrs       range of years to downscale.  example (& default):  [1950,2100]
%       "hist_yrs",   hist_yrs      range of years to use for historical period.  default:  [1950,2000]
%       pdf_map_method "linear"/"clim"  (& a few others;  see comments in ARRM_V2 code)
%                                   default is linear for temperature runs, clim for precip runs.  
%                                       (default should be satisfactory.)
%       "onhost"     "nocona" | "quanah" | "icsf-kmac"  (where to run, if different from current host.) ("hpcc" is alias for "nocona") 
%                                   default:  current host (or any hpcc core, if host begins with "cpu-".
%       "base_dir"    some_folder   base directory of where to write console output, slurm_out and slurm_err files
%                                       default:  on HPCC:  /lustre/scratch/USER
%                                                 on neys, etc: user's home folder 
%       "script_dir", some_folder   Write script to some_folder.  Default:  ~/[hostname]_scripts/[username]
%       "local_script_dir", some_folder  If script will run on a different host (usually hpcc), local folder of where to
%                                           write script files.  If empty and onhost is not current host, default is "./hpcc_scripts" 
%       "email_addr", email_addr    email address to send notifications to, in double quotes.
%                       It should pick up the email address automatically from your login.  See code at end of 
%                       function init_params(...), at end of program.
%       "do_clean",   true/false    If true, deletes any existing *.sh files in script_dir before writing new scripts.  default:  false
%       "exclusive",  true/false    If true, will specify that no other jobs can run on the node(s) requested.  Default:  true
%       "nworkers",   nworkers      Spawn no more than nworkers child processes (# of parallel jobs).  Default:  max for system  (128 Nocona, 36 quanah, 18 neys, 8 icsf-kmac]
%       "nnodes",     nnodes        # of nodes to run on .  Default:  1.  Applies only to HPCC system.
%       "mb_per_core", mb_per_core  Applies only to HPCC runs.  default:  use SLURM default for system (~4 GB*nworkers on nocona, 5.3*nworkers on quanah) 
%                                       NOTE:  max mem per core usage on Nocona, is 515565 MB / (nworkers / nworkers)
%       "use_scratch" true/false    If true,  looks for model input files in /lustre/scratch/iscottfl
%                                   if false, looks for model input files in /lustre/research/khayhoe  default:  true
%       "do_continue",true/false    defaults to false.  If true, will continue a partially completed run.  only runs gridboxes with no output file yet 
%                                       or which don't have completion_status set to 1 (completed) or 2 (completed, but no valid data).
%       "is_tll ",    true/false    defines filename endings (for variable dimension ordering & chunking.)  
%                                       true:  tllc (time/lat/lon, chunked)
%                                       false: llt  (lon/lat/time)
%                                       default:  true (tllc) for cmip6, false (llt) for cmip5
%       "use_dev_folder", true/false or foldername
%                                   by default, uses source code in folder where setup_slurm_run is run from.
%                                   (or /lustre/work/iscottfl/ARRM_V2 if being run on neys, icsf-kmac, etc.)
%                                       set to true to use /lustre/work/iscottfl/ARRM_V2_dev/ARRM_V2
%                                       can also be a string, specifying folder to run from if different. f
%       "use_MPS", true/false       If true, set up run for use with Matlab Parallel Server.
%                                       If false, sets up standard single-node parallel run on HPCC (or Neys)  
%       "do_extended", true/false   If true, Saves z-score values for each data point in a separate variable.
%       "sigma_normalize"           true/false.  Defaults to true for temp, false for precip.
%                                       far_outlier_thresh specifies prob level beyond which a simple scaling is done (outlier-correction)
%                                           (beyond which the statistical mapping is unreliable)
%                                       Set to 0 to turn off outlier-correction algorithm.
%                                   far_outlier_anchor_pt is the point use for scaling far outliers.
%                                       (ignored if far_outlier_thresh is 0) 
%       "far_outlier_thresh"        prob (0-0.5 or [0-.5, .5-1.0]).  Defaults to 1e-4  (.0001, .9999) for temp, .0228 (2-sigma) (.0228, .9772) for precip
%       "far outlier_anchor_pt"     prob (0-0.5 or [0-.5, .5-1.0]).  Defaults to .0228 (.0228, .9772) for temp, .1587 (1-sigma) (.1587, .8413) for precip
%                                  
%
%           these need updating, ian!
%   example:
%
%       setup_slurm_run("livneh","CCSM4","r1i1p1","tasmin","rcp45","nocona");
%   or
%       setup_slurm_run("livneh","CCSM4","r1i1p1","tasmin","rcp45", "nocona", "do_clean", true, "gridsize", [2,2], "mdl_yrs", [1950,2100]);
%
%           creates 2 sbatch scripts (just sbatch the first one):
%
%               down.livneh_CCSM4_r1i1p1_tasmin_rcp45.01.sh	
%               rejoin.livneh_CCSM4_r1i1p1_tasmin_rcp45.sh                  
%
%   or
%       setup_hpcc_run("nrcan","CCSM4","r1i1p1","tasmin","rcp45","quanah");
%   or
%       setup_hpcc_run("nrcan","CCSM4","r1i1p1","tasmin","rcp45", "quanah", "do_clean", true, "gridsize", [2,2], "mdl_yrs", [1950,2100]);
%
%           creates 3 sbatch scripts (just sbatch the first one):
%
%               down.nrcan_CCSM4_r1i1p1_tasmin_rcp45.01.sh	
%               down.nrcan_CCSM4_r1i1p1_tasmin_rcp45.02.sh	
%               rejoin.nrcan_CCSM4_r1i1p1_tasmin_rcp45.sh
%
%   slurm output (matlab's console output) will be written to folder: ~/slurm_out
%   slurm error messages will be written to folder:                   ~/slurm_err
%
%   ARRM_V2's log files (with similar output to the slurm output file) is
%   written to the same folder as the intermediate and final netcdf files.
%
%   The ARRM_V2 output folder and the slurm output and error folders will
%   be created if they do not exist;  an error is thrown if any of the
%   folders cannot be created.
%
%   9/12/2021:  Added support for nclimgrid gridded obs
%   3/10/2022:  added path change to use R2021b, and removed try_next option.  R2021b doesn't have the limit of 4 jobs max at a time 
%   2022-09-18  added sheffield as option for gridded obs data.
%_________________________________________________________________________


    if (nargin == 0), help(mfilename); return; end

    [obs_src, model, ensemble, varname, scenario, grgn, model_set, prcp_distrib, region, gridsize, mdl_yrs, hist_yrs, do_clean, job_lon_ranges, exclusive, nworkers, use_scratch, sourcecode_dir, base_outdir, nstns, do_continue, do_extended, is_tll, pdf_map_method, far_outlier_thresh, far_outlier_anchor_pt, sigma_normalize, use_MPS, nnodes, mb_per_core, onhost, email_addr, script_dir, local_script_dir, unMatched] = init_params(obs_src, model, ensemble, varname, scenario, varargin{:});    

    if (do_clean)
        for s=1:length(scenario)
            lsdir=fullfile(local_script_dir,scenario);
            if (~isfolder(lsdir))
                lsdir=local_script_dir;
            end
            downs = fullfile(lsdir, "down.*.sh");
            joins = fullfile(lsdir, "rejoin.*.sh");
            system(sprintf("rm %s %s 2> /dev/null", downs, joins));
        end  
    end

    if (isempty(scenario) || all(strlength(scenario)==0))
        [runtbl, csvname] = valid_models("models",model,"varnames", varname, "scenarios","-historical","ensembles", ensemble, "grgn", grgn, "model_set", model_set);
    else
        [runtbl, csvname] = valid_models("models",model,"varnames", varname, "scenarios",scenario,"ensembles", ensemble, "grgn", grgn, "model_set", model_set);
    end
    nsets = size(runtbl,1); % for now, should only be 1.
    
    if (nsets == 0)
        fprintf(2, "No matching %s sets found in %s for %s %s %s %s %s", model_set, csvname, model, varname, scenario, ensemble, grgn);
        fprintf(2, "to update %s with all %s combinations, run generate_valid_model_csv_files.sh\n", csvname, model_set)
        return;
    end
    
    for i=1:nsets        
        if (~isfolder(local_script_dir)), error("error:  local script folder %s does not exist", local_script_dir); end
        my_scenario = string(runtbl.scenario{i});
        my_script_dir = fullfile(script_dir,my_scenario);
        my_local_script_dir = fullfile(local_script_dir,my_scenario);
        if (~isfolder(my_local_script_dir))
            my_script_dir = script_dir;
            my_local_script_dir = local_script_dir;
        end

        
        if (strncmp(obs_src,"station",7))
            if (strncmp(obs_src,"station",7)), error("error:  station setup needs to be updated.  Contact Ian"); end    %abort for now if station run...this needs to be updated/tested. 
            [join_name, script_names] = create_station_run_file(runtbl(i,:), gridsize, mdl_yrs, hist_yrs, job_lon_ranges, obs_src, model_set, prcp_distrib, pdf_map_method, far_outlier_thresh, far_outlier_anchor_pt, sigma_normalize, exclusive, nworkers, use_scratch, email_addr, sourcecode_dir, base_outdir, nstns, do_continue, do_extended, is_tll, use_MPS, nnodes, mb_per_core, onhost, my_script_dir, my_local_script_dir, unMatched);
            create_station_join_file(runtbl(i,:), region, obs_src, gridsize, mdl_yrs, join_name, model_set, email_addr, do_extended, sourcecode_dir, base_outdir, onhost, my_local_script_dir);
        else
            [join_name, script_names] = create_gridded_run_file(runtbl(i,:), gridsize, mdl_yrs, hist_yrs, job_lon_ranges, obs_src, model_set, prcp_distrib, pdf_map_method, far_outlier_thresh, far_outlier_anchor_pt, sigma_normalize, exclusive, nworkers, use_scratch, email_addr, sourcecode_dir, base_outdir, do_continue, do_extended, is_tll, use_MPS, nnodes, mb_per_core, onhost, my_script_dir, my_local_script_dir, unMatched);
            create_gridded_join_file(runtbl(i,:), region, obs_src, gridsize, mdl_yrs, join_name, model_set, email_addr, do_extended, sourcecode_dir, base_outdir, onhost, my_local_script_dir);
        end
        fprintf("\n");
        for j=1:length(script_names)
            if (j==1)
                fprintf("sbatch script name:     %s\n", script_names(j));
            elseif (j==length(script_names))
                fprintf("join script name:       %s\n", script_names(j));
            else
                fprintf("next script name:       %s\n", script_names(j));
            end
        end
        
        [~, on] = get_hostname;
        if (on.hpcc_system && strcmp(onhost,"nocona"))
            fprintf("\nTo submit entire job:   sbatch %s\n", script_names(1));
        else
            fprintf("\nTo submit:  \n");
            fprintf("\tcopy scripts to hpcc: %s\n", my_script_dir)
            fprintf("\tand sbatch %s\n", fullfile(my_script_dir, basename(script_names(1)))); 
        end
    end
    
end

function     [join_name, script_names] = create_gridded_run_file(   tbl,      gridsize, mdl_yrs, hist_yrs,     lon_ranges, obs_src, model_set, prcp_distrib, pdf_map_method, far_outlier_thresh, far_outlier_anchor_pt, sigma_normalize, exclusive, nworkers, use_scratch, email_addr, sourcecode_dir, base_outdir, do_continue, do_extended, is_tll, use_MPS, nnodes, mb_per_core, onhost, script_dir, local_script_dir, unMatched)

    if (~isempty(unMatched)), error("can't handle unMatched yet"); end

    model = string(tbl.model(1));
    ensemble = string(tbl.ensemble(1));
    varname = string(tbl.varname(1));
    scenario = string(tbl.scenario(1));
    if (strcmp(model_set,"cmip6"))
        grgn     = string(tbl.grgn(1));
    else
        grgn     = "";
    end
    
    if (isempty(lon_ranges))
        njobs = 1;
    else
        njobs = size(lon_ranges, 1);
    end
    
    if (strcmp(onhost,"quanah"))  % we're on quanah nodes
        run_on_hpcc = true;
        partition = "quanah";
        queue="quanah";
        cluster = "quanah";
    elseif (strcmp(onhost,"nocona"))    % this needs updating with numerical list.  some cpu- systems are on quanah.
        run_on_hpcc = true;
        partition = "nocona";
        queue=strings(0);
        if (use_MPS)
            cluster = "redraider R2020b";
        else
            cluster="local";
        end
    else
        run_on_hpcc = false;
        partition="";
        queue=strings(0);
        cluster = "local";
    end
            % these are output folders for files generated by SLURM
            % redirection consold output and errors.  They are only used in
            % sbatrch scripts, not passed to the ARRM_V2 code.
    outfiles_dir = fullfile(base_outdir, "ARRM_V2_out");            
    slurmout_dir = fullfile(base_outdir, "slurm_out");
    slurmerr_dir = fullfile(base_outdir, "slurm_err");
    if (run_on_hpcc)
        output_dirs = [outfiles_dir, slurmout_dir, slurmerr_dir];
    else
        output_dirs = outfiles_dir;
    end
        
    
    script_names = strings(0,0);
    for i=1:njobs                           % njobs is 1 if we can run all gridpoints within SLURM's time limit (48 hrs on TTU HPCC).  
                                            % This is currently the case if we're running on nocona system, with 128
                                            % cores per node.
%         if (strcmp(model_set,"cmip5"))
        if (strlength(grgn)==0)
            run_lbl    = sprintf("%s_%s_%s_%s_%s_%s",    model_set, obs_src, model, ensemble, varname, scenario);
        else
            run_lbl    = sprintf("%s_%s_%s_%s_%s_%s_%s", model_set, obs_src, model, ensemble, varname, scenario, grgn);
        end

                        
        run_name       = fullfile(      script_dir, sprintf("down.%s.%02d.sh",   run_lbl, i));
        local_run_name = fullfile(local_script_dir, sprintf("down.%s.%02d.sh",   run_lbl, i));
        if (i==njobs)
            next_job_name       = fullfile(      script_dir, sprintf("rejoin.%s.sh", run_lbl));
            local_next_job_name = fullfile(local_script_dir, sprintf("rejoin.%s.sh", run_lbl));
        else
            next_job_name       = fullfile(      script_dir, sprintf("down.%s.%02d.sh",   run_lbl, i+1));
            local_next_job_name = fullfile(local_script_dir, sprintf("down.%s.%02d.sh",   run_lbl, i+1));
        end
        
        [~,job_name,~] = fileparts(run_name);

        fid = fopen(local_run_name, "w");
        system(sprintf("chmod +x %s", local_run_name));
        fprintf(fid,"#!/bin/bash\n");

        script_names(end+1) = local_run_name; %#ok<AGROW>
        
        if (run_on_hpcc)
            fprintf(fid,"#SBATCH --job-name=%s\n", job_name);
            fprintf(fid,"#SBATCH --mail-user=%s\n", email_addr);
            fprintf(fid,"#SBATCH --mail-type=BEGIN,END,FAIL\n");
            fprintf(fid,"#SBATCH --output=%s/%%x.o%%j\n", slurmout_dir);
            fprintf(fid,"#SBATCH --error=%s/err.%%x.%%j\n", slurmerr_dir);
            fprintf(fid,"#SBATCH --partition=%s\n", partition);
            if (~isempty(queue))
                fprintf(fid,"#SBATCH --queue=%s\n", queue);
            end
            if (use_MPS)
                fprintf(fid,"#SBATCH --nodes=%d\n", 1);         % if using MPS (cluster redraider R2020b), then we'll set the SLURM parameters in ARRM_V2_wrapper.
                fprintf(fid,"#SBATCH --ntasks=%d\n", 1);        % if using MPS (cluster redraider R2020b), then we'll set the SLURM parameters in ARRM_V2_wrapper.
            else
                fprintf(fid,"#SBATCH --nodes=%d\n", nnodes);    % if not using MPS. then we're using the local cluster, and we set the SLURM params here.         
                if (exclusive)
                    fprintf(fid, "#SBATCH --exclusive\n");
                    fprintf(fid, "#SBATCH --ntasks=%d\n", nworkers+1);       % add 1 for the parent process.
                else
                    mb_per_core = floor(min(mb_per_core, 5200 * nworkers/(nworkers+1)));    % this needs to be changed if each core needs more than 5.2 GB mem to do its work.
                    fprintf(fid, "#SBATCH --mem-per-cpu=%d\n", mb_per_core);
                    fprintf(fid, "#SBATCH --ntasks=%d\n", nworkers+1);       % add 1 for the parent process.
                end
            end
            
            fprintf(fid,"#SBATCH --chdir=%s\n", sourcecode_dir);

                % Make sure we have all the modules loaded that we need
%             fprintf(fid,"\n# Load the modules needed.  \n#The load_my_modules script loads the appropriate modules depending on the architecture\n\n");
%             fprintf(fid,"#source /home/iscottfl/bin/load_my_modules   # skipping this for now.  Matlab module may be sufficient...\n");
%             fprintf(fid,"#module load matlab\n");
            fprintf(fid, "\nexport PATH=$PATH:/opt/apps/snfs/RedRaider/matlab/R2021b/bin  # for now, to get R2021b, we're setting the path explicitly.  HPCC will create a module for it shortly.\n\n");
            fprintf(fid,"module list\n\n");
            fprintf(fid,"echo ""running on node $HOSTNAME""\n\n");

                % create the output folders if they don't exist
%             fprintf(fid,"\n# create the output folders if they don't exist\n");
%             fprintf(fid,"if [[ ! -d %s ]]; then mkdir %s; fi\n", base_dir, base_dir);
%             fprintf(fid,"if [[ ! -d %s ]]; then mkdir %s; fi\n", outfiles_dir, outfiles_dir);
%             fprintf(fid,"if [[ ! -d %s ]]; then mkdir %s; fi\n", slurmout_dir, slurmout_dir);
%             fprintf(fid,"if [[ ! -d %s ]]; then mkdir %s; fi\n", slurmerr_dir, slurmerr_dir);
%             fprintf(fid,"#\n\n");

            fprintf(fid,"curdir=$(pwd)\n");
            fprintf(fid,"\ncd %s\n\n", sourcecode_dir);

            if (isempty(lon_ranges))
                if (strcmp(model_set, "cmip5"))
                    fprintf(fid,"matlab -nodisplay -r 'cd %s ; retval=model_gridded_downscaling(""%s"",""%s"",""%s"",""%s"",""%s"",""%s"", ""prcp_distrib"",""%s"",                   ""gridsize"", [%d,%d], ""mdl_yrs"", [%d,%d], ""hist_yrs"", [%d,%d], ""pdf_map_method"", ""%s"", ""far_outlier_thresh"", %f, ""far_outlier_anchor_pt"", %f,  ""sigma_normalize"", %d, ""use_scratch"",%d, ""do_continue"", %d, ""do_extended"", %d, ""is_tll"", %d, ""onhost"", ""%s"", ""cluster"", ""%s"", ""nworkers"", %d, ""nnodes"", %d, ""mb_per_core"", %d ); quit(retval); ' > %s/out.${SLURM_JOB_ID}.${SLURM_JOB_NAME}\n", ...
                            sourcecode_dir, obs_src, model, ensemble, scenario, varname, model_set, prcp_distrib,       gridsize, mdl_yrs, hist_yrs, pdf_map_method, far_outlier_thresh, far_outlier_anchor_pt, sigma_normalize, use_scratch, do_continue, do_extended, chunked, onhost, cluster, nworkers, nnodes, mb_per_core, outfiles_dir);
                else
                    fprintf(fid,"matlab -nodisplay -r 'cd %s ; retval=model_gridded_downscaling(""%s"",""%s"",""%s"",""%s"",""%s"",""%s"", ""prcp_distrib"",""%s"", ""grgn"", ""%s"", ""gridsize"", [%d,%d], ""mdl_yrs"", [%d,%d], ""hist_yrs"", [%d,%d], ""pdf_map_method"", ""%s"", ""far_outlier_thresh"", %f, ""far_outlier_anchor_pt"", %f,  ""sigma_normalize"", %d, ""use_scratch"",%d, ""do_continue"", %d, ""do_extended"", %d, ""is_tll"", %d, ""onhost"", ""%s"", ""cluster"", ""%s"", ""nworkers"", %d, ""nnodes"", %d, ""mb_per_core"", %d); quit(retval); ' > %s/out.${SLURM_JOB_ID}.${SLURM_JOB_NAME}\n", ...
                            sourcecode_dir, obs_src, model, ensemble, scenario, varname, model_set, prcp_distrib, grgn, gridsize, mdl_yrs, hist_yrs, pdf_map_method, far_outlier_thresh, far_outlier_anchor_pt, sigma_normalize, use_scratch, do_continue, do_extended, is_tll, onhost, cluster, nworkers, nnodes, mb_per_core, outfiles_dir);
                end
            else
                if (strcmp(model_set, "cmip5"))
                    fprintf(fid,"matlab -nodisplay -r 'cd $s ; retval=model_gridded_downscaling(""%s"",""%s"",""%s"",""%s"",""%s"",""%s"", ""%s"", ""prcp_distrib"",""%s"",                   ""gridsize"",[%d,%d], ""mdl_yrs"", [%d,%d], ""hist_yrs"", [%d,%d], ""lonrange"",[%6g,%6g], ""pdf_map_method"", ""%s"", ""far_outlier_thresh"", %f, ""far_outlier_anchor_pt"", %f,  ""sigma_normalize"", %d, ""use_scratch"",%d, ""do_continue"", %d, ""do_extended"", %d, ""is_tll"", %d, ""onhost"", ""%s"", ""cluster"", ""%s"", ""nworkers"", %d, ""nnodes"", %d, ""mb_per_core"", %d); quit(retval); ' > %s/out.${SLURM_JOB_ID}.${SLURM_JOB_NAME}\n", ...
                            sourcecode_dir, obs_src, model, ensemble, scenario, varname, model_set, prcp_distrib,       gridsize, mdl_yrs, hist_yrs, lon_ranges(i,:), pdf_map_method, far_outlier_thresh, far_outlier_anchor_pt, sigma_normalize, use_scratch, do_continue, do_extended, is_tll, onhost, cluster, nworkers, nnodes, mb_per_core, outfiles_dir);
                else
                    fprintf(fid,"matlab -nodisplay -r 'cd %s ; retval=model_gridded_downscaling(""%s"",""%s"",""%s"",""%s"",""%s"",""%s"", ""%s"", ""prcp_distrib"",""%s"", ""grgn"", ""%s"", ""gridsize"",[%d,%d], ""mdl_yrs"", [%d,%d], ""hist_yrs"", [%d,%d], ""lonrange"",[%6g,%6g], ""pdf_map_method"", ""%s"", ""far_outlier_thresh"", %f, ""far_outlier_anchor_pt"", %f,  ""sigma_normalize"", %d, ""use_scratch"",%d, ""do_continue"", %d, ""do_extended"", %d, ""is_tll"", %d, ""onhost"", ""%s"", ""cluster"", ""%s"", ""nworkers"", %d, ""nnodes"", %d, ""mb_per_core"", %d); quit(retval); ' > %s/out.${SLURM_JOB_ID}.${SLURM_JOB_NAME}\n", ...
                            sourcecode_dir, obs_src, model, ensemble, scenario, varname, model_set, prcp_distrib, grgn, gridsize, mdl_yrs, hist_yrs, lon_ranges(i,:), pdf_map_method, far_outlier_thresh, far_outlier_anchor_pt, sigma_normalize, use_scratch, do_continue, do_extended, is_tll, onhost, cluster, nworkers, nnodes, mb_per_core, outfiles_dir);
                end
            end
            fprintf(fid,"\nretval=$?      # save return value from model_gridded_downscaling\n");
            fprintf(fid, 'echo "model_gridded_downscaling returned $retval"\n');
            
                % If we succeeded (retval==0), then run join_minifiles to create the combined netcdf file.
            
            fprintf(fid,"\nif (( $retval == 0 ))\n");
            fprintf(fid,"then\n");
            fprintf(fid,"#          current job succeeded, so submit script to join minifiles\n");
            fprintf(fid,"     echo ""running join script %s""\n", next_job_name);
            fprintf(fid,"     sbatch %s\n", next_job_name);
            
            fprintf(fid,"else\n");
%
%               No longer looping on failure.  Now, if first step fails, just output error message and exit.
            if (~run_on_hpcc)
                fprintf(fid,"     echo ""problems running model_gridded_downscaling.  Please check output.  follow-on script would be %s""\n", next_job_name);
            else
                fprintf(fid,"     echo ""problems running model_gridded_downscaling.  Please check output.  follow-on script would be %s""\n", next_job_name);
%                 THIS IS A KLUDGE TO DEAL WITH HPCC INFINIBAND PROBLEMS
%                 fprintf(fid,"#            current job failed.  Try resubmitting (up to max of 3 times)\n");
%                 maxtries=3;
%                 logfile="/home/iscottfl/nocona_scripts/submitted.log";
%                 fprintf(fid,"     echo ""current job failed.  Trying to resubmit"" \n");
%                 fprintf(fid,"     /home/iscottfl/bin/try_resubmit.sh %s %d %s\n\n", run_name, maxtries, logfile);       % minor bug here.  resubmit.sh doesn't seem to return the number of tries it has been run.  So it resubmits, regardless of how many tries.
%                 fprintf(fid,"     ntries=$?\n");
%                 fprintf(fid,"     echo ""try_resubmit returned $ntries .  Now sending text"" \n");
%                 
%                 fprintf(fid,"     /home/iscottfl/bin/sendtext ""job %s failed;  trying resubmit for $ntries time""\n", run_name);
            end
            fprintf(fid,"fi\n");
            fprintf(fid,"exit $retval\n");

            fclose(fid);
            
            if (isempty(lon_ranges))
                if (strcmp(model_set,"cmip5"))
                    model_gridded_downscaling(obs_src, model, ensemble, scenario, varname, model_set, "prcp_distrib",prcp_distrib, "dry_run",true,               "gridsize", gridsize, "mdl_yrs", mdl_yrs, "hist_yrs", hist_yrs,                             "pdf_map_method", pdf_map_method, "use_scratch", use_scratch, "do_continue", do_continue, "do_extended", do_extended, "is_tll", is_tll, "onhost", onhost, "cluster", cluster, "nworkers", nworkers, "nnodes", nnodes, "mb_per_core", mb_per_core, "output_dirs", output_dirs);
                else
                    model_gridded_downscaling(obs_src, model, ensemble, scenario, varname, model_set, "prcp_distrib",prcp_distrib, "dry_run",true,               "gridsize", gridsize, "mdl_yrs", mdl_yrs, "hist_yrs", hist_yrs,                             "pdf_map_method", pdf_map_method, "use_scratch", use_scratch, "do_continue", do_continue, "do_extended", do_extended, "is_tll", is_tll, "onhost", onhost, "cluster", cluster, "nworkers", nworkers, "nnodes", nnodes, "mb_per_core", mb_per_core, "output_dirs", output_dirs);
                end
            else
                if (strcmp(model_set,"cmip5"))
                    model_gridded_downscaling(obs_src, model, ensemble, scenario, varname, model_set, "prcp_distrib",prcp_distrib, "dry_run",true, "grgn", grgn, "gridsize", gridsize, "mdl_yrs", mdl_yrs, "hist_yrs", hist_yrs, "lonrange",lon_ranges(i,:), "pdf_map_method", pdf_map_method, "use_scratch", use_scratch, "do_continue", do_continue, "do_extended", do_extended, "is_tll", is_tll, "onhost", onhost, "cluster", cluster, "nworkers", nworkers, "nnodes", nnodes, "mb_per_core", mb_per_core, "output_dirs", output_dirs);
                else
                    model_gridded_downscaling(obs_src, model, ensemble, scenario, varname, model_set, "prcp_distrib",prcp_distrib, "dry_run",true, "grgn", grgn, "gridsize", gridsize, "mdl_yrs", mdl_yrs, "hist_yrs", hist_yrs, "lonrange",lon_ranges(i,:), "pdf_map_method", pdf_map_method, "use_scratch", use_scratch, "do_continue", do_continue, "do_extended", do_extended, "is_tll", is_tll, "onhost", onhost, "cluster", cluster, "nworkers", nworkers, "nnodes", nnodes, "mb_per_core", mb_per_core, "output_dirs", output_dirs);
                end
            end
            
        elseif (any(strcmp(onhost, ["neys","icsf-kmac","icsf-lmac"])))
            fprintf(fid,"cd %s\n", sourcecode_dir);
            if (isempty(lon_ranges))
                if (strcmp(model_set, "cmip5"))
                    fprintf(fid,"matlab -nodisplay -r 'cd %s ; retval=model_gridded_downscaling(""%s"",""%s"",""%s"",""%s"",""%s"",""%s"", ""prcp_distrib"",""%s"",                    ""gridsize"", [%d,%d], ""mdl_yrs"", [%d,%d], ""hist_yrs"", [%d,%d], ""pdf_map_method"", ""%s"", ""do_continue"", %d, ""do_extended"", %d, ""is_tll"", %d, ""onhost"", ""%s""); quit(retval); ' \n", ...
                            sourcecode_dir, obs_src, model, ensemble, scenario, varname, model_set, prcp_distrib,       gridsize, mdl_yrs, hist_yrs, pdf_map_method, do_continue, do_extended, is_tll, onhost);
                else
                    fprintf(fid,"matlab -nodisplay -r 'cd %s ; retval=model_gridded_downscaling(""%s"",""%s"",""%s"",""%s"",""%s"",""%s"", ""prcp_distrib"",""%s"", ""grgn"", ""%s"", ""gridsize"", [%d,%d], ""mdl_yrs"", [%d,%d], ""hist_yrs"", [%d,%d], ""pdf_map_method"", ""%s"", ""do_continue"", %d, ""do_extended"", %d, ""is_tll"", %d, ""onhost"", ""%s""); quit(retval); '\n", ...
                            sourcecode_dir, obs_src, model, ensemble, scenario, varname, model_set, prcp_distrib, grgn, gridsize, mdl_yrs, hist_yrs, pdf_map_method, do_continue, do_extended, is_tll, onhost);
                end
            else
                if (strcmp(model_set, "cmip5"))
                    fprintf(fid,"matlab -nodisplay -r 'cd %s ; retval=model_gridded_downscaling(""%s"",""%s"",""%s"",""%s"",""%s"",""%s"", ""%s"", ""prcp_distrib"",""%s"",                   ""gridsize"",[%d,%d], ""mdl_yrs"", [%d,%d], ""hist_yrs"", [%d,%d], ""lonrange"",[%6g,%6g], ""pdf_map_method"", ""%s"", ""do_continue"", %d, ""do_extended"", %d, ""is_tll"", %d, ""onhost"", ""%s""); quit(retval); '\n", ...
                            sourcecode_dir, obs_src, model, ensemble, scenario, varname, model_set, prcp_distrib,       gridsize, mdl_yrs, hist_yrs, lon_ranges(i,:), pdf_map_method, do_continue, do_extended, is_tll, onhost);
                else
                    fprintf(fid,"matlab -nodisplay -r 'cd %s ; retval=model_gridded_downscaling(""%s"",""%s"",""%s"",""%s"",""%s"",""%s"", ""%s"", ""prcp_distrib"",""%s"", ""grgn"", ""%s"", ""gridsize"",[%d,%d], ""mdl_yrs"", [%d,%d], ""hist_yrs"", [%d,%d], ""lonrange"",[%6g,%6g], ""pdf_map_method"", ""%s"", ""do_continue"", %d, ""do_extended"", %d, ""is_tll"", %d, ""onhost"", ""%s""); quit(retval); '\n", ...
                            sourcecode_dir, obs_src, model, ensemble, scenario, varname, model_set, prcp_distrib, grgn, gridsize, mdl_yrs, hist_yrs, lon_ranges(i,:), pdf_map_method, do_continue, do_extended, is_tll, onhost);
                end
            end
            fprintf(fid,"\nretval=$?      # save return value from model_gridded_downscaling\n");
            
            fprintf(fid,"\nif (( $retval > 0 ))\n");
            fprintf(fid,"then\n");
            fprintf(fid,"     source %s\n", next_job_name);
            fprintf(fid,"fi\n");
            
            fclose(fid);
            
            if (isempty(lon_ranges))
                if (strcmp(model_set,"cmip5"))
                    model_gridded_downscaling(obs_src, model, ensemble, scenario, varname, model_set, "prcp_distrib",prcp_distrib, "dry_run",true,               "gridsize", gridsize, "mdl_yrs", mdl_yrs, "hist_yrs", hist_yrs,                             "pdf_map_method", pdf_map_method, "do_continue", do_continue, "do_extended", do_extended, "is_tll", is_tll, "onhost", onhost, "output_dirs", output_dirs);
                els
                    model_gridded_downscaling(obs_src, model, ensemble, scenario, varname, model_set, "prcp_distrib",prcp_distrib,  "dry_run",true, "grgn", grgn, "gridsize", gridsize, "mdl_yrs", mdl_yrs, "hist_yrs", hist_yrs,                             "pdf_map_method", pdf_map_method, "do_continue", do_continue, "do_extended", do_extended, "is_tll", is_tll, "onhost", onhost, "output_dirs", output_dirs);
                end
            else
                if (strcmp(model_set,"cmip5"))
                    model_gridded_downscaling(obs_src, model, ensemble, scenario, varname, model_set, "prcp_distrib",prcp_distrib, "dry_run",true,               "gridsize", gridsize, "mdl_yrs", mdl_yrs, "hist_yrs", hist_yrs, "lonrange",lon_ranges(i,:), "pdf_map_method", pdf_map_method, "do_continue", do_continue, "do_extended", do_extended, "is_tll", is_tll, "onhost", onhost, "output_dirs", output_dirs);
                else
                    model_gridded_downscaling(obs_src, model, ensemble, scenario, varname, model_set, "prcp_distrib",prcp_distrib, "dry_run",true, "grgn", grgn, "gridsize", gridsize, "mdl_yrs", mdl_yrs, "hist_yrs", hist_yrs, "lonrange",lon_ranges(i,:), "pdf_map_method", pdf_map_method, "do_continue", do_continue, "do_extended", do_extended, "is_tll", is_tll, "onhost", onhost, "output_dirs", output_dirs);
                end
            end
            
        end            

    end
    join_name = next_job_name;
    script_names(end+1) = local_next_job_name;

end

%            create_gridded_join_file(runtbl(i,:), region, obs_src, gridsize, mdl_yrs, join_name, model_set, email_addr, do_extended, sourcecode_dir, base_outdir, onhost, local_script_dir);
function     create_gridded_join_file(   tbl,      region, obs_src, gridsize, mdl_yrs, join_name, model_set, email_addr, do_extended, sourcecode_dir, base_outdir, onhost, local_script_dir)

    model    = string(tbl.model(1));
    ensemble = string(tbl.ensemble(1));
    varnames = string(tbl.varname(1));
    scenario = string(tbl.scenario(1));
    if (strcmp(model_set,"cmip6"))
        grgn = string(tbl.grgn(1));
    else
        grgn = [];
    end
        
        % add list of extended-run vars to list of varnames to copy over.
    if (do_extended)
        extended_vars = [sprintf("%s_zvals", tbl.varname{1})];   %#ok<NBRAK> % this will eventually be a list of variables to include.
        for i=1:length(extended_vars)
            varnames(end+1) = extended_vars(i); %#ok<AGROW>
        end
    end
         % put together the list of varnames to copy to the output file.
    if (length(varnames) == 1)
        varnames_str = sprintf(" ""%s"" ", varnames);
    else
        varnames_str = sprintf("[ ""%s"" ", varnames(1));
        for i=2:length(varnames)
            varnames_str = sprintf("%s; ""%s"" ", varnames_str, varnames(i));
        end
        varnames_str = strcat(varnames_str, "]");
    end
         
    [~,job_name,ext] = fileparts(join_name);
    local_join_name = fullfile(local_script_dir, sprintf("%s%s", job_name, ext));
    fid = fopen(local_join_name, "w");
    system(sprintf("chmod +x %s", local_join_name));
       
    outfiles_dir = fullfile(base_outdir, "ARRM_V2_out");            
    slurmout_dir = fullfile(base_outdir, "slurm_out");
    slurmerr_dir = fullfile(base_outdir, "slurm_err");    
    
    if (any(strcmp(onhost, ["nocona","quanah"])))
        partition = onhost;
        if (strcmp(partition,"quanah"))
            queue = "quanah";
        else
            queue = strings(0);
        end
        fprintf(fid,"#!/bin/bash\n");
        fprintf(fid,"#SBATCH --mail-user=%s\n", email_addr);
        fprintf(fid,"#SBATCH --mail-type=BEGIN,END,FAIL\n");
        fprintf(fid,"#SBATCH --job-name=%s\n", job_name);
        fprintf(fid,"#SBATCH --output=%s/%%x.o%%j\n", slurmout_dir);
        fprintf(fid,"#SBATCH --error=%s/err.%%x.%%j\n", slurmerr_dir);
        fprintf(fid,"#SBATCH --partition=%s\n", partition);
        if (~isempty(queue))
            fprintf(fid,"#SBATCH --queue=%s\n", queue);
        end
        fprintf(fid,"#SBATCH --nodes=1\n");
        fprintf(fid,"#SBATCH --ntasks=%d\n", 1);
        fprintf(fid,"#SBATCH --mem=16G\n");
        fprintf(fid,"#SBATCH --chdir=%s\n", sourcecode_dir);
        fprintf(fid,"\ncd %s\n\n", sourcecode_dir);
        fprintf(fid,"echo ""running on node $HOSTNAME""\n\n");
        
        
                % Make sure we have all the modules loaded that we need
%         fprintf(fid,"\n# Load the modules needed.  \n#The load_my_modules script loads the appropriate modules depending on the architecture\n\n");
%         fprintf(fid,"#source /home/iscottfl/bin/load_my_modules\n");
%         fprintf(fid,"#module load matlab\n");
        fprintf(fid, "\nexport PATH=$PATH:/opt/apps/snfs/RedRaider/matlab/R2021b/bin  # for now, to get R2021b, we're setting the path explicitly.  HPCC will create a module for it shortly.\n\n");
        fprintf(fid,"module list\n\n");
                        
        fprintf(fid,"matlab -nodisplay -r 'cd %s ; rejoin_downscaled_minifiles(""%s"", ""%s"", ""%s"", %s, ""%s"", ""%s"", ""%s"", %s, ""%s"", [%d,%d], [%d,%d], false); quit; ' > %s/out.${SLURM_JOB_ID}.${SLURM_JOB_NAME}\n", ...
                    sourcecode_dir, region, model, ensemble, varnames_str, scenario, grgn, model_set, string(true), obs_src, gridsize, mdl_yrs, outfiles_dir);
    elseif (strcmp(onhost, "neys"))
        fprintf(fid,"matlab -nodisplay -r 'cd %s ; rejoin_downscaled_minifiles(""%s"", ""%s"", ""%s"", %s, ""%s"", ""%s"", ""%s"", %s, ""%s"", [%d,%d], [%d,%d], false); quit; ' | tee %s/out.${SLURM_JOB_ID}.${SLURM_JOB_NAME}\n", ...
                    sourcecode_dir, region, model, ensemble, varnames_str, scenario, grgn, model_set, string(true), obs_src, gridsize, mdl_yrs, outfiles_dir);
    end
        
    fclose(fid);
end    

% These routines (for station runs) need updating still.  icsf 2/11/2021
%        [join_name, script_names] = create_station_run_file(runtbl(i,:), gridsize, mdl_yrs, hist_yrs, job_lon_ranges, obs_src, model_set, pdf_map_method, onhost, nworkers, use_scratch, email_addr, script_dir, sourcecode_dir, base_outdir, nstns, do_continue, do_extended, is_tll, use_MPS, nnodes, mb_per_core, unMatched);
function [join_name, script_names] = create_station_run_file(   tbl,      gridsize, mdl_yrs, hist_yrs,     lon_ranges, obs_src, model_set, pdf_map_method, far_outlier_thresh, far_outlier_anchor_pt, sigma_normalize, exclusive, nworkers, use_scratch, email_addr, script_dir, sourcecode_dir, base_outdir, nstns, do_continue, do_extended, is_tll, USE_MPS, nnodes, mb_per_core, onhost, local_script_dir, unMatched) %#ok<INUSL>

    if (~isempty(tbl)), error("station stuff hasn't been debugged and tested yet! \n"); end
    model = string(tbl.model(1));
    ensemble = string(tbl.ensemble(1));
    varname = string(tbl.varname(1));
    scenario = string(tbl.scenario(1));
    if (strcmp(model_set,"cmip6"))
        grgn = string(tbl.grgn(1));
    else
        grgn = [];
    end
  
    if (~isempty(unMatched)), error("can't handle unMatched yet"); end
    
    if (isempty(lon_ranges))
        njobs = 1;
    else
        njobs = size(lon_ranges, 1);
    end
        
    if (strcmp(onhost,"quanah"))  % we're on quanah nodes
        partition="quanah";
        run_on_hpcc = true;
        queue="quanah";
    elseif (strcmp(onhost,"nocona"))
        partition="nocona";       % we're on nocona nodes.  
        run_on_hpcc = true;
        queue=strings(0);
    else
        partition = "";
        run_on_hpcc = false;
        queue = strings(0);
    end
    
    outfiles_dir = fullfile(base_outdir, "ARRM_V2_out");            
    slurmout_dir = fullfile(base_outdir, "slurm_out");
    slurmerr_dir = fullfile(base_outdir, "slurm_err");

    script_names = strings(0,0);
    for i=1:njobs
                        
        if (strlength(grgn)==0)
            run_lbl    = sprintf("stn_%s_%s_%s_%s_%s_%s",    model_set, obs_src, model, ensemble, varname, scenario);
        else
            run_lbl    = sprintf("stn_%s_%s_%s_%s_%s_%s_%s", model_set, obs_src, model, ensemble, varname, scenario, grgn);
        end
                        
        run_name       = fullfile(      script_dir, sprintf("down.%s.%02d.sh",   run_lbl, i));
        local_run_name = fullfile(local_script_dir, sprintf("down.%s.%02d.sh",   run_lbl, i));
        if (i==njobs)
            next_job_name  = fullfile(script_dir, sprintf("rejoin.%s.sh", run_lbl)); 
        else
            next_job_name  = fullfile(script_dir, sprintf("down.%s.%02d.sh",   run_lbl, i+1));
        end
        
        [~,job_name,~] = fileparts(run_name);
        
        fid = fopen(local_run_name, "w");
        system(sprintf("chmod +x %s", run_name));

        if (run_on_hpcc)
            fprintf(fid,"#!/bin/bash\n");
            fprintf(fid,"#SBATCH --job-name=%s\n", job_name);
            fprintf(fid,"#SBATCH --mail-user=%s\n", email_addr);
            fprintf(fid,"#SBATCH --mail-type=BEGIN,END,FAIL\n");
            fprintf(fid,"#SBATCH --output=%s/%%x.o%%j\n", slurmout_dir);
            fprintf(fid,"#SBATCH --error=%s/err.%%x.%%j\n", slurmerr_dir);
            fprintf(fid,"#SBATCH --partition=%s\n", partition);
            if (~isempty(queue))
                fprintf(fid,"#SBATCH --queue=%s\n", queue);
            end
            if (use_MPS)
                fprintf(fid,"#SBATCH --ntasks=%d\n", 1);
            else
                fprintf(fid,"#SBATCH --nodes=%d", nnodes);
                if (exclusive)
                    fprintf(fid,"#SBATCH --exclusive\n");
                    fprintf(fid, "#SBATCH --ntasks=%d\n", nworkers+1);      % add 1 for the parent process.
                else
                    fprintf(fid, "#SBATCH --mem-per-cpu=%d\n", mb_per_core);
                    fprintf(fid,"#SBATCH --ntasks=%d\n", nworkers+1);       % add 1 for the parent process.
                end
            end
            fprintf(fid,"#SBATCH --chdir=%s\n", sourcecode_dir);
            
            fprintf(fid,"\ncd %s\n\n", sourcecode_dir);
            
                % Make sure we have all the modules loaded that we need
%             fprintf(fid,"\n# Load the modules needed.  \n#The load_my_modules script loads the appropriate modules depending on the architecture\n\n");
%             fprintf(fid,"source /home/iscottfl/bin/load_my_modules\n");
%             fprintf(fid,"#module load matlab\n");
            fprintf(fid, "\nexport PATH=$PATH:/opt/apps/snfs/RedRaider/matlab/R2021b/bin  # for now, to get R2021b, we're setting the path explicitly.  HPCC will create a module for it shortly.\n\n");
            fprintf(fid,"module list\n\n");
            fprintf(fid,"echo ""running on node $HOSTNAME""\n\n");
                        
                % create the output folders if they don't exist
%             fprintf(fid,"\n# create the output folders if they don't exist\n");
%             fprintf(fid,"if [[ ! -d %s ]]; then mkdir %s; fi\n", base_dir, base_dir);
%             fprintf(fid,"if [[ ! -d %s ]]; then mkdir %s; fi\n", outfiles_dir, outfiles_dir);
%             fprintf(fid,"if [[ ! -d %s ]]; then mkdir %s; fi\n", slurmout_dir, slurmout_dir);
%             fprintf(fid,"if [[ ! -d %s ]]; then mkdir %s; fi\n", slurmerr_dir, slurmerr_dir);
%             fprintf(fid,"#\n\n");

            if (isempty(lon_ranges))
                if (strcmp(model_set, "cmip5"))
                    fprintf(fid,"matlab -nodisplay -r 'cd /lustre/work/iscottfl/ARRM_V2 ; model_gridded_downscaling(""%s"",""%s"",""%s"",""%s"",""%s"",""%s"",                    ""gridsize"", [%d,%d], ""mdl_yrs"", [%d,%d], ""hist_yrs"", [%d,%d], ""pdf_map_method"", ""%s"", ""far_outlier_thresh"", %f, ""far_outlier_anchor_pt"", %f, ""sigma_normalize"", %f, ""use_scratch"",%d, ""do_continue"", %d, ""do_extended"", %d, ""is_tll"", %d, ""onhost"", ""%s""); quit; ' > %s/out.${SLURM_JOB_ID}.${SLURM_JOB_NAME}\n", ...
                            obs_src, model, ensemble, scenario, varname, model_set,       gridsize, mdl_yrs, hist_yrs, pdf_map_method, far_outlier_thresh, far_outlier_anchor_pt, sigma_normalize, use_scratch, do_continue, do_extended, is_tll, onhost, outfiles_dir);
                else
                    fprintf(fid,"matlab -nodisplay -r 'cd /lustre/work/iscottfl/ARRM_V2 ; model_gridded_downscaling(""%s"",""%s"",""%s"",""%s"",""%s"",""%s"", ""grgn"", ""%s"", ""gridsize"", [%d,%d], ""mdl_yrs"", [%d,%d], ""hist_yrs"", [%d,%d], ""pdf_map_method"", ""far_outlier_thresh"", %f, ""far_outlier_anchor_pt"", %f, ""sigma_normalize"", %f, ""%s"", ""use_scratch"",%d, ""do_continue"", %d, ""do_extended"", %d, ""is_tll"", %d, ""onhost"", ""%s""); quit; ' > %s/out.${SLURM_JOB_ID}.${SLURM_JOB_NAME}\n", ...
                            obs_src, model, ensemble, scenario, varname, model_set, grgn, gridsize, mdl_yrs, hist_yrs, pdf_map_method, far_outlier_thresh, far_outlier_anchor_pt, sigma_normalize, use_scratch, do_continue, do_extended, is_tll, onhost, outfiles_dir);
                end
            else
                if (strcmp(model_set, "cmip5"))
                    fprintf(fid,"matlab -nodisplay -r 'cd /lustre/work/iscottfl/ARRM_V2 ; model_gridded_downscaling(""%s"",""%s"",""%s"",""%s"",""%s"",""%s"", ""%s"",                    ""gridsize"",[%d,%d], ""mdl_yrs"", [%d,%d], ""hist_yrs"", [%d,%d], ""lonrange"",[%6g,%6g], ""pdf_map_method"", ""%s"", ""far_outlier_thresh"", %f, ""far_outlier_anchor_pt"", %f, ""sigma_normalize"", %f, ""use_scratch"",%d, ""do_continue"", %d, ""do_extended"", %d, ""is_tll"", %d, ""onhost"", ""%s""); quit; ' > %s/out.${SLURM_JOB_ID}.${SLURM_JOB_NAME}\n", ...
                            obs_src, model, ensemble, scenario, varname, model_set,       gridsize, mdl_yrs, hist_yrs, lon_ranges(i,:), pdf_map_method, far_outlier_thresh, far_outlier_anchor_pt, sigma_normalize, use_scratch, do_continue, do_extended, is_tll, onhost, outfiles_dir);
                else
                    fprintf(fid,"matlab -nodisplay -r 'cd /lustre/work/iscottfl/ARRM_V2 ; model_gridded_downscaling(""%s"",""%s"",""%s"",""%s"",""%s"",""%s"", ""%s"", ""grgn"", ""%s"", ""gridsize"",[%d,%d],""mdl_yrs"", [%d,%d], ""hist_yrs"", [%d,%d], ""lonrange"",[%6g,%6g], ""pdf_map_method"", ""far_outlier_thresh"", %f, ""far_outlier_anchor_pt"", %f, ""sigma_normalize"", %f, ""%s"", ""use_scratch"",%d, ""do_continue"", %d, ""do_extended"", %d, ""is_tll"", %d, ""onhost"", ""%s""); quit; ' > %s/out.${SLURM_JOB_ID}.${SLURM_JOB_NAME}\n", ...
                            obs_src, model, ensemble, scenario, varname, model_set, grgn, gridsize, mdl_yrs, hist_yrs, lon_ranges(i,:), pdf_map_method, far_outlier_thresh, far_outlier_anchor_pt, sigma_normalize, use_scratch, do_continue, do_extended, is_tll, onhost, outfiles_dir);
                end
            end
            fprintf(fid,"sbatch %s\n", next_job_name);
            
            script_names(end+1) = run_name; %#ok<AGROW>
        elseif (strcmp(onhost, "neys") || strcmp(onhost,"icsf-kmac"))
            fprintf(fid,"cd /Users/iscottfl/Desktop/atmos2/ARRM_V2/\n");
            if (isempty(lon_ranges))
                if (strcmp(model_set, "cmip5"))
                    fprintf(fid,"matlab -nodisplay -r 'cd /Users/iscottfl/Desktop/atmos2/ARRM_V2/ ; model_gridded_downscaling(""%s"",""%s"",""%s"",""%s"",""%s"",""%s"",                    ""gridsize"", [%d,%d], ""mdl_yrs"", [%d,%d], ""hist_yrs"", [%d,%d], ""pdf_map_method"", ""%s"", ""far_outlier_thresh"", %f, ""far_outlier_anchor_pt"", %f, ""sigma_normalize"", %f, ""do_continue"", %d, ""do_extended"", %d, ""is_tll"", %d, ""onhost"", ""%s""); quit; ' \n", ...
                            obs_src, model, ensemble, scenario, varname, model_set,       gridsize, mdl_yrs, hist_yrs, pdf_map_method, far_outlier_thresh, far_outlier_anchor_pt, sigma_normalize, do_continue, do_extended, is_tll, onhost);
                else
                    fprintf(fid,"matlab -nodisplay -r 'cd /Users/iscottfl/Desktop/atmos2/ARRM_V2/ ; model_gridded_downscaling(""%s"",""%s"",""%s"",""%s"",""%s"",""%s"", ""grgn"", ""%s"", ""gridsize"", [%d,%d], ""mdl_yrs"", [%d,%d], ""hist_yrs"", [%d,%d], ""pdf_map_method"", ""%s"", ""far_outlier_thresh"", %f, ""far_outlier_anchor_pt"", %f, ""sigma_normalize"", %f, ""do_continue"", %d, ""do_extended"", %d, ""is_tll"", %d, ""onhost"", ""%s""); quit; '\n", ...
                            obs_src, model, ensemble, scenario, varname, model_set, grgn, gridsize, mdl_yrs, hist_yrs, pdf_map_method, far_outlier_thresh, far_outlier_anchor_pt, sigma_normalize, do_continue, do_extended, is_tll, onhost);
                end
            else
                if (strcmp(model_set, "cmip5"))
                    fprintf(fid,"matlab -nodisplay -r 'cd /Users/iscottfl/Desktop/atmos2/ARRM_V2/ ; model_gridded_downscaling(""%s"",""%s"",""%s"",""%s"",""%s"",""%s"", ""%s"",                    ""gridsize"",[%d,%d], ""mdl_yrs"", [%d,%d], ""hist_yrs"", [%d,%d], ""lonrange"",[%6g,%6g], ""pdf_map_method"", ""%s"", ""far_outlier_thresh"", %f, ""far_outlier_anchor_pt"", %f, ""sigma_normalize"", %f, ""do_continue"", %d, ""do_extended"", %d, ""is_tll"", %d, ""onhost"", ""%s""); quit; '\n", ...
                            obs_src, model, ensemble, scenario, varname, model_set,       gridsize, mdl_yrs, hist_yrs, lon_ranges(i,:), pdf_map_method, far_outlier_thresh, far_outlier_anchor_pt, sigma_normalize, do_continue, do_extended, is_tll, onhost);
                else
                    fprintf(fid,"matlab -nodisplay -r 'cd /lustre/work/iscottfl/ARRM_V2 ;           model_gridded_downscaling(""%s"",""%s"",""%s"",""%s"",""%s"",""%s"", ""%s"", ""grgn"", ""%s"", ""gridsize"",[%d,%d], ""mdl_yrs"", [%d,%d], ""hist_yrs"", [%d,%d], ""lonrange"",[%6g,%6g], ""pdf_map_method"", ""%s"", ""far_outlier_thresh"", %f, ""far_outlier_anchor_pt"", %f, ""sigma_normalize"", %f, ""do_continue"", %d, ""do_extended"", %d, ""is_tll"", %d, ""onhost"", ""%s""); quit; '\n", ...
                            obs_src, model, ensemble, scenario, varname, model_set, grgn, gridsize, mdl_yrs, hist_yrs, lon_ranges(i,:), pdf_map_method, far_outlier_thresh, far_outlier_anchor_pt, sigma_normalize, do_continue, do_extended, is_tll, onhost);
                end
            end
            fprintf(fid,"./%s\n", next_job_name);
            
            script_names(end+1) = run_name; %#ok<AGROW>
        end            
        fclose(fid);

    end
    join_name = next_job_name;
    script_names(end+1) = join_name;

    if (isempty(lon_ranges))
        model_gridded_downscaling(obs_src, model, ensemble, scenario, varname, model_set, "grgn", grgn, "gridsize", gridsize, "mdl_yrs", mdl_yrs, "hist_yrs", hist_yrs, "dry_run",onhost, "pdf_map_method", pdf_map_method,                             "use_scratch", use_scratch, "do_continue", do_continue, "do_extended", do_extended, "is_tll", is_tll, "onhost", onhost);
    else
        model_gridded_downscaling(obs_src, model, ensemble, scenario, varname, model_set, "grgn", grgn, "gridsize", gridsize, "mdl_yrs", mdl_yrs, "hist_yrs", hist_yrs, "dry_run",onhost, "pdf_map_method", pdf_map_method, "lonrange",lon_ranges(i,:), "use_scratch", use_scratch, "do_continue", do_continue, "do_extended", do_extended, "is_tll", is_tll, "onhost", onhost);
    end
end

function create_station_join_file(tbl,         region, obs_src, gridsize, mdl_yrs, join_name, model_set, prcp_distrib, email_addr, do_extended, sourcecode_dir, base_outdir, onhost, local_script_dir) %#ok<INUSD,INUSL>

    error("this function needs updating.  do_extended %d, sourcecode_dir %s, base_outdir %s\n", do_extended, sourcecode_dir, base_outdir);
    model    = string(tbl.model(1)); %#ok<UNRCH>
    ensemble = string(tbl.ensemble(1));
    varname  = string(tbl.varname(1));
    scenario = string(tbl.scenario(1));
    if (strcmp(model_set,"cmip6"))
        grgn = string(tbl.grgn(1));
    else
        grgn = [];
    end
  
    
    if (strlength(grgn)==0)
        run_lbl   = sprintf("%s_%s_%s_%s_%s", obs_src, model, ensemble, varname, scenario);
    else
        run_lbl   = sprintf("%s_%s_%s_%s_%s_%s", obs_src, model, ensemble, varname, scenario, grgn);
    end
    
    outfiles_dir = fullfile(base_outdir, "ARRM_V2_out");            
    slurmout_dir = fullfile(base_outdir, "slurm_out");
    slurmerr_dir = fullfile(base_outdir, "slurm_err");

    fid = fopen(join_name, "w");
    system(sprintf("chmod +x %s", join_name));
    
            % this still needs to be changed to SLURM code, Ian!
    if (strcmp(onhost, "hpcc")) 
        fprintf(fid,"#!/bin/sh\n");
        fprintf(fid,"#$ -V\n");
        fprintf(fid,"#$ -cwd\n");
        fprintf(fid,"#$ -S /bin/bash\n");
        fprintf(fid,sprintf("#$ -M %s\n", email_addr));
        fprintf(fid,"#$ -m be\n");    
        fprintf(fid,"#$ -N join_%s\n", run_lbl);
        fprintf(fid,"#$ -o %s/%%x.o%%j\n", slurmout_dir);
        fprintf(fid,"#$ -e %s/err.%%x.%%j\n", slurmerr_dir);
        fprintf(fid,"#$ -q omni\n");
        fprintf(fid,"#$ -pe sm 1\n");
        fprintf(fid,"#$ -l h_vmem=16G\n");
        fprintf(fid,"#$ -P quanah\n");
        fprintf(fid,"cd /lustre/work/iscottfl/ARRM_V2\n");
        fprintf(fid,"#module load matlab\n");
        fprintf(fid, "\nexport PATH=$PATH:/opt/apps/snfs/RedRaider/matlab/R2021b/bin  # for now, to get R2021b, we're setting the path explicitly.  HPCC will create a module for it shortly.\n\n");
        fprintf(fid,"matlab -nodisplay -r 'cd /lustre/work/iscottfl/ARRM_V2 ;          rejoin_downscaled_minifiles(""%s"", ""%s"", ""%s"", ""%s"", ""%s"", ""%s"", ""%s"", true, ""%s"", [%d,%d], ""single"", [%d,%d]); quit; ' > %s/out.${SLURM_JOB_ID}.${SLURM_JOB_NAME}\n",region, model, ensemble, varname, scenario, grgn, model_set, obs_src, gridsize, mdl_yrs, outfiles_dir);
    elseif (strcmp(onhost, "neys"))
        fprintf(fid,"cd /Users/iscottfl/Desktop/atmos2/ARRM_V2\n");
        fprintf(fid,"matlab -nodisplay -r 'cd /Users/iscottfl/Desktop/atmos2/ARRM_V2 ; rejoin_downscaled_minifiles(""%s"", ""%s"", ""%s"", ""%s"", ""%s"", ""%s"", ""%s"", true, ""%s"", [%d,%d], ""single"", [%d,%d]); quit; '\n",                                                           region, model, ensemble, varname, scenario, grgn, model_set, obs_src, gridsize, mdl_yrs);
    end
        
    fclose(fid);
end    
%
function     [obs_src, model, ensemble, varname, scenario, grgn, model_set, prcp_distrib, region, gridsize, mdl_yrs, hist_yrs, do_clean, job_lon_ranges, exclusive, nworkers, use_scratch, sourcecode_dir, base_outdir, nstns, do_continue, do_extended, is_tll, pdf_map_method, far_outlier_thresh, far_outlier_anchor_pt, sigma_normalize, use_MPS, nnodes, mb_per_core, onhost, email_addr, script_dir, local_script_dir, unMatched] = init_params(obs_src, model, ensemble, varname, scenario, varargin)

% this is all a bit kludged/hardwired to get params from command line in specific order, but using matlab's parser to
% handle missing data.

    p = inputParser;
    p.StructExpand = true;
    p.KeepUnmatched = true;
    
    addParameter(p,"obs_src","",                @(s) any(strcmp(s,["livneh","nrcan","nclimgrid","stations","sheffield", "stations_25","stations_100","stations_1000", "stations_2400"])));
    addParameter(p,"model",strings(0),          @(s) isstring(s));
    addParameter(p,"ensemble",strings(0),       @(s) isstring(s));
    addParameter(p,"varname","",                @(s) isstring(s));
    addParameter(p,"scenario",                  @(s) isstring(s));
    addParameter(p,"grgn","",                   @(s) isstring(s));
    addParameter(p,"model_set","cmip6",         @(s) (isstring(s)) && any(strcmp(s,["","cmip5","cmip6"])));
    addParameter(p,"region",strings(0),         @(s) isstring(s));
    addParameter(p,"gridsize",[  2,   2],       @(s) isnumeric(s) && length(s)==2);
    addParameter(p,"mdl_yrs",[1950,2100],       @(s) isnumeric(s) && length(s)==2);
    addParameter(p,"hist_yrs",[1950,2000],      @(s) isnumeric(s) && length(s)==2);
%   addParameter(p,"job_lon_ranges",[],         @(s) isnumeric(s) && (isempty(s) || size(s,2)==2));
    addParameter(p,"prcp_distrib",strings(0),   @(s) isstring(s) && (isempty(s) || any(strcmp(s,["generalizedpareto", "inversegaussian", "loglogistic", "lognormal","pwr","log"]))));
    addParameter(p,"user",strings(0),           @(s) isstring(s));
    addParameter(p,"base_dir",strings(0),       @(s) isstring(s));
    addParameter(p,"script_dir",strings(0),     @(s) isstring(s));
    addParameter(p,"local_script_dir", strings(0), @(s) isstring(s));
    addParameter(p,"email_addr",strings(0),     @(s) isempty(s) || isstring(s));
    addParameter(p,"do_clean",false,            @(s) islogical(s));
    addParameter(p,"use_dev_folder",false,      @(s) islogical(s) || isstring(s) || ischar(s));
    addParameter(p,"exclusive",true,            @(s) islogical(s)  || any(s==[0,1]));
    addParameter(p,"nworkers",96,               @(s) s>0 && s<=128);
    addParameter(p,"use_scratch",false,         @(s) islogical(s));
    addParameter(p,"do_continue",true,          @(s) islogical(s));
    addParameter(p,"do_extended",false,         @(s) islogical(s));
    addParameter(p,"is_tll",[],                 @(s) isempty(s) || islogical(s) || any(s==[0,1]));
    addParameter(p,"pdf_map_method",strings(0), @(s) isstring(s));
    addParameter(p,"far_outlier_thresh",[],     @(s) isempty(s) || isnumeric(s));
    addParameter(p,"far_outlier_anchor_pt",[],  @(s) isempty(s) || isnumeric(s));
    addParameter(p,"sigma_normalize",[],        @(s) isempty(s) || islogical(s));
    addParameter(p,"use_MPS",true,              @(s) islogical(s));
    addParameter(p,"nnodes",1,                  @(s) isnumeric(s) && s < 8);
    addParameter(p,"mb_per_core",5200*96/97,    @(s) isempty(s) || (s>=0 && s<=515565));    % applies only to HPCC runs.  default will allow 48 workers in 250GB or 96 workers in 500G
    addParameter(p,"onhost", "nocona",          @(s) isstring(s) && any(strcmp(s,["nocona","hpcc","quanah","icsf-kmac","icsf-lmac","laptop"])));
    
    parse(p,"obs_src", obs_src, "model",model, "ensemble",ensemble,"varname",varname,"scenario",scenario,varargin{:});
    obs_src         = p.Results.obs_src;
    model           = string(p.Results.model);
    ensemble        = string(p.Results.ensemble);
    varname         = string(p.Results.varname);
    scenario        = string(p.Results.scenario);
    grgn            = string(p.Results.grgn);
    gridsize        = p.Results.gridsize;
    mdl_yrs         = p.Results.mdl_yrs;
    hist_yrs        = p.Results.hist_yrs;
    prcp_distrib    = p.Results.prcp_distrib;
    do_clean        = p.Results.do_clean;
    region          = p.Results.region;
%   job_lon_ranges  = p.Results.job_lon_ranges;
    user            = string(p.Results.user);
    email_addr      = p.Results.email_addr;
    model_set       = p.Results.model_set;
    exclusive       = p.Results.exclusive;
    nworkers        = p.Results.nworkers;
    mb_per_core     = floor(p.Results.mb_per_core);
    use_scratch     = p.Results.use_scratch;
    base_outdir     = p.Results.base_dir;
    do_continue     = p.Results.do_continue;
    do_extended     = p.Results.do_extended;
    is_tll          = p.Results.is_tll;
    pdf_map_method  = p.Results.pdf_map_method;
    far_outlier_thresh =  p.Results.far_outlier_thresh;
    far_outlier_anchor_pt =  p.Results.far_outlier_anchor_pt;
    sigma_normalize = p.Results.sigma_normalize;
    use_dev_folder  = p.Results.use_dev_folder;     % for dev code testing.  will select /lustre/work/iscottfl/ARRM_V2_dev/ARRM_V2 as sourcecode_dir
    use_MPS         = p.Results.use_MPS;
    nnodes          = p.Results.nnodes;
    onhost          = p.Results.onhost;
    script_dir      = p.Results.script_dir;
    local_script_dir= p.Results.local_script_dir;
    
    unMatched = [];         % figure out how to pass on unmatched, ian!
    
    if (any(strcmp(onhost,["nocona","hpcc","quanah"])))
        run_on_hpcc = true;
        if (strcmp(onhost,"hpcc"))
            onhost = "nocona";
        end
    else
        run_on_hpcc = false;
    end
        
    thishost = get_hostname(true);
    if (strncmp(thishost,"cpu-",4) || strncmp(thishost,"login-",6)), thishost = "nocona"; end
    if (strncmp(  onhost,"cpu-",4) || strncmp(  onhost,"login-",6)),   onhost = "nocona"; end

    if (strcmp(onhost,"quanah") && nworkers > 36)
        nworkers = 36;
    elseif (strcmp(onhost,"neys"))
        nworkers = 18;
    elseif (strcmp(onhost,"icsf-kmac"))
        nworkers = 8;
    elseif (strcmp(onhost, "icsf-lmac"))
        nworkers = 6;
    end
    if (strcmp(onhost,"nocona") && nworkers <= 48 && exclusive) % setting nworkers to 48 or less will allow 2 downscaling jobs to run on the same machine.
        exclusive = false;
        fprintf(2, "NOTE:  exclusive changed (from default of true) to false because nworkers <= 48\n");
    end
    
    if (~isempty(mb_per_core) && mb_per_core ~= 0)
        mb_min = ceil(2000/(nworkers/nnodes));
        mb_max = floor(515565/((nworkers+1)/nnodes));
        if (mb_per_core < mb_min || mb_per_core > mb_max)
            error("error:  mb_per_core (%d) should be in Mbytes/core, min %d  max %d", mb_per_core, mb_min, mb_max);
        end
    end
        
            % skip setting default mem for now.
%     if (on.hpcc_system && isempty(mb_per_core))
%     mworkers = nworkers + ~use_MPS;     % calculate memory usage including 1 for main routine if using local cluster.
%         if (mworkers < 24)
%             mb_per_node = 250000;
%         elseif ( nworkers <=64)
%             mb_per_node = 250000/mworkers;
%         else
%             mb_per_node = 500000/mworkers;
%         end
%     end
    
            % if user specified the sourcecode_dir, use it.
    if (isstring(use_dev_folder))
        sourcecode_dir = string(use_dev_folder);
    elseif (use_dev_folder)
        sourcecode_dir = "/lustre/work/iscottfl/ARRM_V2_dev/ARRM_V2";
    elseif (~run_on_hpcc)
            sourcecode_dir = "/Users/iscottfl/Desktop/atmos2/ARRM_V2";
    else
        % we'll assume we should be running ARRM_V2 from the current
        % working directory, as long as it ends in .../ARRM_V2
        % This will need changing to run elsewhere, Ian!

        sourcecode_dir = "/lustre/work/iscottfl/ARRM_V2";
    end
    mydir=strsplit(string(sourcecode_dir),filesep);
    if (~strcmp(mydir(end),"ARRM_V2"))
        error("error:  setup_slurm_run(...) must be run from an ARRM_V2 source code folder");
    end        

        % make sure we're not mixing variables of different types, as their
        % run parameters differ.  For example, precip should be run with
        % "pdf_map_method"="clim", while temperature should be run with
        % "pdf_map_method"="linear".
        %
        %   Also, set the pdf_map_method appropriately if it wasn't set by the user.
    isTempRun = [];
    for i=1:length(varname)
        if (any(strcmp(varname(i),["pr","PR","precip","prcp","precipitation"])))
            if (length(varname)>1)
                error("error:  Precip must be run separately, as run parameters are different");
            end
            if (isempty(pdf_map_method)), pdf_map_method = "clim"; end
            if (isempty(sigma_normalize)), sigma_normalize = false; end
            if(isempty(far_outlier_thresh)), far_outlier_thresh = .0228; end
            if (isempty(far_outlier_anchor_pt)), far_outlier_anchor_pt =.1587; end
        elseif (any(strcmpi(varname(i),["tmax","tmin","tas","tavg","tasmax","tasmin"])))
            if (~isempty(isTempRun) && ~isTempRun)
                error("error:  Temperature variables cannot be run with non-temperature variables, as run parameters are different");
            else
                isTempRun = true;
                if (isempty(pdf_map_method)), pdf_map_method = "linear"; end
                if(isempty(far_outlier_thresh)), far_outlier_thresh = 1e-4; end
                if (isempty(far_outlier_anchor_pt)), far_outlier_anchor_pt =.0228; end
                if (isempty(sigma_normalize)), sigma_normalize = true; end
            end
        else
            if (~isempty(isTempRun) && isTempRun)
                error("error:  Temperature variables cannot be run with non-temperature variables, as run parameters are different");
            else
                isTempRun = false;
                if (isempty(pdf_map_method)), pdf_map_method = "linear"; end
            end
        end
    end
            
    obs_src = lower(obs_src);
    nstns = 0;
    if (strcmp(obs_src, "sheffield"))
        if (isempty(region)), region = "global"; end
        job_lon_ranges = [];                % do entire sheffield area in one job.  MIGHT NEED TO CHANGE THIS!
    elseif (strcmp(obs_src, "livneh"))
        if (isempty(region)), region = "conus_livneh"; end
        job_lon_ranges = [];                % do entire livneh area in one job.
    elseif (strcmp(obs_src, "nclimgrid"))
        if (isempty(region)), region = "nclimgrid"; end
        job_lon_ranges = [];                % do entire nclimgrid area in one job.
    elseif (strcmp(obs_src, "nrcan"))
        if (isempty(region)), region = "nrcan"; end
        job_lon_ranges = [];                % do entire nrcan area in one job.
%       job_lon_ranges = [-141.01,-96.49; -96.49, -51.99]; % split area into two regions so processing doesn't run beyond 48-hour limit.
    elseif (strcmp(obs_src, "stations"))
        nstns = "all";
    elseif (strcmp(obs_src, "stations_25"))
        nstns = 25;        
    elseif (strcmp(obs_src, "stations_100"))
        nstns = 100;        
    elseif (strcmp(obs_src, "stations_1000"))
        nstns = 1000;        
    elseif (strcmp(obs_src, "stations_2400"))
        nstns = 2400;        
    else
        error("error:  bad obs_src:  %s\nmust be ""sheffield"", ""livneh"", ""nclimgrid"" or ""nrcan""\n", obs_src);
    end
    
            % find the model set and make sure model, scenario, ensemble (and grgn) are all for the same model_set
    if (strlength(model_set)==0)
        model_set = find_model_set([to_row(model); to_row(ensemble);to_row(scenario); to_row(grgn)]);
        if (length(model_set) ~= 1)
            error("error:   cannot determine model_set from input");
        end
    end    
    if (isempty(is_tll))
        if (strcmp(model_set,"cmip5"))
            is_tll = false;
        else
            is_tll = true;
        end
    end

    if (is_prcp_variable(varname))
        if (isempty(prcp_distrib) || strlength(prcp_distrib) == 0)
%           prcp_distrib = "loglogistic";   % default:  use loglogistic mapping to gaussian for precip.
            prcp_distrib = "pwr";           % default:  use pwr scaling for precip.
        end
    else
        prcp_distrib = strings(0);  % make sure it's empty if not doing a precip run.
    end
    
    
    if (~isempty(fieldnames(p.Unmatched))) 
        fprintf(2, "error: unexpected input parameters: ");
        disp(p.Unmatched);
        error("unexpected input parameters");
    end         
     
    if (isempty(user))
        user = getusername();
    end
 
    if (any(strcmp(onhost,["nocona","quanah"])))
        script_base = "/home/iscottfl";
    else
        script_base = ".";
    end

    if (isempty(base_outdir) || strlength(base_outdir)==0)
        if (any(strcmp(onhost,["nocona","quanah"])))
            base_outdir = fullfile("/lustre/scratch",user);
        else
            base_outdir = fullfile("/home",user);
        end
    end
                
    
    if (isempty(script_dir))
        script_dir = fullfile(script_base,sprintf("%s_scripts", onhost), user);
    elseif (strcmp(script_dir,"."))
        script_dir = fullfile(sourcecode_dir,  sprintf("%s_scripts", onhost), user);
    end
    
    if (isempty(local_script_dir))
        if (strcmp(onhost,thishost))
            local_script_dir = script_dir;
        else
            local_script_dir = sprintf("./%s_scripts/%s",onhost, user);
        end
    end
    
    if (~isfolder(local_script_dir))
        error("error:  local_script_dir %s does not exist\n", local_script_dir);
    end
                    
        % find the email address to use from user's login:
    if (isempty(email_addr))        
        if (any(strcmp(user,["iscottfl","icsf"])))
            email_addr = "ian.scott-fleming@ttu.edu";
        elseif (strcmp(user, "astoner"))
            email_addr = "anne.stoner@ttu.edu";
        elseif (strcmp(user, "khayhoe"))
            email_addr = "katharine.hayhoe@ttu.edu";
        elseif (strcmp(user, "hameibra"))
            email_addr = "hibrahim@illinois.edu";
        else
            error("error:  please specify an email_addr for user %s", user);
        end
    end
end
