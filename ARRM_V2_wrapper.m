function [outnames, retval] = ARRM_V2_wrapper( model,ensemble,varname, scenario, base_yrs, hist_yrs, mdl_yrs, varargin)
%        [outnames, retval] = ARRM_V2_wrapper( model,ensemble,varname, scenario, base_yrs, hist_yrs, mdl_yrs, varargin)
%
%   Program to run ARRM_V2_run on multiple locations, optionally using
%   parfor to run on multiple cores
%
%   Outputs:
%       cell array of output filenames (if not running in parallel)
%       empty cell array (if running in parallel)
%   Inputs:
%
%                           First 7 should be slef-explanatory
%       model, ensemble, varname, scenario      specs to select the model data to use.
%                                                   varname must be a single string value.
%                                                   model, ensemble & scenario can be arrays of values, which
%                                                   will run multiple sets.
%                                                   
%       base_yrs, hist_yrs, mdl_yrs     date ranges to run, given as  [start_yr, end_yr], as years only
%
%
%       varargin            Name/Value pairs for special parameters. 
%                           for station run:  
%                               region          "ncamerica", "samerica", "asia", "africa", "europe", "spacific"
%                           for gridded run:    
%                               nothing
%
%                   ----Required varargin pairs:
%
%       "gridded_obsname",...  For Gridded downscaling
%                               If running more than one variable, specify one gridded obsname for each variable
%                                   [obs_netcdf_var1_filename, obs_netcdf_var2_filename] 
%               OR                          
%       "siteTable", ...    siteTable (read via QC_get_site_table(...) ), or name of station file to read site table from.
%                               If it's a filename, must be the name of a QC'd netcdf station file.
%                               If not specified, then netcdf filename will be assembled from variable, ensemble,
%                               model and scenario, based on the filename prototype in the local DataParams.  
%
%       "testsites", ...    string array of test site names corresponding to latrange & lonrange arrays
%                           for testing gridded code at specific test locations instead of gridded area.
%                           Output is 1 matfile for each site
%                       For limiting the lat/lons and/or stations included in run
%                       lats & lons or latrange & lonrange are used to identify the minifile gridboxes that are run.
%                       lats & lons specify which lats & lons get run, while latrange and lonrange specify only the
%                       range of lats and lons, and the code uses the lat & lon locations from the observation data.
%       "lats".[...]        list of output lats & lons.  For stations, will find closest station to each lat/lon pair.
%       "lons",[...]        For gridded, lats & lons specifies the range of obs locations which will be downscaled. 
%                           if testsites specified, will find closest lat/lon in obs file, and interpolate model data
%                               for testsites, lat & lon are paired.
%
%               AND/OR
%       "stations",[...]    For Station downscaling.  
%                               list of sites (or "all" or region name ("ncamerica","caribbean","samerica", etc.)
%                               list of sites can be "all", and then limited by setting "lats" and "lons" to min & max
%                               lat and lon range.
%
%                   ----Optional varargin pairs:
%
%       "do_parfor",t_or_f      true/false.  If true,  loops w/ parfor instead of for (i.e., uses multiple processors on current system)
%                               If false, loops using a regular for loop.  Slower, but helpful for debugging.
%                               Generally, use true for production runs
%
%       "trend_yrs", trend_yrs              Max range of years to use for calculating trend.  
%                                               (Will be limited to actual data years)
%                                           For temperature, defaults to: [base_yrs(1), mdl_yrs(end)]
%                                           For precip,      defaults to: [];
%
%       "nworkers",val    vals is max # of child processes to start running, if Do_Parfor is true. 
%                               default:  0 (which means:  use max numper possible on host)
%                               Leave this parameter off unless you want to limit the number of child processes.
%                               You would do that if you want to be able to run other heavy tasks on your system.
%       "cluster", pcluster  pcluster is either "local" or "redraider R2020b" or "threads" (experimental).
%                               default:  "RedRaider R2020b" on HPCC, "local" on everything else.
%                               NOTE:  threads won't work with functions that read or write to disk, 
%                                       so will not work with current ARRM_V2 code.  icsf 1/22.
%       "nnodes", nnodes        # of nodes to run on (default: 1).  Use > 1 only on systems that use Matlab Parallel Server, such as HPCC. 
%       "mb_per_core", mb_per_core  Amount of memory, in MBytes, to use per core.  range 2000/nworkers - 515625/nworkers.   Applies only for runs using MPS at HPCC.
%                                   (should be set in the SBATCH script for use on the local cluster, not here.
%                                           
%       "globalAtts",global_attributes_filename      
%                           name of file with any global attributes you might want in the netcdf metadata.
%                                   Each line is attribute_name value
%                                   This info will be written to the output netcdf's global attributes file, so this
%                                   can be used to add extra identifying information to the file
%                                   Use [] if no additional global attributes are desired.
%       "do_log",val         val is 0 to 3.     
%                                   default:  3
%                                   
%                                   0: write log info to console only.  
%                                   1: append log info to log file '...log'
%                                   2: append log info to log file, then copy to console at end
%                                   3: write  log info to log file (overwrite), then copy to console at end
%                                       generally, for production without parfor, use 3.
%                                       
%                                   4: write  log info to log file (overwrite), then copy to console at end.
%                                               generally, if using parfor, use 4.
%
%       "gridbox", [...]    provide the lat/lon range for 1 or more gridboxes, rather than have ARRM_V2_wrapper
%                               calculate them.  Useful for restarting a lat/lon run.
%       "do_continue", t/f  boolean.  If true, only runs gridboxes with no output file yet or which don't have completion_status set to 1 (completed) or 2 (completed, but no valid data).
%       "do_extended", t/f  boolean.  If true, Saves z-score values for each data point in a separate variable.
%
%       % this will need rethinking, Ian.
% %       "start_ix", start_ix    starting and ending indexes of stnIDs to process.
% %       "end_ix", end_ix
% %                               Use with caution with "do_continue" set to true.  If do_continue is true, code will
% %                               look for stnID of start_ix in output file, and start overwriting at that point.
% %                               If start_ix missing, or stnID of start_ix not found, will simply append to existing file.
% %                                   NOTE:  if "stations" supplied, then start_ix and end_ix are indexes within list of
% %                                          stations, not within entire siteTable.
% %       "start_stnID", stnID    station ID of where to begin or continue a run  (alternative to start_ix)
% %       "end_stnID, stnID       ending stnID (alternative to end_ix)
%
%
%       "llgrid_size"           grid size (in degrees) to cover for minifiles.  
%
%       "exit_on_error", false  boolean.  If true, run is aborted if any error is trapped
%                               If absent or false, error message is reported, but operation continues on to next loop
%                                   (station or gridcell)
%                               Generally, leave off (or set to false), except for testing/debugging.
%       "plotFlags",[...]       vector of flags to select plots to create for each site.  Do not use for gridded data
%                                   or for station runs with more than a few stations.
%                                   See function plot_ARRM_results(...) for info on what plots can be generated.
%       "figbase",figbase       starting figure number to use if plotting results.  Default:  1
%
%       Plus, any additional  optional varargins for ARRM_V2_run.  Those parameters are listed in ARRM_V2_run.m and
%       ARRM_V2_run_params_TTUCSC.m
%
%   Return values for retval:
%       0       success
%       1       thrown error caught, and exit_on_error was true.
%       2
%       3
%       4       
%       5       error in ARRM_V2_run_gridded(...) || ARRM_V2_run_stations(...)
%       6       thrown error caught, but exit_on_error was false.

    if (nargin < 7), help(mfilename('fullpath')); outnames=strings(0); retval=6; return; end
    
    % make sure the Matlab search path is set to find all the code we need.    
%    myhost = ARRM_V2_setpath(mfilename('fullpath'));   % don't need the hostname anymore in this program...
              ARRM_V2_setpath(mfilename('fullpath')); 

    [RP, DSP_base, Parms, logname, Unmatched] = initParams(model,ensemble,varname, scenario, base_yrs, hist_yrs, mdl_yrs, varargin{:});
    
    ARRM_V2_version_info(true);     % update the ARRM_V2 version info for the child processes.
    
    fprintf("%s\n", Parms.run_description);
    wrapper_log(logname, "%s\n", Parms.run_description);
    
    errors=false; 
         
    try
        reported = false;   % flag for error reporting...
        
        % generate list of lats, lons from range & stepsize
        
        outnames = cell(Parms.ngrids,1);
        DSP_run = cell(Parms.ngrids,1);         % this is only to suppress a Matlab warning related to temporary variables inside a parfor loop...
        
        if (Parms.do_parfor)      % run in parallel

            
            try
                    % get a jobname from the output filename, removing "downscaled", gridbox and extension ".nc"
                jobname=strsplit(string(DSP_base.outname),".");
                if (length(jobname)>10)
                    jobname=join(jobname(2:end-3),".");
                else
                    jobname=strings(0);
                end
%               ppool = ARRM_V2_start_workers2(par_info.nworkers, par_info.nnodes, par_info.cluster);
                ppool = ARRM_V2_start_workers2(Parms.nworkers, Parms.nnodes, Parms.cluster, jobname, Parms.mb_per_core);      % start parallel pool w/ max # of workers if not running yet.
            catch me
                report_me_error(me);
                error("error:  cannot start workers:  nworkers: %d   cluster:  %s   nnodes: %d mb_per_core: %d", Parms.nworkers, Parms.cluster, Parms.nnodes, Parms.mb_per_core);
                
            end
            
            fprintf("\n------------------------\n\n");
            fprintf(             "ARRM_V2_wrapper():  successfully started %d workers on %d nodes on %s\n", Parms.nworkers, Parms.nnodes, Parms.cluster);
            wrapper_log(logname, "ARRM_V2_wrapper():  successfully started %d workers on %d nodes on %s\n", Parms.nworkers, Parms.nnodes, Parms.cluster);
            fprintf("ppool is:\n");
            disp(ppool);
            fprintf("\n------------------------\n\n");
            
            lastgrid = false;
            start_tic = tic;
            retvals = zeros(Parms.ngrids,1);  % return values from each gridbox.
            parfor g=1:Parms.ngrids       % for each grid box.
                try
                    t = getCurrentTask(); 
                    workerID = t.ID;
                                    % tailor data params to specific run and gridbox
                    [DSP_run{g}, runID] = setup_run(Parms, DSP_base, g, lastgrid); %#ok<PFOUS>
                    outnames{g} = DSP_run{g}.outname;
                    llgrid_lbl = DSP_run{g}.llgrid_lbl;
                    if (islogical(DSP_base.debug_flag))
                        DSP_run{g}.debug_flag = get_debug_info(DSP_run{g}.gridbox);
                    end
                    wrapper_log(logname, "grid start worker %3d: %6d of %6d  %s: (%4d, %4d) to (%4d, %4d) %s %s\n", workerID, g, Parms.ngrids, datestr(now,'yyyy-mm-dd HH:MM:SS'), Parms.gridbox(g,:), llgrid_lbl, outnames{g});
%                   wrapper_log(logname, "grid start: %6d of %6d  %s\n", g, Parms.ngrids, datestr(now,'yyyy-mm-dd HH:MM:SS'));
                    if (~isempty(DSP_run{g}))
                        if (Parms.isStationRun)
                            
                            retvals(g) = ARRM_V2_run_stations('RP', RP, 'DSP', DSP_run{g}, "do_log", Parms.do_log, 'runID',   runID,                                                         "run_description", Parms.run_description, 'do_append',  false, Unmatched);     % PARFOR!  make changes to serial mode as well, Ian!
                        else
                            retvals(g) = ARRM_V2_run_gridded( 'RP', RP, 'DSP', DSP_run{g}, "do_log", Parms.do_log, 'runID',   runID, "latrange", Parms.latrange, "lonrange", Parms.lonrange, "run_description", Parms.run_description, Unmatched);
                        end
                        fids=fopen("all");
                        numfiles=numel(fids);
                        wrapper_log(logname, "grid   end worker %3d: %6d of %6d  %s %s %s numfiles open: %d ARRM_V2_run_gridded returned %d\n", workerID, g, Parms.ngrids, datestr(now,'yyyy-mm-dd HH:MM:SS'), llgrid_lbl, outnames{g}, numfiles, retvals(g));
                    else
                        outnames{g} = sprintf("gridbox %4d: %s:  no data to process %s\n", g, runID, llgrid_lbl);
                        wrapper_log(logname,"\tgridbox %4d: %s:  no data to process %s\n", g, runID, llgrid_lbl);
                    end
                catch me_par
                    fprintf('--------ARRM_V2_wrapper:  Error (parfor) processing %s:  %s\n', DSP_run{g}.runID, mfilename);
                    fprintf('ARRM_V2_wrapper:  error:  %s\t\t%s\n', me_par.identifier, me_par.message);
                    msgtext=getReport(me_par);
                    fprintf('%s\n', msgtext);
                    fprintf('------\n');
                    if (DSP_run{g}.exit_on_error)
 %                      reported = true;
                        rethrow(me_par); 
                    end
                    retvals(g) = 5;
                end
            end
            
            retval = max(retvals);
            
        else        % run in series
            fprintf('running in serial mode (not parfor)\n');
            first_grid = true;
            lastgrid = false;
            start_tic = tic;
            retvals = zeros(Parms.ngrids,1);
            for g=1:Parms.ngrids
                try
                    do_append = Parms.do_append & ~first_grid;
                    if (g==Parms.ngrids), lastgrid = true; end
                    [DSP_run{g}, runID] = setup_run(Parms, DSP_base, g, lastgrid); 
                    outnames{g} = DSP_run{g}.outname;
                    llgrid_lbl = DSP_run{g}.llgrid_lbl;
                    if (islogical(DSP_base.debug_flag))
                        DSP_run{g}.debug_flag = get_debug_info(DSP_run{g}.gridbox);
                    end
                    if (Parms.do_continue)
                        done = check_completion_status(outnames{g});
                        if (done), continue; end
                    end
                    wrapper_log(logname, "grid start: %6d of %6d  %s: (%4d, %4d) to (%4d, %4d), %s %s\n", g, Parms.ngrids, datestr(now,'yyyy-mm-dd HH:MM:SS'), Parms.gridbox(g,:), llgrid_lbl, outnames{g});
                    if (~isempty(DSP_run{g}))
                        if (Parms.isStationRun)
%                           ARRM_V2_run_stations('RP', RP, 'DSP', DSP_run{g}, "do_log", Parms.do_log, "runID",   runID, 'do_continue',Parms.do_continue);   % PARFOR!  make changes to serial mode as well, Ian!
                            ARRM_V2_run_stations(           'RP', RP, 'DSP', DSP_run{g}, "do_log", Parms.do_log, "runID",   runID,                                                                                      "run_description", Parms.run_description, "do_append",  do_append, Unmatched);     % PARFOR!  make changes to serial mode as well, Ian!
                        elseif (~istable(Parms.lonrange))
                            if (~isempty(Parms.testsites))
                                error("this needs updating still, Ian");
                                ARRM_V2_test_gridded_sites( 'RP', RP, 'DSP', DSP_run{g}, 'do_log', Parms.do_log, 'runID',   runID, "lats", Parms.lats,        "lons", Parms.lons,         "testsites", Parms.testsites, "run_description", Parms.run_description, Unmatched); %#ok<UNRCH>
                            else
                                ARRM_V2_run_gridded       ( 'RP', RP, 'DSP', DSP_run{g}, 'do_log', Parms.do_log, 'runID',   runID, "latrange", Parms.latrange, "lonrange", Parms.lonrange,                              "run_description", Parms.run_description,  Unmatched);
                            end
                        else   % hmmm.  not sure what I was planning here...  ability to downscale just a few specific lats & lons, I think.  ARRM_V2_test_gridded_sites(...) should handle that.
%                             stn_tbl = Parms.lonrange;
%                             latrange = stn_tbl.lats;
%                             lonrange = stn_tbl.lons;
%                             ARRM_V2_run_gridded( 'RP', RP, 'DSP', DSP_run{g}, 'do_log', Parms.do_log, 'runID',   runID, "latrange", Parms.latrange, "lonrange", Parms.lonrange, "stn_tbl", stn_tbl, Unmatched);
                            error("oops.  not implemented yet");
                        end
                        fprintf("gridbox netcdf file:  %s\n", outnames{g});
                        first_grid = false;
                        wrapper_log(logname, "grid   end: %6d of %6d  %s %s %s\n", g, Parms.ngrids, datestr(now,'yyyy-mm-dd HH:MM:SS'), llgrid_lbl, outnames{g});
                    else
                        outnames{g} = sprintf("gridbox %4d: %s:  no data to process %s\n", g, runID, llgrid_lbl);
                        wrapper_log(logname,"\tgridbox %4d: %s:  no data to process %s\n", g, runID, llgrid_lbl);
                    end
                 catch me
                    fprintf('ARRM_V2_wrapper:  error:  %s\t\t%s\n', me.identifier, me.message);
                    msgtext=getReport(me);
                    fprintf('error traceback:\n%s\n', msgtext);
                    fprintf('--------ARRM_V2_wrapper:  Error processing %s:  %s\n', DSP_run{g}.runID, mfilename);
                    fprintf('------\n');
                    if (DSP_run{g}.exit_on_error)
                        reported = true;
                        rethrow(me); 
                    end
                    retvals(g) = 4;
                end
            end
            retval = max(retvals);
        end
        
            % probably could add some code here to merge the output minifiles into 1 netcdf file
            % and delete the minifiles.  Later!
            
    catch me
        if (~reported || ~contains(me.identifier,'ARRM'))
            fprintf('ARRM_V2_wrapper initialization failure:--------%s: caught exception----------\n', mfilename);
            fprintf('ARRM_V2_wrapper:  error:  %s\t\t%s\n', me.identifier, me.message);
            fprintf('traceback:\n');
            msgtext=getReport(me);
            fprintf('%s\n', msgtext);
            fprintf('------\n');
            errors=true;
        end
        retval = 1;
    end
    if (exist("start_tic","var"))
        elapsed = toc(start_tic);
        elapsed_str = timestr(elapsed);
    else
%       elapsed = 0;
        elapsed_str = "(unknown)";
    end
    if (errors)
        fprintf('%s:  run completed, errors encountered.  total elapsed time:  %s for %4d output files\n', mfilename, elapsed_str, Parms.ngrids);
        wrapper_log(logname, '%s:  run completed, errors encountered.  total elapsed time:  %s for %4d output files\n', mfilename, elapsed_str, Parms.ngrids);
        retval = 6;
    else
        fprintf('%s:  total elapsed time:  %s for %4d output files\n', mfilename, elapsed_str, Parms.ngrids);
        wrapper_log(logname, '%s:  run completed.  total elapsed time:  %s for %4d output files\n', mfilename, elapsed_str, Parms.ngrids);
    end
    % close any open netcdf or log files
    fids=fopen('all');
    if (~isempty(fids))
        for j=1:length(fids)
            fname = fopen(fids(j));
            [~,~,ext] = fileparts(fname);
            if (strcmpi(ext,'.nc') || strcmpi(ext,'.log'))
                fprintf('-------------- ARRM_V2_wrapper:  closing unexpected open file %s\n', fname);
                fclose(fids(j));
            end
        end
    end
    
    if (exist("ppool","var") && ~isempty(ppool))
        ppool.delete();
    end
        % quit instead of returning if we're on hpcc system.
        % all their names end with ".local'.  however, so do my computers
        % when it is at home.  so make sure we're not on any of my computers first

%    if (~contains(myhost,'icsf') )
%        if (contains(myhost,'.local') || contains(myhost,'compute') )
%            quit;
%        end
%    end

end


function debug_info = get_debug_info(gridbox)

    if (all(gridbox(1:2)==[25,mod(-100,360)]))
        debug_info = [1,1; 1,5; 25, -100];
    elseif (all(gridbox(1:2)==[25,mod(-105,360)]))
        debug_info = [1,1; 12,16; 25, -105];
%     elseif (all(gridbox(1:2)==[25,-105]))
%         debug_info = [19,19; 12,16; 35, -115];
    else
        debug_info=[];
    end
end



function [stninfo, Parms, DSP] = get_station_info(Parms, DSP)
%  returns stninfo for all stations selected by input parameters.
%
%   Parms.siteTable is one of:
%       Parms.siteTable  (if it's a QC netcdf station table)
%       stninfo read from file named Parms.siteTable (if it's the name of a netcdf QC station file)
%       stninfo read from DSP's obsvname file, after parameter substitution.
%       name of a QC netcdf station table
%
%   Parms.stations is one of:
%       "all"
%       name of a csv file with column "stnID"
%       vector of strings with stnIDs.
%
%   If Parms.stations is not empty & != "all", then only stations in the list and in the stninfo are kept
%   if Parms' start_ix, end_ix, start_stnID or end_stnID are specified, list of stations is limited to given range.
%   If Parms.lats & Parms.lons are not empty, then only stations in specified range are kept.

    siteTable = Parms.siteTable;
    if (isempty_s(siteTable) && isstring(Parms.stations) && length(Parms.stations)==1 && isfile(Parms.stations))
        siteTable = Parms.stations;
        Parms.stations = strings(0,0);
    end
    
    if (isempty(siteTable))
        siteTable = fullfile(DSP.obsdir, DSP.obsnames);
    end
 
    if (ischar_s(siteTable))
        if (~isQCnetcdf(siteTable)), error('siteTable %s is not a QC netcdf station file', siteTable); end
    end
    
    if (isempty(Parms.latrange)), Parms.latrange = [min(Parms.lats), max(Parms.lats)+1e-6]; end
    if (isempty(Parms.lonrange)), Parms.lonrange = [min(Parms.lons), max(Parms.lons)+1e-6]; end

    siteTable = QC_get_site_table(siteTable, ...
                                  Parms.base_yrs(1), Parms.base_yrs(2), "fullYears", false,  "minYears",Parms.min_obs_yrs, ...
                                  "stnID", Parms.stations, "latrange", Parms.latrange, "lonrange", Parms.lonrange ...
                                 );
                             
     siteTable.lon = mod(siteTable.lon,360);    % make sure lon's are in range 0-360, not -180 to 180.
        
    if (isQCstntbl(siteTable))
        stninfo = siteTable; 
    else
        error("can't determine stninfo for run");
    end
    
        % What follows is to allow running a subset of stations, mainly for restarting a partially complete run.
        % User can specify start_stnID, etc..  Otherwise, we'll use all the stations specified by Parms.stations.
                            % if user specified start_ix, end_ix, start_stnID or end_stnID, start by limiting to that
                            % range of stations
    if (~isempty_s(Parms.start_stnID))
        ix = find(strcmp(Parms.stations, Parms.start_stnID),1);
        if (isempty(ix)), error("error:  starting station ID %s not found in station list", Parms.start_stnID); end        
        if (~isempty(Parms.start_ix))
            Parms.start_ix = max(Parms.start_ix, ix);
        else
            Parms.start_ix = ix;
        end
    end
    if (~isempty_s(Parms.end_stnID))
        ix = find(strcmp(Parms.stations, Parms.end_stnID),1);
        if (isempty(ix)), error("error:  ending station ID %s not found in station table", Parms.end_stnID); end        
        if (~isempty(Parms.end_ix))
            Parms.end_ix = min(Parms.end_ix, ix);
        else
            Parms.end_ix = ix;
        end
    end
            % keep only the keepers!
    if ((~isempty(Parms.start_ix) && Parms.start_ix > 0) || (~isempty(Parms.end_ix) && Parms.end_ix > 0))
        keepers = true(size(Parms.stations,1),1);
        if (Parms.start_ix > 0), keepers(1:(Parms.start_ix-1)) = false; end
        if (Parms.end_ix > 0),   keepers((Parms.end_ix+1):end) = false; end
        keepix = find(keepers);
        [~,~,keepix] = intersect(keepix, stninfo.index);
        if (isempty(keepix)), error("no stations in start/end range"); end
        stninfo = stninfo(keepix);
    end
            
    if (size(stninfo,1)==0), error("no stations in lat/lon range"); end
    
%                   keep only stations with enough data to start with.
    keepers = check_years_available(stninfo, Parms.base_yrs, DSP.monthly_valid_count);
    
    if (isempty(keepers))
        error("no stations with sufficient monthly data in input specs");
    else
        stninfo = stninfo(keepers,:);
    end
    
        % put the stninfo and stations list back into Parms.
    Parms.stninfo = stninfo;
    Parms.stations = stninfo.stnID;
    
    DSP = DSP.update('stninfo', stninfo, DSP.Unmatched); 
end

function [Parms, DSP] = get_gridded_info(Parms, DSP)

    % empty placeholder 
    % in case we need to read some stuff from the input files before we go too far...

end


function [RP, DSP_base, Parms, logname,  Unmatched] = initParams(model,ensemble,varname, scenario, base_yrs, hist_yrs, mdl_yrs, varargin)

    % returns an ARRM_V2_RunParams and an ARRM_V2_DataParams object with settings from input arguments, along with
    % parameters for use locally.
        
    model = string(model);
    ensemble = string(ensemble);
%   if (mod(length(varargin),2) ~= 0), error("%s:  error: missing argument for keyword/argument parameters", mfilename); end
        
                    % parse input for DA_title
    p = inputParser;
    p.StructExpand = true;
    p.KeepUnmatched = true;        
    
%     dv = datevec(now);
%     curyear = dv(1);
         
                    % these are the params we want to handle in ARRM_V2_wrapper
    addRequired(p,"model",      @(s)ischars(s));
    addRequired(p,"ensemble",   @(s)ischars(s));
    addRequired(p,"varname",    @(s)ischars(s) && length(string(varname)) == 1);
    addRequired(p,"scenario",   @(s)ischars(s));
    addRequired(p,"base_yrs",   @(s)isnumeric(s) && length(s)==2 && s(1) >= 1850 && s(2) <= 2038);      % 2038 so we can run GFDL perfect model data.
    addRequired(p,"hist_yrs",   @(s)isnumeric(s) && length(s)==2 && s(1) >= 1900 && s(2) <= 2038);
    addRequired(p,"mdl_yrs",    @(s)isnumeric(s) && length(s)==2 && s(1) >= 1900 && s(2) <= 2115);      % 2115 so we can run perfect model data

    addParameter(p,"latrange",[],     @(s) isempty(s) ||  isnumeric(s));       % numeric, for lat/lon limited run
    addParameter(p,"lonrange",[],     @(s) isempty(s) ||  isnumeric(s) || istable(lonrange));       % numeric, for lat/lon limited run
    addParameter(p,"lats",[],     @(s) isempty(s) ||  isnumeric(s));       % numeric, for gridded
    addParameter(p,"lons",[],     @(s) isempty(s) ||  isnumeric(s));       % numeric, for gridded
    addParameter(p,'llgrid_size',[], @(s) isnumeric(s) && (isempty(s) || (numel(s)==2 && ~any(s<1) && ~any(s>45))));
    addParameter(p,"obsvname",[]);                      % if not a standard obs. varname.
    addParameter(p, "testsites",strings(0));            % string array of test sites to test gridded on. (very non-standard run...results stored in mat files.)
    addParameter(p,"stations",[]);                      % "stations", for list of sites (for station run)
    addParameter(p,"siteTable",[]);                     % or site table or name of station netcdf file
    addParameter(p,"gridTable",[]);                     % or grid table or name of station netcdf file
    addParameter(p,"min_obs_yrs",30, @(s) isnumeric(s) && s >= 0 && mod(s,1)==0); % keep only stations with min # of years of obs data. 
    addParameter(p,"gridbox", [], @(s) isempty(s) || (isnumeric(s) && size(s,2)==4) || ischars(s));   %specify one or more gridboxes.  [lat1,lon1, lat2, lon2;...]
    
    addParameter(p,"do_parfor",false, @(s) islogical(s) || (isnumeric(s) && (s==0 || s==1)));
    addParameter(p,"do_log",3, @(s) isnumeric(s));
    addParameter(p,"nnodes",1, @(s) isnumeric(s) && s>0 && s <=4);
    addParameter(p,"nworkers",0, @(s) isnumeric(s));
    addParameter(p,"mb_per_core",0, @(s) isnumeric(s));
    addParameter(p,"cluster",strings(0), @(s) isempty(s) || strlength(s)==0 || any(strcmp(s, ["local","redraider R2020b","threads"])));
    addParameter(p,"exit_on_error", false, @(s) islogical(s) || (isnumeric(s) && any(s==[0,1])));
    addParameter(p,"do_append",false);                  % only valid for station runs.  Check that it still works, Ian!
    
    addParameter(p,"do_continue",false, @(s) islogical(s) || (isnumeric(s) && s >= 0 && s <= 1));   % set to true to continue a previous run.
                                                                                                % will check completion status of all gridboxes and remove gridboxes which are flagges as complete.
                                                                                                %  USE WITH CAUTION!
    addParameter(p,"do_extended",false, @(s) islogical(s) || (isnumeric(s) && s >= 0 && s <= 1));   % set to true to write z-score values to output netcdf file.
    addParameter(p, "run_description",strings(0));  % text to write at start of log files
    addParameter(p, "pdf_map_method",strings(0));   % will default to "clim" for precip and "linear" otherwise.

%     addParameter(p,"start_stnID","", @(s) ischar_s(s));
%     addParameter(p,"end_stnID","", @(s) ischar_s(s));
%     addParameter(p,"start_ix",0, @(s) isnumeric(s) && s > 0 && mod(s,1)==0);
%     addParameter(p,"end_ix",0, @(s) isnumeric(s) && s >= 0 && mod(s,1)==0);
                % all other params will be passed on in Unmatched struct.
                
%   IAN:  you need to add in some defaults for:
%       pdf_map_method  ("clim"), clim_nterms (5), clim_sig_terms(2), anom_nterms (5), anom_sig_terms (2) 
%       
        
    parse(p, model, ensemble, varname, scenario, base_yrs, hist_yrs, mdl_yrs, varargin{:});
    Parms = p.Results;
    Unmatched = p.Unmatched;        % save rest of input params for later
    
    Parms.model = string(Parms.model);
    Parms.ensemble = string(Parms.ensemble);
    Parms.varname = string(Parms.varname);
    Parms.scenario = string(Parms.scenario);
    Parms.obsvname = string(Parms.obsvname);
    Parms.do_append = ~Parms.do_parfor && Parms.do_append; % isempty(Parms.gridbox);
    
        % hard-coded for now.
    Parms.start_stnID = "";
    Parms.end_stnID = "";
    Parms.start_ix = 0;
    Parms.end_ix = 0;
    Parms.do_append = false;
    
    if (isempty(Parms.cluster) || strlength(Parms.cluster) == 0)
        [~,on] = get_hostname(true);
        if (on.hpcc_system)
            Parms.cluster = "redraider R2020b";
        else
            Parms.cluster = "local";
        end
        if (isempty(Parms.nworkers) || Parms.nworkers == 0)
            if (on.hpcc_system)
                Parms.nworkers = 128*Parms.nnodes;
            elseif (on.neys)
                Parms.nworkers = 18;
            elseif (on.dev_system)
                Parms.nworkers = 8;
            end
        end
    end
    
        % for now, error if hist_yrs ~= base_yrs
        
    if (Parms.hist_yrs ~= Parms.base_yrs)
        error("error: for now, hist_yrs ([%d,%d]) must match base_yrs ([$d,%d])", Parms.hist_hrs, Parms.base_yrs);
    end
    
%     if (~isempty(Parms.lons)) 
%         Parms.lons = mod((Parms.lons),360);      %make sure lons are in range [0,360]);
%         if (Parms.lons(end)==0), Parms.lons(end)=360; end
%     end
%     

    Parms.runType = "Downscaling";
    if (~isempty(Parms.stations) || ~isempty(Parms.siteTable))
        Parms.isStationRun = true;
        Parms.runType(2)="Station";
        if (isempty(Parms.llgrid_size)), Parms.llgrid_size = [5,5]; end
%         if (~isempty(Parms.stations))
%             Parms.siteTable = QC_get_site_table(Parms.siteTable,"stnID", Parms.stations, "loadTblVars",false);
%         end
        if (isempty(Parms.obsvname)), Parms.obsvname = QC_station_varname(Parms.varname); end
        DSP_base = ARRM_V2_DownscalingParams("runType",Parms.runType, "exit_on_error", Parms.exit_on_error,...
                                            "do_parfor",Parms.do_parfor, "do_append", Parms.do_append, ...
                                            "model",Parms.model,"ensemble",Parms.ensemble,"scenario",Parms.scenario, ...
                                            "varname", Parms.varname, 'obsvname',Parms.obsvname, ...
                                            "base_yrs", Parms.base_yrs, "hist_yrs",Parms.hist_yrs,"rolling_yrs",Parms.mdl_yrs,...
                                            "dspType","global", "llgrid_size", Parms.llgrid_size, ...
                                            "extended", Parms.do_extended, ...
                                            Unmatched); 
                                
                    % set up basic data params.
%       Parms.obsvname = QC_station_varname(Parms.varname);        % might not need this, Ian...
        [stninfo, Parms, DSP_base] = get_station_info(Parms, DSP_base);
        
                    % figure out what stations to keep.
                    
                        % if input specified latrange and lonrange, find stations within that range.
        if (~isempty(Parms.latrange) || ~isempty(Parms.lonrange))
            if (isempty(Parms.latrange)), error("error:  lonrange specified, but latrange not given"); end
            if (isempty(Parms.lonrange)), error("error:  latrange specified, but lonrange not given"); end
            
            Parms.lonrange = longitude_360(Parms.lonrange);
            
            if (Parms.latrange(1) <= Parms.latrange(end))
                latkeepers = stninfo.lat >= Parms.latrange(1) & stninfo.lat <= Parms.latrange(end);
            else
                latkeepers = stninfo.lat <= Parms.latrange(1) | stninfo.lat >= Parms.latrange(end);
            end
            if (Parms.lonrange(1) <= Parms.lonrange(end))
                lonkeepers = stninfo.lon >= Parms.lonrange(1) & stninfo.lon <= Parms.lonrange(end);
            else
                lonkeepers = stninfo.lon <= Parms.lonrange(1) | stninfo.lon >= Parms.lonrange(end);
            end
            keepers = latkeepers & lonkeepers;
            Parms.lats = stninfo.lat(keepers);
            Parms.lons = stninfo.lon(keepers);
        else    
                    % input specified lats & lons.  Set latrange and lonrange from range from last & lons.
                    
            if (isempty(Parms.lats))
                if (~isempty(Parms.lons))
                    error("lats & lons must be same length (or both empty)");
                end
                Parms.lats = stninfo.lat; 
                Parms.lons = stninfo.lon;
            elseif (length(Parms.lats) ~= length(Parms.lons))
                error("lats & lons must be same length (or both empty)");
            end
            Parms.latrange = [min(Parms.lats), max(Parms.lats)];
            Parms.lonrange = [min(Parms.lats), min(Parms.lons)];
        end
        
        
    else        % for gridded run
        
        Parms.runType(2)="Gridded";
        Parms.isStationRun = false;
        Parms.do_append = false;
        if (isempty(Parms.llgrid_size)), Parms.llgrid_size = [5,5]; end
        DSP_base = ARRM_V2_DownscalingParams("runType",Parms.runType, "exit_on_error", Parms.exit_on_error, ...
                                            "do_parfor",Parms.do_parfor, ...
                                            "model",Parms.model,"ensemble",Parms.ensemble,"scenario",Parms.scenario, ...
                                            "varname", Parms.varname, 'obsvname',Parms.obsvname, ...
                                            "base_yrs", Parms.base_yrs, "hist_yrs",Parms.hist_yrs,"rolling_yrs",Parms.mdl_yrs,...
                                            "dspType","global", "llgrid_size", Parms.llgrid_size, ...
                                            "extended", Parms.do_extended, ...
                                            Unmatched); 
                                        
                    % set up basic data params.
        [Parms, DSP_base] = get_gridded_info(Parms, DSP_base);
        
                    % get list of lats/lons to use.  If user specified lats & lons, use them.
                    % Otherwise, scan input obs files for full list of lats & lons
                    % Include only lats & lons within range if latrange or lonrange is specified.
                    % we need the obs lats & lons to set up the gridboxes.
                    
                    % get the observation & model lats & lons from the first obs & model file
        [obslats, obslons] = ncdf_get_latlons(fullfile(DSP_base.obsdir, DSP_base.obsnames{1}));
        obslons = mod(obslons, 360);
        
        if (isempty(Parms.lats) && isempty(Parms.latrange))
            Parms.latrange = [-90,90];
            [~,Parms.lats] = lat_region(Parms.latrange, obslats,[-1,-1], true);     % get all the obslats within the range.  Allow wrapping (if latrange(1) > latrange(2), then gets the poles and excludes the middle latitudes)
        elseif (isempty(Parms.latrange))
            Parms.latrange = [Parms.lats(1), Parms.lats(end)];
        else
            [~,Parms.lats] = lat_region(Parms.latrange, obslats,[-1,-1], true);     % get all the obslats within the range.  Allow wrapping (if latrange(1) > latrange(2), then gets the poles and excludes the middle latitudes)
        end
        
        if (isempty(Parms.lons) && isempty(Parms.lonrange))
            Parms.lonrange = [0,360];
           [~, Parms.lons] = lon_region(Parms.lonrange, obslons, [-1,-1], true);  % get all the obslons within the range, allow wrapping (across 360->0 meridian)
        elseif (isempty(Parms.lonrange))
            Parms.lonrange = [Parms.lons(1), Parms.lons(end)];
        else
           [~, Parms.lons] = lon_region(Parms.lonrange, obslons, [-1,-1], true);  % get all the obslons within the range, allow wrapping (across 360->0 meridian)
        end
        
    end
    
        % make sure lats & lons are double
    Parms.lats = double(Parms.lats);
    Parms.lons = double(Parms.lons);
    Parms.latrange = double(Parms.latrange);
    Parms.lonrange = double(Parms.lonrange);
    
    Unmatched = DSP_base.Unmatched;
    
    RP = ARRM_V2_RunParams(Parms.varname, "runType",Parms.runType,'do_problines',false, 'do_daily_trends',false, Unmatched);   % will need to turn off do_pdfs for gridcells...
    Unmatched = RP.Unmatched;
%     if (~isempty(fieldnames(Unmatched)))
%         fields=fieldnames(Unmatched);
%         fprintf("error:  unknown input parameters:\n");
%         for i=1:length(fields)
%             fprintf("\t%s\n", fields{i});
%         end
%         error("ARRM_V2_wrapper aborting");
%     end

%     if (DSP_base.isPrecipRun)
%         probs = [0.5000 0.8413 0.9000 0.9500 0.9772 0.9900 0.9938 0.9990 0.9999 1.0000 1];
%         RP = RP.update("trend_yr_flags", [], "trend_thresh", 0, "trend_order",[], "probs",probs);
%         if (isempty(RP.pdf_map_method)), RP.pdf_map_method = "clim"; end
%     else
%         if (isempty(RP.pdf_map_method)), RP.pdf_map_method = "linear"; end
%     end
    
    if (isempty(Parms.testsites))    
        if (isempty(Parms.gridbox))
    %       Parms.gridbox = make_gridboxes(Parms.latrange, Parms.lonrange, Parms.lats, Parms.lons, DSP_base.llgrid_size, Parms.isStationRun, Parms.do_parfor);
            Parms.gridbox = make_gridboxes(Parms.latrange, Parms.lonrange, Parms.lats, Parms.lons, DSP_base.llgrid_size, Parms.isStationRun);
        elseif (ischars(Parms.gridbox))
            grstring = string(Parms.gridbox);
            Parms.gridbox = zeros(length(grstring),4);
            for i=1:length(grstring)
                Parms.gridbox(i,:) = parse_gridbox_string(grstring(i), DSP_base.llgrid_size);
            end
        else
            Parms.gridbox = Parms.gridbox;
        end
    else
            % testing gridded run with a list of stations, so we'll treat it as one large gridbox.
        Parms.do_continue = false;
        Parms.do_parfor = false;
        Parms.gridbox = [floor(min(Parms.latrange)), floor(min(Parms.lonrange)), ceil(max(Parms.latrange)), ceil(max(Parms.lonrange))];
    end
    
    outdir=DSP_base.replace_keywords([], DSP_base.outdir);
%   logname=fullfile(outdir, sprintf("ARRM_V2_wrapper.%s.log",datestr(now,"yyyy-mm-dd.HH.MM.SS")));
    logname=fullfile(outdir, "ARRM_V2_wrapper.log");
    
    Parms.ngrids = size(Parms.gridbox,1);
    if (Parms.do_continue)
        Parms.gridbox = check_gridboxes_for_completion(Parms, DSP_base);
        Parms.ngrids = size(Parms.gridbox,1);
    end
    
    wrapper_log(logname, "ARRM_V2_wrapper run, started %s.  %d gridboxes\n\n", datestr(now,"yyyy-mm-dd HH:MM"), Parms.ngrids);
    
            % if continue flag set, keep only the gridboxes that aren't incomplete...
end

function gbox = parse_gridbox_string(gridbox,  llgrid_size)
    % converts strings like N05W110 to a gridbox.

        if (isempty_s(gridbox) || strncmpi(gridbox,"all",3))
            gbox=[-90,-180,90,180];
        else
            if (strncmpi(gridbox,"grid_",5))
                gridbox = extractAfter(gridbox, 5);
            end
            lat1=str2double(extractBetween(gridbox,2,3));
            c=char(extractBefore(gridbox,2));
            if (lower(c)=='s'), lat1 = -lat1; elseif (lower(c) ~= 'n'), error("error:  bad grid specification: %s", gridbox); end
            lon1=str2double(extractAfter(gridbox,4));
            c=char(extractBetween(gridbox,4,4));
            if (lower(c)=='w'), lon1 = -lon1; elseif (lower(c) ~= 'e'), error("error:  bad grid specification: %s", gridbox); end
            gbox = [lat1, lon1, lat1+llgrid_size(1), lon1+llgrid_size(2)];
        end
end


%        [DSP_run, runID] = setup_run(Parms, DSP, i);
function [DSP,     runID] = setup_run(Parms, DSP, g, lastgrid)

    if (Parms.isStationRun)
        [DSP, runID] = setup_station_run(Parms, DSP,  g, lastgrid);
    else
        [DSP, runID] = setup_gridded_run(Parms, DSP, g);
    end

end

function [DSP, runID] = setup_station_run(Parms, DSP_base, gridix, lastgrid)
% finds all stations within gridbox(gridix,:);  returns DSP with stninfo for only those inside the gridbox.
%  also sets the grid label for output filename, llgrid_lbl , and an updated runID for the log file.

        % already added so far:  
%   DSP_base = ARRM_V2_DownscalingParams("runType",runType, "exit_on_error", Parms.exit_on_error,"do_parfor",Parms.do_parfor,...
%                                   "base_yrs", Parms.base_yrs, "hist_yrs",Parms.hist_yrs,"rolling_yrs",Parms.mdl_yrs,...
%                                   "dspType","global", ...
%                                   Unmatched); 

    runID = sprintf("%s_%s_%s_%s_%s_%d_%d",Parms.model,Parms.ensemble, Parms.varname, Parms.scenario, DSP_base.region, gridix, Parms.ngrids);
    stninfo = Parms.stninfo;
    
    
            % if only 1 gridbox, it should encompass the entire dataset, so only set limit if more than 1 gridbox.
            % exception:  gridbox was set to "all"
        gridbox = Parms.gridbox(gridix,:);
    if (size(gridbox,1)>1 || ~isequal(gridbox,[-90,-180,90,180]))
        latmin = gridbox(1);
        lonmin = gridbox(2);
        latmax = gridbox(3);
        lonmax = gridbox(4);
        keepers = stninfo.lat >= latmin & stninfo.lat < latmax & stninfo.lon >= lonmin & stninfo.lon < lonmax;
%       disp([(stninfo.lat >= latmin & stninfo.lat < latmax)';(stninfo.lon >= lonmin & stninfo.lon < lonmax)';keepers']);
        stninfo = stninfo(keepers,:);
        if (Parms.do_append)
            llgrid_lbl = "combined";
        else
            llgrid_lbl = sprintf("grid_%s", llgrid_string(latmin, lonmin));
        end
%       fprintf('nkeepers:  %d\n', sum(keepers));
   else
        llgrid_lbl = "combined";
    end
    
    if (size(stninfo,1)==0)         % if no stations in gridbox, set DSP to empty and return.
        DSP=[];                     % (shouldn't happen...we weeded out the empty ones while making the gridboxes.)  
        return;
    end
    
    DSP = DSP_base.update("stninfo",stninfo, ... 
                          "llgrid_lbl", llgrid_lbl, ...
                           "lats",stninfo.lat,"lons",stninfo.lon, ...
                          "runID",runID, 'lastgrid', lastgrid);


%           Shouldn't need this, Ian.  S/B done by finalize.
%     if (isempty(DSP.mdlnames))
%         DSP = DSP.make_fnames("mdlnames");
%     end
%     if (isempty(DSP.histnames))
%         DSP = DSP.make_fnames("histnames");
%     end
% 
            % set the trend years if they weren't specified.  
    if (~DSP.isPrecipRun && isempty(DSP.trend_yrs))
        trend_yrs = max_yr_range(DSP.base_yrs, DSP.hist_yrs, DSP.rolling_yrs);
        DSP = DSP.update('trend_yrs',trend_yrs);
    end

    DSP = DSP.finalize("runID", runID);       % do parameter substitution, etc.
    
end

function [DSP, runID] = setup_gridded_run(Parms, DSP_base, gridix)

    % this still needs work, Ian.
        
    runID = sprintf("%s_%s_%s_%s_%d_%d",Parms.model,Parms.ensemble, Parms.varname, Parms.scenario, gridix, Parms.ngrids);

    if (size(Parms.gridbox,1) > 1)
        gridbox = Parms.gridbox(gridix,:);
        latmin = gridbox(1);
        lonmin = gridbox(2);
        latmax = gridbox(3);
        lonmax = gridbox(4);
        dlat   = latmax - latmin;
        dlon   = mod(lonmax - lonmin, 360);
        
        latkeepers = Parms.lats >= latmin & Parms.lats < latmax;
        lonkeepers = Parms.lons >= lonmin & Parms.lons < lonmax;
        lats = Parms.lats(latkeepers);
        lons = Parms.lons(lonkeepers);
    
        gbox = gridbox;
        if (Parms.do_append)
            llgrid_lbl = "";
        else
            llgrid_lbl = sprintf("grid_%s", llgrid_string(latmin, lonmin, dlat, dlon, true));
        end
    else
        lats = Parms.lats;
        lons = Parms.lons;
        gbox = [lats(1), lons(1), lats(end), lons(end)];
        dlat = lats(end) - lats(1);
        dlon = mod(lons(end) - lons(1), 360);
        llgrid_lbl = sprintf("all_%s", llgrid_string(lats(1), lons(1), dlat, dlon, true));
        
    end
    
        % already added so far:  
%   DSP_base = ARRM_V2_DownscalingParams("runType",runType, "exit_on_error", Parms.exit_on_error,"do_parfor",Parms.do_parfor,...
%                                   "base_yrs", Parms.base_yrs, "hist_yrs",Parms.hist_yrs,"rolling_yrs",Parms.mdl_yrs,...
%                                   "dspType","global", ...
%                                   Unmatched); 
    
    DSP = DSP_base.update("lat",lats,"lon",lons, ... 
                          "llgrid_lbl", llgrid_lbl, ...
                          "gridbox", gbox, ...
                          "runID",runID, ...
                          "model",Parms.model,"ensemble",Parms.ensemble,"scenario",Parms.scenario); 
                             
    DSP = DSP.finalize("runID", runID);       % do parameter substitution, etc.
    
end

% function gridboxes = make_gridboxes(latrange, lonrange, lats, lons, llgrid_size, isStationRun, do_randomize)
function   gridboxes = make_gridboxes(latrange, lonrange, lats, lons, llgrid_size, isStationRun)
% calculates the gridboxes, to break the processing down into minifiles
%   if llgrid_size is empty, or is greater than the range of lats & lons, then creates a single gridbox
%   Otherwise, creates a list of gridboxes (lower left lat/lon, upper right lat/lon), 
%   Grid will be located on even multiples of llgrid_size.

    if (isempty(latrange)), latrange = [lats(1), lats(end)]; end
    if (isempty(lonrange)), lonrange = [lons(1), lons(end)]; end
        % if llgrid_size is empty, then just make 1 gridbox
    if (isempty(llgrid_size))
%         gridbox=[minlat, minlon, maxlat, maxlon]; 
        gridboxes = [latrange(1),lonrange(1),latrange(end),lonrange(end)];      % THIS NEEDS TESTING STILL< IAN
    else

            % redo bounds to even multiples of the grid steps.
        latgridstep = llgrid_size(1);
        longridstep = llgrid_size(2);

            % find bounded area to process.


        if (latrange(1) > latrange(end))
            dlat = 180 - (latrange(1) - latrange(end));
        else
            dlat = latrange(end) - latrange(1);
        end
        if(lonrange(1) > lonrange(end))
            dlon = 360 - (lonrange(1) - lonrange(end));   % works for both 0-360 and -180 - 180.
        else
            dlon = lonrange(end) - lonrange(1); 
        end


            % if llgrid_size is greater than the actual area, just make 1 gridbox.
        if (latgridstep >= dlat && longridstep >= dlon)
    %         gridbox=[minlat, minlon, maxlat, maxlon]; 
            gridboxes = [latrange(1),lonrange(1),latrange(end),lonrange(end)];      % THIS NEEDS TESTING STILL< IAN
        else

            nlatgrid = ceil(180/latgridstep);
            nlongrid = ceil(360/longridstep);

            ngrids = nlatgrid * nlongrid;

                % put together the grid boxes.  [lower_left_lat, lower_left_lon, upper_right_lat, upper_right_lon] for each gridbox.
            gridboxes = zeros(ngrids,4);
            ng=0;

                % make sure lats & lons are in range [-90,90) and 90,360).
            lats = mod((lats+90),180)-90;
            lons = mod(lons,360); 

            for i=1:nlatgrid
                lat1 = -90 + (i-1)*latgridstep;
                lat2 = -90 +    i *latgridstep;
                lat2 = min(90, lat2);
                for j=1:nlongrid
                    lon1 = (j-1)*longridstep;
                    lon2 =    j *longridstep;
                    lon2 = min(360, lon2);
                    if (isStationRun)
                        keepers = lats >= lat1 & lats < lat2 & lons >= lon1 & lons < lon2;
                        nsites = sum(keepers);
                        if (nsites > 0)     % keep the gridbox if there are stations to process in it.
                            ng=ng+1;
                            gridboxes(ng,:) = [lat1,lon1,lat2,lon2];
                        end
                    else
                        latkeepers = lats >= lat1 & lats < lat2;
                        lonkeepers = lons >= lon1 & lons < lon2;
                        has_gridpts = any(latkeepers) && any(lonkeepers);
                        if (has_gridpts)     % keep the gridbox if there are stations to process in it.
                            ng=ng+1;
                            gridboxes(ng,:) = [lat1,lon1,lat2,lon2];
                        end
                    end
                end
            end

            gridboxes(ng+1:end,:) = [];     % remove all unused gridboxes.

%             if (do_randomize)           % randomize order of gridbox.  May fix a problem with parfor loop ordering
%                                         % parfor seems to queue up many loops to the processors at once;  if the same processor
%                                         % gets all the near-empty gridboxes, parfor doesn't seem to send it any more jobs
%                                         % when it finishes.  This  randomizes the order of the gridboxes, so that is less likely
%                                         % to happen.
%                 ng = size(gridboxes,1);                        
%                 orig_rng = rng(1,'multFibonacci');  % switch to fibonacci rng and seed it, so we always get the same ordering.
%                 gix=randperm(ng);
%                 gridboxes = gridboxes(gix,:);
%                 rng(orig_rng);
%             end
        end
    end

end

function keepers = check_years_available(stninfo, base_yrs, monthly_valid_count)
% removes any station which can't meet the minimum monthly valid count. 
% Default is 20 years' data for every 1-month period.

    start_dnum = datenum(base_yrs(1),1,1);
    end_dnum   = datenum(base_yrs(2),12,31);
    ndays = min(end_dnum, stninfo.endDate) - max(start_dnum, stninfo.startDate);    % # of years available for each station
    stn_valid_monthly_count = ndays .* stninfo.pctValid / 100.0 /12 ;                         % max possible value for station' smallest monthly valid count
    keepers = stn_valid_monthly_count >= monthly_valid_count;                       % keep only those which might have enough data. 
    
end

function wrapper_log(logname, fmt, varargin)
        % if not running parallel, does fprintf of varargins to console and to file
        % else just does fprintf to file.
    fid = fopen(logname,'a');
    fprintf(fid, fmt, varargin{:});
    fclose(fid);
    fprintf(fmt, varargin{:});
end

function gridboxes = check_gridboxes_for_completion(Parms, DSP_base)
    ng = size(Parms.gridbox,1);
    keepers = true(ng,1);
    ngood=0;
    nerrs=0;
    nmiss=0;
    nincp=0;
    fprintf("checking for completion status on %d files\n", ng);
    for i=1:ng
        DSP = setup_gridded_run(Parms, DSP_base, i);
        fn = DSP.outname;
        [~,bn,bext] = fileparts(fn);
        bn = strcat(bn,bext);
        [completed, completion_status] = check_completion_status(fn);
        if (isempty(completion_status))
            nmiss = nmiss+1;
            if (nmiss < 10)
                fprintf("\nmissing file      %4d:     %s\n", nmiss, fn);
            end
        elseif (completion_status==-1)
            nerrs = nerrs+1;
            fprintf("\nerrors in file %4d %4d:   %s\n", i, nerrs, bn);
        elseif (completion_status== 0)
            nincp = nincp+1;
            fprintf("\nincomplete file %4d %4d:  %s\n", i, nincp, bn);
        elseif (completion_status== 1)
            ngood = ngood+1;
            if (ngood < 10)
                fprintf("\ncompleted file:        %4d %s\n",ngood,  fn);
            end
        end
        if (completed)
            keepers(i) = false;
%           fprintf("completed: %d %d %s\n", i, ng, fn);
        end
        show_progress(i, ng);
    end
    fprintf("%4d files already completed\n", ngood);
    fprintf("%4d files incomplete\n", nincp);
    fprintf("%4d files not present\n", nmiss);
    fprintf("%4d files with error flag\n", nerrs);
    gridboxes = Parms.gridbox(keepers,:);   % keep only the gridboxes whose output files are either not present or not complete.
end

function [completed, completion_status] = check_completion_status(fname)
%   returns true if fname's completion_status attribute is set to 1.
%       also returns the value of completion_status, if present.
%       completion_status
%           1       run completed
%           2       run completed, but no data written.  (e.g., using livneh gridded and location is not over land.)
%           0       run not completed
%           -1      error encountered
%           []      (empty) file not found.
%       
    completion_status = [];
    completed = false;
    try
%         [status, retval] = system(sprintf("ncdump -h %s | grep status | awk '{print $3}'", fname)); % faster than opening the netcdf and reading the attribute...
%         if (status ~= 0 || isempty(retval))
%             completion_status = nan;
%         else
%             completion_status=str2double(strtrim(retval));
%         end
%         if (completion_status > 0)
%             completed=true;
%         end
        nc = ncdf(fname);
        completion_status = nc.getattvalue("completion_status");
        if (completion_status > 0)
            completed = true;
        end
    catch
    end
end
        
