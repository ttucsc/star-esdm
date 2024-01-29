% function [outname, DA_downs_out, DA_obss_out, DA_mdls_out, DA_hists_out, DA_resultss_disag, DA_obss_fut, nprocessed, nfigs] = ARRM_V2_run_gridded(varargin)
function error_flag = ARRM_V2_run_gridded(varargin)
%   COMMENTS NEED UPDATING, IAN!!!

%   Program to run all douwscaled locations given a gridded set of lats/lons

%   Data & Downscaling Parameters (filenames, etc.) are in DSP.
%   Run Parameters are in RP
%   Can be run standalone, or started via ARRM_V2_wrapper(...)
% 
%   Inputs:     
%       'RP',RP             RunParameters object. (optional)       Can be missing, if all DSP params are in varargin.
%       'DSP',DSP           DowscalingParameter object. (optional) Can be missing, if all DSP params are in varargin.
%       'do_log', 1-3       flag for whether to write to log file or not.
%       varargin            any other needed arguments.  See comments in ARRM_V2_wrapper for list of possible arguments
%       
%
%   Gets run params and input file info from control file, and for each lat/lon, reads data from
%   input netcdf files and interpolates to lat/lon from nearest four gridcells.
%

%   2019-10-27 icsf
%   2019-11-27 icsf:  needs to be modified to specify the lonrange properly to cross the meridian if needed. 
%   2022-09-19 icsf:  moved zvals to a separate output file.

%-----------------------------------------------------------------------------------

       % make sure the path is set to find all the code we need.      
    [myhost, myprog] = ARRM_V2_setpath(mfilename('fullpath')); 
    start_tic = tic;   
    error_flag = 1;
    
    if (nargin<2)
        help(mfilename('fullpath'));
        return;
    end
    
%       1.  Parse the data run parameters

	[RP_base, DSP_base, DSP_obss, DSP_hists, DSP_mdls, all_obsdata, all_histdata, all_mdldata, Parms, calendar] = initParams(varargin{:});
    outname  = DSP_base.outname;
%   Parms.nargout = nargout;
    
    if (DSP_base.do_parfor)
        t = getCurrentTask(); 
        workerID = t.ID;
    else
        workerID = 0;
    end
%
%       Downscales gridded model data using gridded obs data

    log_header_info(DSP_base, myhost, myprog, workerID, Parms.run_description);
    
%       For model grid rectangle
    mdllats = DSP_mdls.file_lats;
    mdllons = DSP_mdls.file_lons;
    nmdllats = length(mdllats);
    nmdllons = length(mdllons);
    nwritten      = 0;
    ninsufficient = 0;
    nproblems     = 0;
    nout          = 0;
    nprocessed    = 0;
    
    jgrid = 0;
    ngrids = (nmdllats-1)*(nmdllons-1);
    jloc = 0;
    nlocs = length(DSP_obss.lats)*length(DSP_obss.lons);
    nfigs = 0;
    is_appended = false;
    
            % create the output minifile file
    if (Parms.save_output)
        nc_out   = ncdf_init_grid(RP_base, DSP_base, myprog, calendar, false);
        if (DSP_base.extended)
            zvals_name = sprintf("%s_zvals", DSP_base.varname);
            [fdir,fname,fext] = fileparts(DSP_base.outname);
            fname = strrep(fname, DSP_base.varname, zvals_name);
            zvals_outname = fullfile(fdir, strcat(fname, fext));
            nc_zvals = ncdf_init_grid(RP_base, DSP_base.update("outname", zvals_outname), myprog, calendar, true);
        else
            nc_zvals = [];
        end
    else
        nc_out   = [];        
        nc_zvals = [];
    end

            % some arrays that were useful for testing/debugging.  Stuff is
            % commented out now.  icsf 6/21
%     DA_obss_out  = cell(nlocs,1);
%     DA_hists_out = cell(nlocs,1);
%     DA_mdls_out  = cell(nlocs,1);
%     DA_downs_out = cell(nlocs,1);
%     DA_obss_fut  = cell(nlocs,1);
%     DA_resultss_disag = cell(nlocs,1);
                               
    runLbl  = DSP_base.runLbl; % sprintf('%s %s %s %s', DSP.model, DSP.ensemble, DSP.varname, DSP.scenario);

    site_error_flag = false;
    grid_error_flag = false;
    run_gridded_error_flag = false;
    
    myfig = 1; %#ok<NASGU>
    
    try
                % work through each quadrant of the model gridpoints.
                % We'll do some of the disaggregation of the model data for each corner once (up to the histogramming),
                % then find all the obs locations inside the rectangle and downscale them
                % That saves the comput time of disaggregating the model gridpoints each time.
        for mdllat_ix = 1:(nmdllats-1)      % for each model lat point
    %       mlats = DSP_mdls.file_lats(latix:latix+1);       %lats for the rectangle we're looking at

                % get indexes of obs lats inside current gridbox
            looplats_ix = find(DSP_obss.lats >= mdllats(mdllat_ix) & DSP_obss.lats <  mdllats(mdllat_ix+1));
            nllats = length(looplats_ix);
            

            for mdllon_ix=1:(nmdllons - 1)      % for each model lon point
                 jgrid = jgrid+1;
    %           mlons = DSP_mdls.file_lons(lonix:lonix+1);    % lons for the rectangle we're looking at.

                    % get indexes of lons inside current gridcell quadrant.
                    
                if (mdllons(mdllon_ix+1) > mdllons(mdllon_ix))  % do we cross 360->0 or 180->-180?  
                    looplons_ix =  find(DSP_obss.lons >= mdllons(mdllon_ix) & DSP_obss.lons <  mdllons(mdllon_ix+1));
                else
                    looplons_ix =  find(DSP_obss.lons >= mdllons(mdllon_ix) | DSP_obss.lons <  mdllons(mdllon_ix+1));
                end                    
                nllons = length(looplons_ix);

                         
                nsites = nllats*nllons;
                if (nsites==0) 
                    if (DSP_base.do_parfor) 
                        fprintf(    "-------------------ARRM_V2_run_gridded %s grid %5d of %5d:   no obs sites in gridcell\n",   DSP_base.llgrid_lbl, jgrid, ngrids); 
                    end
                    DSP_base.print_log("-------------------ARRM_V2_run_gridded %s grid %5d of %5d: no obs sites in gridcell\n", DSP_base.llgrid_lbl, jgrid, ngrids);
                    continue; 
                end    
                

                try
                            % extract the four corners of the model gridcell quadrant
                    grid_error_flag = false;
                    [mdl_DAs_quadrant, hist_DAs_quadrant, DSP_mdl, insufficient_data] = get_quadrant_DAs(mdllat_ix, mdllon_ix, mdllats, mdllons, all_mdldata, all_histdata, RP_base, DSP_base, DSP_mdls, DSP_hists); 
                    
                    if (insufficient_data), continue; end
                    
                catch me

                    grid_error_flag = true;
                    nproblems = nproblems + 1;
                    DSP_base.warn_log('error: ARRM_V2_run_gridded:  on latix %d, lonix %d, run id %s\n', mdllat_ix, mdllon_ix, runLbl);
                    DSP_base.warn_log('\t\t******** error:  %s\t\t%s\n', me.identifier, me.message);
                    
                    if (DSP_base.exit_on_error), rethrow(me); end % ,------------bail out if we're supposed to exit on error.
                    
                    msgText = getReport(me);
                    DSP_base.warn_log('\nerror caught:\n%s\n',msgText);
                    continue 
                end
                
                            % now process each obs site within the model gridcell quadrant.
                for illat = 1:nllats
                    obslatix = looplats_ix(illat);
                    lat_pt = DSP_obss.lats(obslatix);
                    for illon = 1:nllons
                        jloc = jloc+1;
                        loop_tic = tic();
                        obslonix = looplons_ix(illon);
                        lon_pt = DSP_obss.lons(obslonix);
%                       merge_seed = latlon_seed(lat_pt, lon_pt);       % create a unique seed based on the lat & lon.
                        gridptID = sprintf("[%4d,%4d]", obslatix, obslonix);
                        gridptName = sprintf("( %8.4f, %9.4f )", lat_pt, mod(lon_pt+180,360)-180);
                        gridpt_loc = [lat_pt, mod(lon_pt+180,360)-180];
                        gridptLbl = sprintf("%s, %s", gridptName, gridptID);
                        if (strcmp(myhost,"neys")) 
                            lbl2 = ""; %sprintf("%03d_%03d", obslatix,obslonix);
                        else
                            lbl2="";
                        end
%                           For debugging.  illat & illon may not be exact location for gridptLbl.  If using actual lat
%                           & lon, make the code check for "close".
% %                         if (lat_pt == 24.5625 && (any(lon_pt == [-81.8125, 360-81.8125])))
%                           if (illat==4 && illon==7)
%                         if (jloc > 500)
%                             fprintf("here %d %d     %f %f\n", illat, illon, lat_pt, lon_pt-360);
%                             fprintf("\n");
%                         else
%                             continue
%                         end
                        try
                            site_error_flag = false;
                            obs_data = all_obsdata(:,obslatix, obslonix);

                            runLbl  = sprintf('%s %s %s %s : gridptID %s %s', DSP_base.model, DSP_base.ensemble, DSP_base.varname, DSP_base.scenario, gridptID, gridptName);
                            DSP_obs = DSP_obss.update('stnID', gridptID, 'stnName',gridptName, 'fnames',DSP_base.obsnames, 'dspType','obs', "runLbl", runLbl, "gridpt_loc", gridpt_loc);
                            if (DSP_base.do_parfor) 
                                fprintf(               "---worker %3d: ARRM_V2_run_gridded %s grid %5d of %5d, site %5d of %5d:  %s\n",   workerID, DSP_base.llgrid_lbl, jgrid, ngrids, jloc, nlocs, runLbl); 
                            end
                            DSP_base.print_log("\n\n-------worker %3d: ARRM_V2_run_gridded %s grid %5d of %5d, site %5d of %5d:  %s\n\n", workerID, DSP_base.llgrid_lbl, jgrid, ngrids, jloc, nlocs, runLbl);

                            DSP_obs = DSP_obs.update('lats',lat_pt,'lons',lon_pt,'figbase',Parms.figbase);  % set a low threshold for precip so we keep as many obs precips as we can.
                            run_lbl = sprintf('DA_run_%s_%02d_%02d_%04d_%04d_%04d',DSP_base.llgrid_lbl, mdllat_ix, mdllon_ix, obslatix, obslonix, jloc);
                            DA_title = sprintf("DA_obs: %s %s", gridptLbl, run_lbl);
                            RP_obs = RP_base.update('do_anomalies',true, 'do_pdfs',true, 'weights',1); %, 'merge_seed', merge_seed);
%                           DSP_base.print_log("obs  %s: ");
                            
                            [mylatix, mylonix, mdl_weights] = find_closest(lat_pt, lon_pt, DSP_mdl.file_lats, DSP_mdl.file_lons, DSP_mdl.interp_method);        % lonix, latix s/b 1,2...should never have more than 2 in file_lats and file_lons!
                            mwts=mdl_weights(:);
%                           DSP_base.print_log(": latix, lonix:  %4d %4d,  %4d %4d: ", mylatix, mylonix);
%                           DSP_base.print_log("weights:  %6.4f %6.4f %6.4f %6.4f \n", mwts);
                            DSP_base.print_log("latpt, lonpt:  %9.4f %9.4f\n", lat_pt, lon_pt);
                            km_per_degree = 2*pi*6371/360;
                            for jj=1:2
                                for ii = 1:2
                                    mydist = distance(lat_pt, lon_pt, DSP_mdl.file_lats(ii), DSP_mdl.file_lons(jj))*km_per_degree;
                                    DSP_base.print_log("ix: %2d %2d  latlon: %9.4f %9.4f  wt %6.4f  distance %8.3f\n", ii, jj, DSP_mdl.file_lats(ii), DSP_mdl.file_lons(jj), mdl_weights(ii,jj), mydist);
                                end
                            end
                            if (any(mwts<0) || any(mwts>1))
                                DSP_base.warn_log("error:  gridpt;  %s  model outside range 0 to 1:  %8.4f %8.4f %8.4f %8.4f\n", gridptName, mwts);
                                nproblems = nproblems+1;
                                continue
                            elseif (abs(sum(mwts)-1) > 1e-12)   % see if weights are too far off.
                                if (abs(sum(mwts)-1) < 5e-7)    % for nclimgrid, lat & lon are float, not double, and with 1/24 degree gridpoints, sum of weights is often off by about 1e-7. 
                                    mwts = mwts/sum(mwts); % normalize so they sum to 1.
                                else
                                    DSP_base.warn_log("error:  gridpt;  %s  model weights don't sum to 1:  %8.4f %8.4f %8.4f %8.4f  sum %22.20f\n", gridptName, mwts, sum(mwts));
                                    nproblems = nproblems+1;
                                    continue
                                end
                            end
                            if (any(mylonix~=[1,2]) || any(mylatix~=[1,2])) 
                                DSP_base.warn_log("error:  something's wrong, ian!  gridpt:  %s    mylonix:  %d %d     mylatix: %d %d", gridptName, mylonix, mylatix);
                                nproblems = nproblems+1;
                                continue;
                            end
                            
                            DA_obs = ARRM_V2_DisaggregateSignal(obs_data, RP_obs, DSP_obs, DA_title, false, false, false, false, false, false);  % don't do any disaggregation;  we'll do the other steps later. 
                            if (DA_obs.insufficient_data)
                                DSP_base.warn_log("Insufficient Obs data:  %s %s\n", runLbl, gridptName);
                                ninsufficient = ninsufficient+1;
                                continue
                            end
                            
                            if (~Parms.do_compare)
                                obs_fut = [];
                            elseif (isempty(Parms.obsname_fut))
                                obs_fut = obs_data;
                            else
                                [~,obs_fut] = ncdf_read_gridded(Parms.obsname_fut, DSP_base.obsvname, lat_pt, lon_pt, DSP_base.rolling_yrs, calendar, false);  % read & interpolate to lat/lon.
                                if (DA_obs.DP.isPrecipRun)
                                    units = ncdf_getvar_info(Parms.obsname_fut, "pr");
                                    obs_fut = jc_units_conversion(obs_fut, units, "mm");
                                end
                            end
                            
                            
                            [DA_mdl, DA_hist] = single_stream_DAs(mdl_DAs_quadrant, hist_DAs_quadrant, mwts, run_lbl, gridpt_loc);
                            
                            DA_hist.DP = DA_hist.DP.update("figbase", 2);
                            DA_mdl.DP  = DA_mdl.DP.update("figbase", 3);
                            
                            all_files = union(DSP_base.mdlnames, DSP_base.histnames);
                            if (length(all_files) == length(DSP_base.mdlnames))            
                                all_mdl_raw_data = DA_mdl.raw_data;
                            elseif (length(all_files) == length(DSP_base.histnames))
                                all_mdl_raw_data = DA_hist.raw_data;
                            else
                                all_mdl_raw_data = [DA_hist.raw_data;DA_mdl.raw_data];  % This assumes hist & mdl data are contiguous, with no overlap.  This might need fixing, Ian!
                    %             if (DSP_mdl.varname ~= DSP_hist.varname), error("error:  can't read joint data...varnames don't agree:  mdl: %s  hist: %s", DSP_mdl.varname, DSP_hist.varname); end
                    %             all_raw_data = ncdf_read_gridded_data(all_files,  DSP_mdl.varname, lat_pt, lon_pt, DSP_mdl.rolling_yrs, calendar);
                            end
        
                            ok = check_nans(DA_obs,'base_yrs');
                            if (ok)
                                                % run ARRM_V2
                    %                   run ARRM_V2 or ARRM_V2_precip
                                if (DSP_base.isPrecipRun)
                                    
                                        % doesn't actually calc anoms...just moves the raw_data into anoms
                                        % Also sets the pdf_yrlen and calculates the wet-day climatology.
                                    DA_obs.calc_anomalies();
                                    DA_hist.calc_anomalies();
                                    DA_mdl.calc_anomalies();
                                    DA_obs.wet_day_clim_0  = DA_obs.wet_day_clim;
                                    DA_hist.wet_day_clim_0 = DA_hist.wet_day_clim;
                                    DA_mdl.wet_day_clim_0  = DA_mdl.wet_day_clim;      % this will need fixing for regular runs.  OK for GFDL PM.
                                              % trim_excess precip
                                            % NOTE:  there may be 'trace' precip left in the model data.  This will
                                            % occur if there are days in the year where there are more obs precips than
                                            % model precips.  We'll bump these trace values up to the first bin in the
                                            % scale_prcp function.

                                    [DA_obs, DA_mdl, DA_hist, obs_fut, all_mdl_anoms] = trim_prcp(DA_obs, DA_mdl, DA_hist, obs_fut, all_mdl_raw_data, gridptName, Parms.show_figs, gridptLbl, lbl2);
                                    
%                                   obs_fut_orig = obs_fut;     % save for later...
                                    
                                    % and recalculate the climatology after trimming the precip:
                            
                                    DA_obs.calc_precip_climatology("anoms");
                                    DA_hist.calc_precip_climatology("anoms");
                                    DA_mdl.calc_precip_climatology("anoms");
                                   
                                    obs_in   = DA_obs.anoms;        % save original obs input for later comparisons
                                    [bstart, bend, ~        ] = DA_mdl.using_range('rolling_yrs');                    
                                    mdl_in   = DA_mdl.anoms(bstart:bend);        % and save original mdl input for later comparisons
                                    
%                                   prcp_min = DA_obs.DP.prcp_min;

%                                   Parms.show_figs=true;
                                    if (Parms.show_figs)
                                        fignum = 200+2;% *jloc;
                                    else
                                        fignum = [];
                                    end

%                                   if (~isempty(DA_obs.RP.prcp_distrib) && ~any(strcmp(DA_obs.RP.prcp_distrib,["pwr","log"])))
                                    if (DA_obs.RP.isPrecipRun && any(strcmp(DA_obs.RP.prcp_distrib,["generalizedpareto", "inversegaussian", "loglogistic", "lognormal"])))

                                                % calc common binning for all precips.  bin edges are non-uniform, based on log-logistic (prcp_distrib), 
                                                % so precip values are maintained in their original units. 

                                        [~, DA_obs, DA_mdl, DA_hist]  = calc_prcp_binning( DA_obs, DA_mdl, DA_hist, false);
                                    else
%                                         if (~isempty(DA_obs))   % little kludge to avoid matlab warnings!  DA_obs will never be empty here.
%                                             error("error:  you must use loglogistic for the prcp_distrib for now.  icsf 7/23");     % not allowing this version of precip scaling...probably need to remove all code related to this.  icsf 7/23
%                                         end
                                        mnp = [4,1,1;4,1,2;4,1,3; 4,1,4]; 
                                        fidlog = DA_obs.DP.fidlog;
                                        obs_scale_info = get_scaling_info(DA_obs, DA_obs.RP.prcp_distrib, DA_obs.DP.prcp_min, "Obs", gridptLbl, fignum, mnp(1,:), Parms.scaling_offset, fidlog);
                                        if (Parms.show_figs)
                                            get_scaling_info(DA_mdl,  DA_mdl.DP.prcp_min,  "Mdl",  gridptLbl, fignum, mnp(2,:), Parms.scaling_offset, fidlog); 
                                            get_scaling_info(DA_hist, DA_hist.DP.prcp_min, "Hist", gridptLbl, fignum, mnp(3,:), Parms.scaling_offset, fidlog); 
                                        end
                                        mdl_scale_info = get_scaling_info(all_mdl_anoms, DA_obs.RP.prcp_distrib, DA_mdl.DP.prcp_min, "Hist+Mdl", gridptLbl, fignum, mnp(4,:), Parms.scaling_offset, fidlog, lbl2);


                                        DA_mdl.DP.prcp_min = DA_obs.DP.prcp_min;
                                        DA_hist.DP.prcp_min = DA_obs.DP.prcp_min;                                        

                                                % calculate best scaling and scale the data
                                        DA_obs  = scale_prcp(DA_obs,  "anoms", "forward", obs_scale_info);       % scaling info will be stored in DA.prcp_scaling struct by scale_data.
                                        DA_mdl  = scale_prcp(DA_mdl,  "anoms", "forward", mdl_scale_info);
                                        DA_hist = scale_prcp(DA_hist, "anoms", "forward", mdl_scale_info);
                                    end

                                else
                                        % temperature run, not precip.  
                                    obs_in   = DA_obs.orig_anoms;        % save original obs input for later comparisons    % BUT: we haven't disaggregated obs data yet, so this is empt, Ian!  Comparisons won't work yet.
                                    [bstart, bend, ~        ] = DA_mdl.using_range('rolling_yrs');                    
                                    mdl_in   = DA_mdl.orig_anoms(bstart:bend);        % and save original mdl input for later comparisons
                                    
                                            % and now disaggregate and calc binning & climatology for everybody.
                                end
                                
  %                             [~, DA_results, DA_obs_out, DA_hist_out, DA_mdl_out, DA_results_disag] = ARRM_V2(DA_obs, DA_hist, DA_mdl, 'title',runLbl); %,...
                                [~, DA_results, DA_obs_out, DA_hist_out, DA_mdl_out                  ] = ARRM_V2(DA_obs, DA_hist, DA_mdl, 'title',runLbl); %,...
                                                                                              %  'obs_weights', obs_weights, 'hist_weights', hist_weights, ...
                                                                                              %  'mdl_weights', mdl_weights);
                                                                                              
                                                                                            
                            else
                                DSP_base.warn_log('ARRM_V2_run_gridded:  too many NAs processing %s %s run id %s\n', gridptID, gridptName, runLbl);
                                ninsufficient = ninsufficient+1;
                                continue;
                            end

                            if (any(DSP_base.plotFlags))
                                [DSP_base.fignum, mfigs] = plot_ARRM_V2_results(DSP_base.plotFlags, DSP_base.fignum, DA_obs_out, DA_hist_out, DA_mdl_out, DA_results, runLbl);
                                nfigs = nfigs + mfigs;
                            end
                            
                            if (DSP_base.isPrecipRun)
%                               if (isempty(DA_obs.RP.prcp_distrib) || strlength(DA_obs.RP.prcp_distrib) == 0)
                                if (any(strcmp(DA_obs.RP.prcp_distrib, ["pwr","log"])))
                                    DA_results = scale_prcp(DA_results, "mapped_output", 'reverse', DA_results.prcp_scaling);
                                end
                                DA_results.reset_zeros("all");
%                                 DA_results_disag = scale_prcp(DA_results_disag, "mapped_output", 'reverse', DA_results_disag.prcp_scaling);
%                                 DA_results_disag.reset_zeros("all");
                            end


                            % Compare and plot the differences

                            if (Parms.do_compare)%  || Parms.show_figs)
                                mdl_out = DA_results.mapped_output(bstart:bend);

                                cmpr(obs_in, mdl_in, obs_fut, mdl_out, DA_obs.prcp_scaling.prcp_min, 1001+10*jloc, DSP_base, gridptName, Parms.do_compare, Parms.show_figs);
                            end
                            
                            nprocessed = nprocessed+1;
                            if (Parms.save_output)
                                nout = nout+1;
                                nwritten = nwritten+1;
                                
                                
                                save_results(nc_out, nc_zvals, obslatix, obslonix, nwritten, ninsufficient, nlocs, nproblems, DA_results);  % save results to netcdf file
                            end

                            if (~DSP_base.isPrecipRun)
                                report_trends({  DA_obs_out, DA_hist_out, DA_mdl_out, DA_results}, DSP_base, runLbl);
                            end
                            report_counts({  DA_obs_out, DA_hist_out, DA_mdl_out, DA_results}, {'base_yrs','base_yrs','rolling_yrs','rolling_yrs'}, DSP_base, runLbl);
                            report_outliers({DA_obs_out,                          DA_results}, {'base_yrs',                         'rolling_yrs'}, ["obs","map"],DSP_base, runLbl);

            %                     run_ok = true;

                                    % Commented out...no longer keeping
                                    % these arrays of outputs.  Was useful
                                    % for debugging and testing.
                                    % hang on to results if calling code is asking for them.
                                    % only hang on to the ones being asked for, though!
%                             if (Parms.nargout > 1)
%                                 DA_downs_out{nprocessed} = DA_results.clone();
%                                 if (Parms.nargout > 2)
%                                     DA_obss_out{nprocessed}  = DA_obs_out.clone();
%                                     if (Parms.nargout > 3)
%                                         DA_mdls_out{nprocessed} = DA_mdl_out.clone();
%                                         if (Parms.nargout > 4)
%                                             if (isempty(DA_hist_out))
%                                                 DA_hists_out{nprocessed} = [];
%                                             else
%                                                 DA_hists_out{nprocessed}  = DA_hist_out.clone();
%                                             end
%                                             if (Parms.nargout > 5)
%                                                 if (isempty(DA_results_disag))                                    
%                                                     DA_resultss_disag{nprocessed} = [];
%                                                 else
%                                                     DA_resultss_disag{nprocessed} = DA_results_disag.clone();
%                                                 end
%                                                 if (Parms.nargout > 6)
%                                                     if (~isempty(obs_fut))
%                                                         DP_fut = DA_results_disag.DP;
%                                                         RP_fut = DA_results_disag.RP;
%                                                         DA_obs_fut = ARRM_V2_DisaggregateSignal(obs_fut, RP_fut, DP_fut,"GFDL_obs_fut", false, false, false, false, false, false);
%                                                         DA_obs_fut.prcp_scaling = DA_obs.prcp_scaling;  % this scaling might not be perfect, Ian
%                                                         DA_obs_fut = scale_prcp(DA_obs_fut, "raw_data", 'forward', DA_obs_fut.prcp_scaling);
%                                                         DA_obs_fut.run_disaggregation(true, true, true, false, true, true); 
%                                                         DA_obs_fut.mapped_output = obs_fut_orig;
% 
%                                                         DA_obss_fut{nprocessed} = DA_obs_fut.clone();
%                                                     else
%                                                         DA_obss_fut{nprocessed} = [];
%                                                     end
%                                                 end
%                                             end
%                                         end
%                                     end
%                                 end
%                             end                                    

                        catch me
                            site_error_flag = true;
                            nproblems = nproblems + 1;
                            DSP_base.warn_log('error: worker %3d:  ARRM_V2_run_gridded:  on %s %s , run id %s\n', workerID, gridptID, gridptName, runLbl);
                            msgText = getReport(me);
                            DSP_base.warn_log('\t\t error:  %s\t\t%s\n', me.identifier, me.message);
                            DSP_base.warn_log('\nerror caught:\n%s\n',msgText);
                            if (DSP_base.exit_on_error) 
                                rethrow(me);
                            end
                       end
                        elapsed = toc(loop_tic);
                        avg_elapsed = toc(start_tic)/jloc;
                        if (DSP_base.do_parfor)
                            if (Parms.save_output)
                                fprintf(             '---worker %3d: end of run:   %4.1f s (avg %4.1f).        grid %5d of %5d, site %5d of %5d:  %s .  Location [%4d,%4d] in %s\n', workerID, elapsed, avg_elapsed, jgrid, ngrids, jloc, nlocs, runLbl, obslatix, obslonix, outname);
                            else
                                fprintf(             '---worker %3d: end of run:   %4.1f s (avg %4.1f).        grid %5d of %5d, site %5d of %5d:  %s .  %5d sites processed\n',      workerID, elapsed, avg_elapsed, jgrid, ngrids, jloc, nlocs, runLbl, nout);
                            end                            
                        end
                        if (Parms.save_output)
                            DSP_base.print_log('\n-------worker %3d: end of run:   %4.1f s (avg %4.1f).        grid %5d of %5d, site %5d of %5d:  %s .  Location [%4d,%4d] in %s\n', workerID, elapsed, avg_elapsed, jgrid, ngrids, jloc, nlocs, runLbl, obslatix, obslonix, outname);
                        else
                            DSP_base.print_log('\n-------worker %3d: end of run:   %4.1f s (avg %4.1f).        grid %5d of %5d, site %5d of %5d:  %s .  %5d sites_processed\n',      workerID, elapsed, avg_elapsed, jgrid, ngrids, jloc, nlocs, runLbl, nout);
                        end

                        %           if (~keepraw), obs_data{ilat, ilon} = []; end
                        if (~isempty(Parms.figbase)), Parms.figbase = Parms.figbase + 100; end
                        
                    end
                end
            end

                    % make sure we release memory;  shouldn't be needed, but DAs are handles, so may not go away on
                    % their own if matlab isn't quick about garbage collection.
            if (exist('DA_obs', 'var')), clear('DA_obs', 'DA_obs_out'); end
            if (exist('DA_mdl', 'var')), clear('DA_mdl', 'DA_mdl_out'); end
            if (exist('DA_hist','var')), clear('DA_hist','DA_hist_out'); end
        end
    
    catch me        
        if (~site_error_flag && ~grid_error_flag)
            run_gridded_error_flag = true;
            DSP_base.warn_log('\t\terror:  ******** ARRM_V2_run_gridded:  worker %3d:  %s\t\t%s\n', workerID, me.identifier, me.message);
        end
    end

%           myfig=nmdllons*(mdllat_ix-1)+mdllon_ix;
%     myfig = myfig + 100*workerID + randi(99);
%     [mydir,myfn,myext]=fileparts(outname);
%     myfn = strcat(myfn, myext);
%     myttl = sprintf("%s %s", DSP_base.llgrid_lbl, DSP_base.runID);
%     animate_temps(myfn,  DSP_base.varname, myttl, "fignum",myfig,"mode",0,"stepsize",1,"filter",1, "pauselen",.5, "dir",mydir,"year_range",[1979,1979], "matname", "test.mat", "minmax",[])

    DSP_base.print_log('\n');

    nc_total_time = toc(start_tic);
    
    if (nout > 0)
        DSP_base.print_log('output file:  %s, %d gridpoints written\n', outname, nout);
    elseif (Parms.save_output)
        DSP_base.warn_log('ARRM_V2_run_gridded run end.  No output written to %s\n', outname);
    end
    
    error_flag = run_gridded_error_flag || grid_error_flag || site_error_flag;
    
    fids=fopen("all");
    numfiles=numel(fids);
    if (~DSP_base.do_parfor)
        fprintf(       'ARRM_V2_run_gridded: Run completed:  worker %3d: %4s :  %6d sites written, %6d total in file.  total run time:  %10.2f sec.  num open files: %d returning err code %d\n', workerID, runLbl, nwritten, nout, nc_total_time, numfiles, error_flag); 
    end
    DSP_base.print_log('ARRM_V2_run_gridded: Run completed:  worker %3d: %4s :  %6d sites written, %6d total in file.  total run time:  %10.2f sec.  num open files: %d returning err code %d\n', workerID, runLbl, nwritten, nout, nc_total_time, numfiles, error_flag);    
    flag_completion(nc_out, error_flag, nout, nwritten, nlocs, ninsufficient, nproblems, is_appended);
    if (DSP_base.extended)
        flag_completion(nc_zvals, error_flag, nout, nwritten, nlocs, ninsufficient, nproblems, is_appended);
    end
    DSP_base.print_log('log file   :  %s\n', DSP_base.logname);
               
    if (~isempty(DSP_base.fidlog) && DSP_base.fidlog > 2)
        fclose(DSP_base.fidlog);
        if (DSP_base.do_parfor)
            pause(.1);
            type(DSP_base.logname);
            pause(.25);     % give it time to output to screen.  drawnow() won't work here, because I'm using type(...)
        end
    end
    if (exist('me','var'))           
        fprintf("error caught:  worker %3d\n", workerID);
        
%             msgText = getReport(me);
%             DSP_base.warn_log('\nerror caught:\n%s\n',msgText);
        if (DSP_base.exit_on_error), rethrow(me); end         % if exception caught, rethrow exception so wrapper can clean up before abort as well.
    end
    pause(1);
end

function ok = check_nans(DA, yr_type)
%   counts the # of valid points and nas in the data for each of the DAs.

    if (~iscell(DA))
        ok = DA.check_nans(yr_type);
    else
        ndas = numel(DA);

        for k=1:ndas
            if (~isempty(DA{k}))
                ok = DA{k}.check_nans(yr_type);
                if(~ok), return; end
            end
        end    
    end
end

% function grp = make_param_group(grpname, grp_description, obj, excludes, progname)
%         
%     grp = Group(grpname,'Comment', grp_description);
%     grp.putatt('ARRM_V2_Run_Params', progname);
%     props=properties(obj);
%     for i=1:length(props)
%         prop=props{i};
%         if (~ismember(prop,excludes))
%             if (islogical(obj.(props{i})))
%                 grp.putatt(props{i},uint8(obj.(props{i})));
%             elseif (ischars(obj.(props{i})))
%                 grp.putatt(props{i},char(join(obj.(props{i}),'|')));
%             else
%                 grp.putatt(props{i},obj.(props{i}));
%             end
%         end
%     end
% end
% 
% function put_downscaling_params(nc, DSP, progname)
%         
%     nc.putatt('ARRM_V2_Downscaling_Params', progname);
%     props=DSP.properties;
%     excludes=DSP.ncExcludeList;
%     for i=1:length(props)
%         prop=props{i};
%         if (~ismember(prop,excludes))
%             nc.putatt(props{i},RP.(props{i}));
%         end
%     end
% end
% 

% function put_run_params(nc, RP, progname)
%         
%     nc.putatt('ARRM_V2_Run_Params', progname);
%     nc.putatt('run_yrs', 'Years used in this run:  obs_yrs=DSP.obs_yrs(all available observation data), training_yrs=DSP.base_yrs, downscaled_yrs=DSP.model_yrs');
%     nc.putatt('obs_yrs', DSP.base_yrs);
%     nc.putatt('training_yrs', DSP.base_yrs);
%     nc.putatt('downscaled_yrs', DSP.rolling_yrs);
%     nc.putatt('clim_nterms', RP.clim_nterms);
%     nc.putatt('clim_sig_terms', RP.clim_sig_terms);
%     nc.putatt('anom_nterms', RP.anom_nterms);
%     nc.putatt('anom_sig_terms', RP.anom_sig_terms);
%     nc.putatt('thresh', RP.cdf_thresh);
%     nc.putatt('pdf_yrstep', RP.pdf_yrstep);
%     if (~DSP.isPrecipRun)
% %        nc.putatt('mapping_step', RP.mapping_step);
%         nc.putatt('sigrange', RP.sigrange);
%         nc.putatt('nclim_yrs', RP.nclim_yrs);
%         nc.putatt('n_ext_yrs', RP.n_ext_yrs);
%         nc.putatt('pdf_yrs', RP.pdf_yrs);
%         if (~DSP.isStationRun)
%             nc.putatt('trend_order', RP.trend_order);
%             nc.putatt('monthly_valid_count', DSP.monthly_valid_count);
%         end
%     end
% end

function save_results(nc, ncz, latix, lonix, nwritten, ninsufficient, nsites, nproblems, results)
% function save_results(nc, stn_ix, results)  % save results for writing out later

%     nc = ncdf(results.DP.outname);
%     ncid = ncopen_ic(results.DP.outname, 'WRITE');
    mapped = results.mapped_output;  % downscaled results were stored in results.raw_data...
    mapped(isnan(mapped)) = results.DP.FillValue;
    
%   dates = nc.get("Variables/time/data");    
    
    if (results.DP.extended)
        extended = true;
        minz = results.DP.zval_offset;
        maxz = abs(minz);
%         zval_offset = results.DP.zval_offset;
%         zval_scaling = results.DP.zval_scaling;
%         [zpacked, nanval] = ncpack_probs(results.anom_probs_rolling, minz, maxz, 2, zval_offset, zval_scaling);
%         nout_of_range = sum(zpacked == nanval);
%         if (nout_of_range > 0)            
%             results.DP.print_log("warning:  %d probabilities beyond %.1f sigmas\n", nout_of_range, maxz);
%         end
    else
        extended = false;
    end

    nout_of_range = sum(results.anom_probs_rolling < normcdf(minz) | results.anom_probs_rolling > normcdf(maxz));
    if (nout_of_range > 0)            
        results.DP.print_log("warning:  %d probabilities beyond %.4f sigmas\n", nout_of_range, maxz);
    end
    
    start = [1, latix, lonix];    % reminder:  netcdf indexes are 0 based, but matlab's ncwrite uses 1-based indexes.
%     count = [npts, 1];
        
    varname = char(results.DP.varname);
    nc.writevar(varname, mapped, start);
    
    nc.writeatt( "nsites",sprintf("current %6d of %6d sites, %6d insufficient data", nwritten, nsites, ninsufficient));
    
    pct_valid = sum(~isnan(results.mapped_output))/length(results.mapped_output) * 100.0;
    
    results.DP.print_log('pct valid:  loc [%4d, %4d]:  %.2f  \n', latix, lonix, pct_valid)
    if (nproblems > 0)
        nc.writeatt( "problem_sites", int32(nproblems));
    end
    if (extended)
        zvals_name = sprintf("%s_zvals", varname);
        zvals = norminv(results.anom_probs_rolling);
        zvals = max(minz, min(maxz, zvals));    % truncate to zval range.
        ncz.writevar(zvals_name, zvals,      start);
        ncz.writeatt("nsites",sprintf("current %6d of %6d sites, %6d insufficient data", nwritten, nsites, ninsufficient));
        if (nproblems > 0)
            ncz.writeatt("problem_sites", int32(nproblems));
        end
    end
    
end

function flag_completion(nc, error_flag, nout, nwritten, nlocs, ninsufficient, nproblems, is_appended)
%function flag_completion(nc, error_flag, nwritten, ninsufficient, nlocs, is_appended)  % Flags netcdf minifile as completed.

    if (~error_flag)
        if (nwritten == 0)
            nc.writeatt( 'completion_status',int32(2));      % 2 flags "completed, but no data written", usually because no gridpoints have sufficient data to process.  (such as over water livneh gridded obs)
        else
            nc.writeatt( 'completion_status',int32(1));      % 1 flags "completed, with at least some valid data"
        end
    else
        nc.writeatt( 'completion_status',int32(-1));
    end
    if (is_appended)
        nc.writeatt( "nsites",sprintf("%6d sitess in file;  written %6d of %6d sites, %6d insufficient data", nout, nwritten, nlocs, ninsufficient));
    else
        nc.writeatt( "nsites",sprintf("%6d sites in file;    added %6d of %6d sites, %6d insufficient data", nout, nwritten, nlocs, ninsufficient));
    end
    if (nproblems > 0)
        nc.writeatt( "problem_sites", int32(nproblems));
    end
end

function log_header_info(DSP, myhost, myprog, workerID, run_description)

    DSP.print_log( '----------------------- ARRM V2 Gridded Run, %s %s -----------------------host: %s, %s, worker %3d\n%s\n', ...
            datestr(now), DSP.runLbl, myhost, myprog, workerID, run_description);
    
    DSP.print_log( 'obs, hist, model names: \n');
    DSP.print_log(   '        %-12s ', 'obs');
    DSP.print_log( '%s ', DSP.obsnames);
    if (~isempty(DSP.histnames))
        DSP.print_log( '\n        %-12s ', 'hist');
        DSP.print_log( '%s ',DSP.histnames);
    end
    DSP.print_log( '\n        %-12s ', 'model');
    DSP.print_log( '%s ',DSP.mdlnames);
    DSP.print_log( '\n');
    
    DSP.print_log('\nrun log: %s\n\n', DSP.logname);
        
    pause(.1);      % pause long enough to update screen...
end

function [latix, lonix, weights] = find_closest(lat_pt, lon_pt, file_lats, file_lons, interp_method )
% returns index(es) of the four rectangular lat/lon grid point(s) to (lonpt, latpt) in meshgrid(file_lons,file_lats).
% Also returns weights to use for weighted summing data for nearest points.  (These actually aren't necessarily the
% closest, but rather the bounding lat/lon box containing the point (lat_pt, lon_pt).
%
%   Inputs:
%       Lonpt & Latpt are the location to downscale to.  file_lons and file_lats are the actual lats & lons in the file 
%                           covering the area encompassed by the downscaling lats & lons..
%       file_lats, file_lons    lat/lon of the grid covered by the data.  grid is defined by meshgrid(file_lons, file_lats)
%       interp_method       'closest', or 'bilinear', or 'inverse_distance'DSP
%   Outputs:
%       latix, lonix:  indexes of lats & lons closest to lat_pt and lon_pt.
%                           these should bracket the point.
%       weights:        weights to use for merging.  Will be bsed on interp_method.
% 
%   NOTE:   for now, this should always return [1,2] for latix and lonix, because file_lats and file_lons are 2-elements long each;  
%           I've left this as function intact in case I go back to needing to find the closest in a wider range of lats
%           and lons.
%
    if (strncmpi(interp_method,'clo',3))        % use closest gridcell
        [latix, lonix, ~, ~] = closest(lat_pt, lon_pt, file_lats, file_lons, 1, 1);
        weights=1;
    else
        [latix, lonix, latout, lonout] = closest(lat_pt, lon_pt, file_lats, file_lons, 2, 2);     % find closest 4 gridpoints of rectangular bounding box

        weights = calc_weights(lat_pt, lon_pt, latout, lonout, interp_method);
    end
end



function [RP, DSP_base, DSP_obs, DSP_hist, DSP_mdl, all_obsdata, all_histdata, all_mdldata, Parms, calendar] = initParams(varargin)


    % returns an ARRM_V2_RunParams and an ARRM_V2_DownscalingParams object with settings from input arguments.
    % 4 DSP's are returned.  A main one, then one tailored each for Obs, Hist and Model data.
    
    p = inputParser;
    p.StructExpand = true;
    p.KeepUnmatched = true;        
    addParameter(p,"RP",[]);
    addParameter(p,"DSP",[]);
%   addParameter(p,"testsites",strings(0));
    addParameter(p,"figbase",[], @(s) isempty(s) || (isnumeric(s) && numel(s)==1 && s > 0));
    addParameter(p,"save_output",true,@(s) islogical(s) || s==0 || s==1);
    addParameter(p,"latrange", [], @(s) isempty(s) || (isnumeric(s) && length(s)==2 && (all(s >=-90) && all(s <=  90.00001))));
    addParameter(p,"lonrange", [], @(s) isempty(s) || (isnumeric(s) && length(s)==2 && (all(s >=-180) && all(s <= 540.00001))));
    addParameter(p,"hist_yrs", [], @(s) isempty(s) || (isnumeric(s) && length(s)==2 && s(1) >= 1850 && s(2) <= 2100 && s(1)<s(2)));
    addParameter(p,"obs_prcp_min", .1, @(s) isnumeric(s) && s >= 0 && s < 10);
    addParameter(p,"mdl_prcp_min", .1, @(s) isnumeric(s) && s >= 0 && s < 10);
    addParameter(p,"scaling_offset", 1e-5, @(s) isnumeric(s) && s >= 0 && s <=1);
    addParameter(p,"do_internal_plots", [], @(s) isnumeric(s));
    addParameter(p,"do_compare", false, @(s) islogical(s));
    addParameter(p,"show_figs", false, @(s) islogical(s));
    addParameter(p, "obsname_fut", strings(0), @(s) isstring(s) || ischar(s));
    addParameter(p, "run_description", strings(0), @(s) isempty(s) || isstring(s) || ischar(s));
    addParameter(p, "DA_title",strings(0))
%   addParameter(p,"debug_rg",true,@(s) islogical(s) || s==0 || s==1);
   
    parse(p, varargin{:});
    RP                      = p.Results.RP;
    DSP_base                = p.Results.DSP;
    obs_prcp_min            = p.Results.obs_prcp_min;
    mdl_prcp_min            = p.Results.mdl_prcp_min;
    hist_yrs                = p.Results.hist_yrs;
    Parms.figbase           = p.Results.figbase;
    Parms.save_output       = p.Results.save_output;
    Parms.obsname_fut       = p.Results.obsname_fut;
    Parms.do_internal_plots = p.Results.do_internal_plots;
    Parms.scaling_offset    = p.Results.scaling_offset;
    Parms.do_compare        = p.Results.do_compare;
    Parms.show_figs         = p.Results.show_figs;
    Parms.run_description   = p.Results.run_description;
    Parms.DA_title          = p.Results.DA_title;
    
    if (isempty(Parms.figbase)), Parms.figbase = DSP_base.figbase; end   % figbase param may have been grabbed while configuring DSP.
    args = p.Unmatched;
    
    latrange = double(p.Results.latrange);
    lonrange = double(p.Results.lonrange);
%   debug_rg = p.Results.debug_rg;
    
    if (Parms.do_compare || Parms.show_figs)
        fid=fopen("gfdl_test.txt", "w");
        fprintf(fid, "%s\n", datestr(now));
        fclose(fid);
    end

    if (isempty(DSP_base))
        DSP_base = ARRM_V2_DownscalingParams(args);
    else
        DSP_base = DSP_base.update(args);
    end
    args = DSP_base.Unmatched;
    
    if (isempty(DSP_base.progname)), DSP_base.progname = basename(mfilename); end
    DSP_base.ARRM_V2_version = ARRM_V2_version_info(false);      % update the version info and store it.    
        
            % set up list of Run-Parameters that we don't want to save in the output netcdf file
    ncExcludes = {'binstep', 'weights','edges'};    % RP fields which shouldn't be saved in netCDF file.  These change with each location.
    if (~DSP_base.isPrecipRun), ncExcludes{end+1}='prcp_min'; end
    
        % create the RunParams object
    if (isempty(RP))
        RP = ARRM_V2_RunParams(DSP_base.varname, args);
    else
        RP = RP.update(args);
    end
    Unmatched = RP.Unmatched;
    
    for i=1:length(ncExcludes)
        prop = ncExcludes{i};
        if (~ismember(prop, RP.ncExcludeList))
            RP.ncExcludeList{end+1}=prop;
        end
    end
    
        % defaults for Precip runs.
    if (DSP_base.isPrecipRun)
        if (~any(strcmp(varargin,"pdf_yrlen"))), RP.pdf_yrlen = 128; end
        if (~any(strcmp(varargin, "anom_nterms")))
%             RP.cdf_thresh = .001;        %  cdf not valid beyond .001.  This should be based on the number of precip events...
%             RP.outlier_thresh = .0228;   % flagged as outliers beyond 97.72% (2 sigma) and far outliers beyond 98%  These will be mapped by a simple scaling This should be based on the number of precip events... 
%             RP.far_outlier_thresh = .01; % flag far outliers beyond 99th percentile            
%             RP.far_outlier_anchor_pt = .1; % and use 1-sigma point for rescaling far outliers.
%             RP.cdf_thresh = .0005;        %  cdf not valid beyond .001.  This should be based on the number of precip events...
%             RP.outlier_thresh = [0,.99];   % flagged as outliers beyond 97.72% (2 sigma) and far outliers beyond 98%  These will be mapped by a simple scaling This should be based on the number of precip events... 
%             RP.far_outlier_thresh = [0,.999]; % flag far outliers beyond 99th percentile            
%             RP.far_outlier_anchor_pt = .1; % and use 1-sigma point for rescaling far outliers.
% xxxxx           DSP_base.cdf_append_pts = 0;
%             RP.na_thresh = 1e-12;
%           RP.pdf_map_method = "clim";
        end
    end
            % set up the log file.
    if (DSP_base.do_log)
        [dr,nam,~] = fileparts(DSP_base.outname);
        if (strlength(dr)==0), dr = DSP_base.outdir; end
        if (isempty(DSP_base.logname))
            DSP_base.logname = fullfile(dr,sprintf("run_gridded.%s.log", nam));
        end
        if (DSP_base.do_log == 1 || DSP_base.do_log == 2)
            DSP_base.fidlog = fopen(DSP_base.logname,'a');
        elseif (DSP_base.do_log == 3)
            DSP_base.fidlog = fopen(DSP_base.logname,'w');
        end
    else
        if (ispc())     % no log output.
            DSP_base.fidlog = fopen('NUL:');        % windows equivalent of Unix's /dev/null
        else
            DSP_base.fidlog = fopen('/dev/null');   % Unix pseudo-file for writing to when you don't want to keep the output.
        end
    end
    
    RP.do_pdfs = false;

    if (~isempty(Unmatched) && ~isempty(fieldnames(Unmatched)))
        fprintf(2,"Error:  unmatched input parameters:\n"); 
        fnms = fieldnames(Unmatched);
        for i=1:length(fnms)
            fprintf(2, "\t%-12s %s\n", fnms{i}, string(Unmatched.(fnms{i})));
        end
        error("unexpected input parameters");        
    end
    
    DSP_base = DSP_base.replace_keywords();   % parses inputs and replaces any keywords in the inputs with appropriate values, such as
                                    % replacing MODEL, SCENARIO and ENSEMBLE with actual model, scenario and ensemble. 
        
            % get date ranges from input files.
            
            
    obs_file_yrs = ncdf_get_date_range(DSP_base.obsnames);
    mdl_file_yrs = ncdf_get_date_range(DSP_base.mdlnames);

    hist_file_yrs = ncdf_get_date_range(DSP_base.histnames);
    
            % get lat/lon ranges.
            
    [mdl_lats, mdl_lons] = ncdf_get_latlons(DSP_base.mdlnames);
    mdl_lons = mod(mdl_lons,360);

    [obs_lats, obs_lons] = ncdf_get_latlons(DSP_base.obsnames);
    obs_lons = mod(obs_lons,360);
    
    if (~isempty(latrange))
        [~, obs_lats] = lat_region(latrange, obs_lats, [-1,-2], true);
    end
    if (~isempty(lonrange))
        [~, obs_lons] = lon_region(lonrange, obs_lons, [-1,-2], true);
    end
    
    
    [hist_lats, hist_lons] = ncdf_get_latlons(DSP_base.histnames);
    hist_lons = mod(hist_lons,360);

    if (~all(hist_lats == mdl_lats) || ~all(hist_lons == mdl_lons)), error("error:  hist lats/lons don't match mdl lats/lons"); end
    
                % If run info didn't specify lat/lon range, then use entire obs range .
    if (isempty(DSP_base.lats)), DSP_base.lats = obs_lats; end
    if (isempty(DSP_base.lons)), DSP_base.lons = obs_lons; end
    
            % get range of lat & lon values for pulling Hist & Mdl data.
    Parms.latrange = [min(DSP_base.lats),max(DSP_base.lats)];
    if (DSP_base.lons(1) <= DSP_base.lons(end))
        Parms.lonrange = [min(DSP_base.lons),max(DSP_base.lons)];
    else
        Parms.lonrange = [Parms.lonrange(1), 359.99999; 0, Parms.lonrange(end)];
    end
    if (Parms.latrange(1)==Parms.latrange(2) && Parms.lonrange(1)==Parms.lonrange(2))
        Parms.latrange=Parms.latrange(1);
        Parms.lonrange=Parms.lonrange(1);
    end
    
%     fix the years, Ian!
    if (isempty(DSP_base.base_yrs)),     error("error:  base_yrs    not set"); end
    if (isempty(DSP_base.rolling_yrs)),  error("error:  rolling_yrs not set"); end
                
        % create obs, hist & model-specific DSP's
            % Create the obs DSP
            % read the Obs ncdf metadata, and set all the date ranges for the Obs data.
    DSP_obs = DSP_base;                        
    nc =  ncdf(DSP_obs.obsnames{1}); 

    DSP_obs.prcp_min = obs_prcp_min;       % update the min. prcp value if passed in, or with default if not
    DSP_obs.file_lats = obs_lats;
    DSP_obs.file_lons = obs_lons;

    DSP_obs.interp_method = DSP_base.obs_interp_method;
    DSP_obs = DSP_obs.update('dspType','obs','rolling_yrs',[]);
    DSP_obs = DSP_obs.set_yr_limits(min_yr_range(DSP_base.base_yrs, obs_file_yrs));          % truncates all yr_ranges to the range of data available in the file.
 
    DSP_obs.base_yrs = min_yr_range(DSP_obs.base_yrs, obs_file_yrs);
    ncvar = nc.get(DSP_obs.obsvname);
    DSP_obs.units = ncvar.getattvalue('units');
    try
        DSP_obs.varlongname = ncvar.getattvalue('long_name');
    catch
        DSP_obs.varlongname = '';
    end
%   DSP_obs.data_final_yrs = min_yr_range(DSP_obs.data_final_yrs, DSP_obs.rolling_yrs);         % shouldn't need this, ian!
%   if (isempty(DSP_obs.rolling_yrs)), DSP_obs.rolling_yrs = obs_file_yrs; end

    DSP_obs = DSP_obs.set_yr_limits(obs_file_yrs);          % truncates all yr_ranges to the range of data available in the file.

    DSP_obs.dspType = "Observations"; 
    
        % Create the model DSP
            % read the Model ncdf metadata, and set all the date ranges for the Model data.

    DSP_mdl = DSP_base;
    DSP_mdl.prcp_min = mdl_prcp_min;       % update the min. prcp value if passed in, or with default if not
    
    if (~isempty(DSP_mdl.fignum)), DSP_mdl.fignum = DSP_mdl.fignum + 25; end  % so we don't tromp on the figures for obs data
    nc =  ncdf(DSP_mdl.mdlnames{1}); 
    
        % make sure model lats & lons bracket the range needed.
    [~,~,~, qlats,qlons] = latlon_region(DSP_mdl.lats, DSP_mdl.lons, mdl_lats, mdl_lons,[2,2]);      % get list of lats & lons needed to bound the run's lats & lons.
    DSP_mdl.lats = qlats;
    DSP_mdl.lons = qlons;
    
            % get the variable units
    ncvar = nc.get(DSP_mdl.varname);
    DSP_mdl.units = ncvar.getattvalue('units');
    try
        DSP_mdl.varlongname = ncvar.getattvalue('long_name');
    catch
        DSP_mdl.varlongname = '';
    end
    
            % set up the rolling years
    DSP_mdl.data_final_yrs = min_yr_range(DSP_mdl.data_final_yrs, DSP_mdl.rolling_yrs);
    DSP_mdl.base_yrs = [];      % why was this commented out, Ian?
    if (isempty(DSP_mdl.rolling_yrs)), DSP_mdl.rolling_yrs = mdl_file_yrs; end

    DSP_mdl = DSP_mdl.set_yr_limits(mdl_file_yrs);          % truncates all yr_ranges to the range of data available in the file.
            % limit the rolling years calculations to only the range for final years.
%     if (~isempty(DSP_base.data_final_yrs))
%         DSP_base.rolling_yrs = [nanmax(DSP_base.data_final_yrs(1),DSP_base.rolling_yrs(1)), nanmin([DSP_base.data_final_yrs(2),DSP_base.rolling_yrs(2)])];
%     end
    DSP_base.rolling_steps = DSP_base.calc_rolling_steps(RP.pdf_yrstep);    
    DSP_mdl.rolling_steps  =  DSP_mdl.calc_rolling_steps(RP.pdf_yrstep);  
    DSP_mdl = DSP_mdl.update('dspType','model', 'stninfo',[], 'stnID',[]);
    DSP_mdl.trend_yrs = DSP_mdl.rolling_yrs;    % Also use rolling years for trend if separate hist.  NOTE:  we'll need to add back the difference between hist's avg_base and model's avg_base later.

        % Create the hist DSP
            % normally, hist will overlap with downscaled data and be pulled from the same files
            
    DSP_hist = DSP_base;
    DSP_hist.prcp_min = mdl_prcp_min;       % update the min. prcp value if passed in, or with default if not
    nc =  ncdf(DSP_hist.histnames{1}); 
    if (~isempty(DSP_hist.fignum)), DSP_hist.fignum = DSP_hist.fignum + 50; end  % so we don't tromp on the figures for obs data
    if (~isempty(hist_yrs))
        DSP_hist.base_yrs = hist_yrs;
    end
    DSP_hist.base_yrs = min_yr_range(DSP_hist.base_yrs, hist_file_yrs);
    DSP_hist.rolling_steps = [];                % so we don't calculate any rolling stuff for separate hist.
    % make sure our hist and model lats & lons bracket the range needed.
    [~,~,~, qlats,qlons] = latlon_region(DSP_hist.lats, DSP_hist.lons, hist_lats, hist_lons,[2,2]);      % get list of lats & lons needed to bound the run's lats & lons.
    DSP_hist.lats = qlats;
    DSP_hist.lons = qlons;

    ncvar = nc.get(DSP_hist.varname);
    DSP_hist.units = ncvar.getattvalue('units');
    try
        DSP_hist.varlongname = ncvar.getattvalue('long_name');
    catch
        DSP_hist.varlongname = '';
    end
%   DSP_obs = DSP_obs.set_yr_limits(min_yr_range(DSP_base.base_yrs, obs_file_yrs));          % truncates all yr_ranges to the range of data available in the file.
    DSP_hist = DSP_hist.set_yr_limits(min_yr_range(DSP_base.base_yrs, hist_file_yrs));          % truncates all yr_ranges to the range of data available in the file.        
    DSP_hist.rolling_yrs = [];
    DSP_hist.data_final_yrs = [];
    DSP_hist = DSP_hist.update('dspType','hist', 'stninfo',[], 'stnID', []);

        % update changed fields in main DSP with changed info.
%     DSP.lats            = DSP_obs.lats;       % already done earlier.
%     DSP.lons            = DSP_obs.lons;
    DSP_base.units           = DSP_obs.units;
%     DSP.stnID           = DSP_obs.stnID;
%     DSP.stninfo         = DSP_obs.stninfo;

%   DSP_base.trend_yrs       = max_yr_range(DSP_obs.trend_yrs, DSP_mdl.trend_yrs);
%   DSP_base.base_yrs        = DSP_obs.base_yrs;
%   DSP_base.rolling_yrs     = DSP_mdl.rolling_yrs;
    DSP_base.data_final_yrs  = DSP_mdl.data_final_yrs;
    DSP_base.rolling_steps   = DSP_mdl.rolling_steps;
%   DSP_base.data_yrs        = max_yr_range(obs_file_yrs, mdl_file_yrs, hist_file_yrs);
    
    DSP_base.varlongname     = DSP_mdl.varlongname;
    
%     DSP_base.fnames = basename([DSP_base.obsnames, DSP_base.histnames, DSP_base.mdlnames]);
%     DSP_base.datadir = dirname([DSP_base.obsnames, DSP_base.histnames, DSP_base.mdlnames]);   
    DSP_base.dspType = 'global';
    
    DSP_base = DSP_base.finalize();
    DSP_obs  = DSP_obs.finalize();
    DSP_mdl  = DSP_mdl.finalize();
    DSP_hist = DSP_hist.finalize();
    
    DSP_base.check_files_exist();     % makes sure all files exist.  Throws error if not.
    
    if (DSP_base.yrlen == 365)
        calendar = '365-day';
    elseif (DSP_base.yrlen == 365.25)
        calendar = 'standard';
    else
        calendar = '360-day';
    end
        
    nc_obs = ncdf_read_files(     DSP_obs.obsnames, DSP_obs.obsvname, Parms.latrange, Parms.lonrange, DSP_obs.base_yrs, calendar, 0, DSP_obs.obsvname, false, [-1,-1]);     % 0:  no random selection.   [-1, -1]: lon, lon:  "up to range";
%   nc_obs = ncdf_read_files(     DSP_obs.obsnames, DSP_obs.obsvname, Parms.latrange, Parms.lonrange, DSP_obs.base_yrs, calendar, 2, DSP_obs.obsvname, false, [-1,-1]);     % 2:  will seed RNG with 2.   [-1, -1]: lon, lon:  "up to range";
    all_obsdata = nc_obs.getvardata(DSP_obs.obsvname);
    DSP_obs.file_lats = nc_obs.getvardata('lat');
    DSP_obs.file_lons = nc_obs.getvardata('lon');
    
    if (DSP_base.isPrecipRun)
        units = nc_obs.getattvalue(sprintf("/Variables/%s/Attributes/units", DSP_base.obsvname));
        all_obsdata = jc_units_conversion(all_obsdata, units,"mm");
        all_obsdata = max(0, all_obsdata,"includenan");     % includenan leaves nan's in place.  omitnan replaces them with max of remaining stuff.
        DSP_obs.units="mm";
    end
    
    nc_hist = ncdf_read_files(     DSP_hist.histnames, DSP_hist.varname, Parms.latrange, Parms.lonrange, DSP_hist.base_yrs, calendar, 0, DSP_hist.varname, false, [1,1] );     % 0:  no random selection.  lat, lon:  "at least range"
%   nc_hist = ncdf_read_files(     DSP_hist.histnames, DSP_hist.varname, Parms.latrange, Parms.lonrange, DSP_hist.base_yrs, calendar, 2, DSP_hist.varname, false, [1,1] );     % 2:  will seed RNG with 2.  lat, lon:  "at least range"
    all_histdata = nc_hist.getvardata(DSP_hist.varname);
    DSP_hist.file_lats = nc_hist.getvardata('lat');
    DSP_hist.file_lons = nc_hist.getvardata('lon');
    if (DSP_base.isPrecipRun)
        units = nc_hist.getattvalue(sprintf("/Variables/%s/Attributes/units", DSP_hist.varname));
        all_histdata = jc_units_conversion(all_histdata, units,"mm");
        all_histdata = max(0, all_histdata,"includenan");
        DSP_hist.units="mm";
    end

    nc_mdl = ncdf_read_files(     DSP_mdl.mdlnames, DSP_mdl.varname, Parms.latrange, Parms.lonrange, DSP_mdl.rolling_yrs, calendar, 0, DSP_mdl.varname, false, [1,1]);     % 0:  no random selection.  lat, lon:  "at least range"
%   nc_mdl = ncdf_read_files(     DSP_mdl.mdlnames, DSP_mdl.varname, Parms.latrange, Parms.lonrange, DSP_mdl.rolling_yrs, calendar, 2, DSP_mdl.varname, false, [1,1]);     % 2:  will seed RNG with 2.  lat, lon:  "at least range"
    all_mdldata = nc_mdl.getvardata(DSP_mdl.varname);
    DSP_mdl.file_lats = nc_mdl.getvardata('lat');
    DSP_mdl.file_lons = nc_mdl.getvardata('lon');
    check_hist_mdl_lat_lons(DSP_hist, DSP_mdl);     % make sure lats, lons match!
    if (DSP_base.isPrecipRun)
        units = nc_mdl.getattvalue(sprintf("/Variables/%s/Attributes/units", DSP_mdl.varname));
        all_mdldata = jc_units_conversion(all_mdldata, units,"mm");
        all_mdldata = max(0, all_mdldata,"includenan");
        DSP_mdl.units="mm";
    end
end

function check_hist_mdl_lat_lons(DSP_hist, DSP_mdl)
% throws an error if the list of lats & lons aren't identical for hist and model files.

    ok = true;
    if (length(DSP_hist.lats) ~= length(DSP_mdl.lats) || length(DSP_hist.lons) ~= length(DSP_mdl.lons)), ok = false; end
    
    if (ok && (~all(DSP_hist.lats == DSP_mdl.lats) || ~all(DSP_hist.lons == DSP_mdl.lons))), ok = false; end
    
    if (~ok), error("Error:  model and hist lats and lons don't match"); end

end

function nc = ncdf_init_grid(RP, DSP, progname, calendar, do_zvals)

    lats = DSP.lats;
    lons = DSP.lons;
    if (~isempty_s(DSP.globatts_name)), DSP.NC_atts = read_global_attributes(DSP.globatts_name, DSP.NC_atts); end
    
    DSP = DSP.replace_keywords();       % update all keyword strings in DSP's fields.
    
    
    if (isAbsolute(DSP.outname))
        out_ncName = DSP.outname;
    else        
        out_ncName = fullfile(DSP,outdir, DSP.outname);
    end
    
    nc = ncdf('','Filename',out_ncName, 'Format','netcdf4','create_ok',true);
    
    put_globals(nc, DSP, progname, do_zvals);  % inserts title, long title, comments, description, date_range, creation date, references, history.  creator  
        
    nc.putatt('completion_status',int32(0));  % 0:  not completed.  1:  completed;  2:  no valid data in gridbox.  -1:  errors encountered.
    if (DSP.do_parfor)
        nc.putatt('single_gridbox',int32(0));
    else
        nc.putatt('single_gridbox',int32(1));
    end   
    
    base_year = DSP.rolling_yrs(1);
         
%   [timedim, latdim, londim, ctimedim] = put_dims(nc, lats, lons, base_year, DSP.model_yrs);   % put_dims(nc,time, meta.model.lats, meta.model.lons, base_year);
    [timedim, latdim, londim] = put_dims(nc, lats, lons, base_year, DSP.data_final_yrs, calendar);   % put_dims(nc,time, meta.model.lats, meta.model.lons, base_year);
 
    dims = {timedim, latdim, londim};
    
    if (~do_zvals)
        attributes = {Attribute('long_name', DSP.varlongname), Attribute('units',DSP.units)};
        add_var(nc, DSP.varname, attributes, dims, single(DSP.FillValue), []);    
        nc.putatt('data_variables',DSP.varname);
        DSP.ncExcludeList = [DSP.ncExcludeList, "zval_offset","zval_scaling"];
    else
        zvals_name = sprintf("%s_zvals", DSP.varname);
        zvalVar = Variable(zvals_name,'Dimensions', [timedim, latdim, londim], 'Datatype','uint16','FillValue',65535);
        zvalVar.putatt('add_offset',DSP.zval_offset);
        zvalVar.putatt('scale_factor', DSP.zval_scaling);
        zvalVar.putatt('description', 'cumulative probability, stored as a z-score-equivalent.  To get probability of a datapoint from its zval in matlab:  prob(zval) = normcdf(zval); or in R:  prob(zval) = normalcdf(0, zval)');
        nc.putvar(zvalVar);
    end
    

    rungrp = make_param_group('RunParams', 'ARRM V2 Run Parameters', RP, RP.ncExcludeList, progname);
    nc.putgrp(rungrp);
    downgrp = make_param_group('DownscalingParams','ARRM V2 Downscaling Parameters', DSP, DSP.ncExcludeList, progname);
    nc.putgrp(downgrp);


%   add_var(nc, 'completed',{Attribute('description','per-site completion flags')}, {latdim, londim}, 0, zeros(length(DSP.lats), length(DSP.lons)));   
    
    src_grp = make_source_group(DSP.obsnames, DSP.histnames, DSP.mdlnames, DSP.saveMetaData);
    nc.putgrp(src_grp);
    
    nc.writeschema(true, 'Dimensions',true,'Variables',false);
end

function put_globals(nc, DSP, progname, do_zvals)
% Inserts title, long title, comments, description, date_range, creation date, references, history.
% Also, progname, data sources.

    filenames = join([string(DSP.obsnames),string(DSP.histnames),string(DSP.mdlnames)]); 
    
    NC_atts = DSP.NC_atts;
    
    fields = fieldnames(NC_atts);
    for i=1:length(fields)
        field = fields{i};
        nc.putatt(field,NC_atts.(field));
    end
    nc.putatt('progname',progname);
    nc.putatt('data_source',sprintf('data sources: %s', filenames));
    if (do_zvals)
        nc.putatt("ZVALS", sprintf("This file contains probability values for %s, encoded as z-vals", DSP.varname));
    end
    
    [~,uname] = getusername();
    nc.putatt('creator',uname);
    
    
%   nc.putatt('interpolation_method',DSP.interp_method)
    
%     if (DSP.isStationRun)
%         nc.putatt('Station_ID', DSP.stnID);
%         nc.putatt('Station_Name', DSP.stnName);
%         nc.putatt('Station_Lat_Lon', [DSP.lats(1),DSP.lons(1)]);
%         nc.putatt('Station_Elevation', DSP.elev);
%     end
    
end

    
% function [timedim, latdim, londim, ctimedim, yrstepsdim, probsdim, monthsdim] =  put_dims(nc, lats, lons, base_year, model_yrs, calendar, ctime, yrsteps, probs, months)
function   [timedim, latdim, londim] =  put_dims(nc, lats, lons, base_year, model_yrs, calendar)
    
    time_model = get_timevals_cal(model_yrs, base_year, calendar);
    time_units = sprintf('days since %4d-01-01,00:00:00', base_year); 
    
    timedim = Dimension('time',length(time_model), false);
    timevar = Variable('time',time_model,'Dimensions',timedim,'standard_name','time','long_name','Time of measurements', 'calendar',calendar,'units',time_units,'Datatype','double');
    nc.putdim(timedim);
    nc.putvar(timevar);
    
    latdim = Dimension('latitude',length(lats),false);
    latvar = Variable('latitude',lats,'Dimensions',latdim,'standard_name','latitude','long_name','latitude','units','degrees_north','Datatype','double');
    nc.putdim(latdim);
    nc.putvar(latvar);
    
    londim = Dimension('longitude',length(lons),false);
    lonvar = Variable('longitude',lons,'Dimensions',londim,'standard_name','longitude','long_name','longitude','units','degrees_east','Datatype','double');
    nc.putdim(londim);
    nc.putvar(lonvar);    
    
%     if (exist('ctime','var') && ~isempty(ctime))
%         ctimedim = Dimension('ctime',length(ctime), false);
%         ctimevar = Variable('ctime',ctime,'Dimensions',ctimedim,'standard_name','ctime','long_name','Basic Climatology Time', 'calendar','365_day','units','day of year','Datatype','double');
%         nc.putdim(ctimedim);
%         nc.putvar(ctimevar);    
%     else
%         ctimedim=[];
%     end
%     if (exist('yrsteps','var') && ~isempty(yrsteps))
%         yrstepsdim = Dimension('yrsteps',length(yrsteps), false);
%         yrstepsvar = Variable('yrsteps',yrsteps,'Dimensions',yrstepsdim,'standard_name','yrsteps','long_name','model year distribution years','units','years since 1900-00-00','Datatype','double');
%         nc.putdim(yrstepsdim);
%         nc.putvar(yrstepsvar);    
%     else
%         yrstepsdim=[];
%     end
%     if (exist('probs','var') && ~isempty(probs))
%         probsdim = Dimension('probs',length(probs), false);
%         probsvar = Variable('probs',probs,'Dimensions',probsdim,'standard_name','probs','long_name','probability values for problines','Datatype','double');
%         nc.putdim(probsdim);
%         nc.putvar(probsvar);    
%     else
%         probsdim=[];
%     end    
%     if (exist('months','var') && ~isempty(months))
%         monthsdim = Dimension('months',length(months), false);
%         monthsvar = Variable('months',months,'Dimensions',monthsdim,'standard_name','months','long_name','month of year (1-12)','units','month','Datatype','double');
%         nc.putdim(monthsdim);
%         nc.putvar(monthsvar);    
%     else
%         monthsdim=[];
%     end
end

function grp = make_param_group(grpname, grp_description, obj, excludes, progname)
        
    grp = Group(grpname,'Comment', grp_description);
    grp.putatt('ARRM_V2_Run_Params', progname);
    props=properties(obj);
    for i=1:length(props)
        prop=props{i};
        if (~ismember(prop,excludes))
            if (islogical(obj.(props{i})))
                grp.putatt(props{i},uint8(obj.(props{i})));
            elseif (ischars(obj.(props{i})))
                grp.putatt(props{i},char(join(obj.(props{i}),'|')));
            else
                grp.putatt(props{i},obj.(props{i}));
            end
        end
    end
end

function add_var(ncdf_out, varName, attributes, dimensions, FillValue, vdata)

    if (strcmp(ncdf_out.Format,'64bit'))
        FillValue_lbl = '_FillValue';
    else
        FillValue_lbl = 'FillValue';
    end
    if (exist('vdata','var'))
        if (isempty(vdata)), vdata = make_vdata(dimensions); end
        myvar = Variable(varName, vdata, 'Attributes',attributes, 'Dimensions',dimensions,FillValue_lbl,FillValue, 'Datatype',class(vdata)) ;
    else       
        myvar = Variable(varName,        'Attributes',attributes, 'Dimensions',dimensions,FillValue_lbl,FillValue, 'Datatype','single') ;
    end
    ncdf_out.putvar(myvar);
end

function vdata = make_vdata(dimensions)
    ndims = length(dimensions);
    siz = zeros(1,ndims);
    for i=1:ndims
        siz(i) = dimensions{i}.Length;
    end
    vdata = single(nan(siz));
end



function    cmpr(obs_in, mdl_in, obs_fut, mdl_out, prcp_min, fignum, DSP_base, siteLbl, do_compare, show_figs)

    obs_in = double(obs_in);
    mdl_in = double(mdl_in);
    obs_fut = double(obs_fut);
    mdl_out = double(mdl_out);
    
    nyrs_obs = length(obs_in)/365;
    nyrs_mdl = length(mdl_out)/365;

    nyrs_obs_fut = length(obs_fut)/365;
    nyrs_mdl_in = length(mdl_in)/365;

%   maxval = nanmax([obs_in; mdl_in; obs_fut; mdl_out]);
    maxval =    max([obs_in; mdl_in; obs_fut; mdl_out]);  % in 2020b & later, max defaults to 'omitnan'
    edges = .1:1:maxval+1;
    step = edges(2)-edges(1);  
    bins = edges(1:end-1) + step/2;
    h_in  = histcounts(obs_in, edges);
    h_mdl = histcounts(mdl_out, edges);
    hobs_fut = histcounts(obs_fut, edges);
    hmdl_out = histcounts(mdl_out, edges);
    

    if (show_figs)
                        prcp_max = double(max(obs_fut));
                        s_obs_fut = sort(obs_fut(obs_fut>0));
                        numo = length(s_obs_fut);
                        s_mdl_out = sort(mdl_out(mdl_out>0));
                        numm = length(s_mdl_out);
                        ixout = round(linspace(1,numm, numo));
                        s_mdl_in  = sort(mdl_in(mdl_in>0));
                        numin = length(s_mdl_in);
                        ixin = round(linspace(1,numin, numo));
                        
                        ix50 = round(.5*numo);
                        ix90 = round(.9*numo);
                        ix95 = round(.95*numo);
                        ix99 = round(.99*numo);
    
                        diffs = s_mdl_out(ixout) - s_obs_fut;
                        diffs0 = diffs;

                        err_min = min(diffs);
                        err_max = max(diffs);
                        err_min_99 = min(diffs(ix50:ix99));
                        err_max_99 = max(diffs(ix50:ix99));
                        pcts=[.5,.9,.95,.99, .995, .999];
                        mymin = min(-12, .85*err_min);
                        mymax = max( 12, .85*err_max);
                        for i=1:length(pcts)
                            ix = round(pcts(i)*numo);
                            pcts2(i) = s_obs_fut(ix); %#ok<AGROW>
                        end

                        h = figure(fignum);
                        h.Position = [1200,5,1100,1330];
                        
                        subplot(5,1,1);
                        plot([0,s_obs_fut(end)],[0,0], "r-", "linewidth",3);
                        hold on; 
                        scatter(s_obs_fut, diffs, 15, "b","filled");
                        for i=1:length(pcts2)
                            plot([pcts2(i),pcts2(i)],[mymin, mymax],"r--","linewidth",2);
                            ypos = mymin + .85*(mymax-mymin) * i/length(pcts2);
                            text(pcts2(i)+.005*prcp_max, ypos ,sprintf(" %.1f %%", pcts(i)*100),"color","k");
                        end
                        hold off;
                        axis("tight")
                        title(sprintf("Error, mm, (Downscaled - Obs), %s", siteLbl),"interpreter","none");
                        grid on;
                        yl=ylim();
                        yl = [min(-20,yl(1)),max(20, yl(2))];
                        ylim(yl);
                        xlabel("precip, mm");
                        ylabel("error, mm");
                        if (range(yl) < 50)
                            yticks(-100:5:100);
                        end
                        
                        
                            % limit the errors to 15 mm.
                        in  = diffs > -20 & diffs < 20;
                        out = ~in;
                        nout = sum(out);
                        diffs = min(40, max(-40, diffs));
                        yl = [min(-20, err_min_99), max(20, err_max_99)];
                        
                        subplot(5,1,2);
                        x = linspace(0,100,numo);
                        plot([0,100],[0,0], "r-", "linewidth",3);
                        hold on;
                        h1 = scatter(x(in),diffs(in), 15, "b","filled");
                        if (nout > 0)
                            h2 = scatter(x(out),diffs(out), 15, "m","filled");
                            legend([h1,h2], "diff","diff (off-scale)");
                        end
                        hold off;
                        title("error, mm, by percentile","interpreter","none");
                        xlabel("percentile")
                        ylabel("error, mm");
                        if (range(yl) < 50)
                            yticks(-100:5:100);
                        end
                        ylim(yl);
                        grid on;
                        
                        wdobs  = sum(obs_in>0)/nyrs_obs;
                        wdmdl  = sum(mdl_out>0)/nyrs_mdl;
                        wdfobs = sum(obs_fut>0)/nyrs_mdl;
                        
                        subplot(5,1,3);
                        plot([0,100],[0,0], "r-", "linewidth",3);
                        hold on;
                        h1 = scatter(x(in),diffs(in), 15, "b","filled");
                        if (nout > 0)
                            h2 = scatter(x(out),diffs(out), 15, "m","filled");
                            legend([h1,h2], "diff","diff (off-scale)");
                        end
                        hold off;
                        title(sprintf("error, mm, by percentile, top 95-100%%   wet days: obs %5.1f    mdl: %5.1f    fut_obs: %5.1f",wdobs, wdmdl, wdfobs), "interpreter","none");
                        xlabel("percentile")
                        ylabel("error, mm");
                        xlim([95,100]);
                        ylim(yl);
                        if (range(yl) < 50)
                            yticks(-100:5:100);
                        end
                        grid on;

                        pct_err = diffs0 ./ s_obs_fut * 100;
                        in = pct_err > -50 & pct_err < 50;
                        out = ~in;
                        nout = sum(out);
                        
                        pct5090 = mean(pct_err(ix50:ix90));
                        pct5095 = mean(pct_err(ix50:ix95));
                        pct9095 = mean(pct_err(ix90:ix95));
                        
                        subplot(5,1,4);
                        plot([0,100],[0,0], "r-", "linewidth",3);
                        hold on;
                        h1 = scatter(x(in),pct_err(in), 15, "b","filled");
                        if (nout > 0)
                            h2 = scatter(x(out),pct_err(out), 15, "b","filled");
                            legend([h1,h2], "pct diff","pct diff (off-scale)");
                        end
                        
                        hold off;
                        title(sprintf("Percent error, 50-100%% percentile   50-90: %.1f    50-95: %.1f    90-95: %.1f", pct5090, pct5095, pct9095),"interpreter", "none");
                        xlabel("percentile")
                        ylabel("percent error");
                        xlim([50,100]);
%                         yl = ylim();
%                         ylim([min(-20,yl(1)),max(20,yl(2))]);
                        ylim(yl);
                        if (range(yl) < 50)
                            yticks(-100:5:100);
                        end 
                        grid on
                        
                        cum_obs_fut = cumsum(s_obs_fut)/nyrs_obs_fut;
                        cum_mdl_out = cumsum(s_mdl_out)/nyrs_mdl;
                        cum_mdl_in  = cumsum(s_mdl_in)/nyrs_mdl_in;
                        
                        tot_yearly = cum_obs_fut(end);

                        diff_out = cum_mdl_out(ixout)  - cum_obs_fut;
                        diff_in  = cum_mdl_in(ixin) - cum_obs_fut;
                        
                        yl = [min(min(diff_out), -20), max(max(diff_out), 20)];
                        out = diff_in > yl(2) | diff_in < yl(1);
                        in = ~out;
                        diff_in = max(yl(1),min(yl(2), diff_in));
                        
%                         diff_out_pct = diff_out ./ cum_obs_fut * 100.0;
%                         diff_in_pct  = diff_in  ./ cum_obs_fut * 100.0;
% 
                        subplot(5,1,5);
                        plot([0,100],[0,0], "r-", "linewidth",3);
                        hold on;
                        h1 = scatter(x(in) ,diff_in(in), 15, [0,.75,0],"filled");
                        h2 = scatter(x(out),diff_in(out), 10, 'm',"filled");
                        h3 = scatter(x,diff_out, 12, "b","filled");
                        hold off;
                        title(sprintf("Diff Cumulative Precip (Mdl - Obs), vs Count-Quantile, mm/yr   (yearly total: %.0f mm)", tot_yearly),"interpreter","none");
                        xlabel("percentile")
                        ylabel("mm/year");
                        legend([h1,h2,h3], "mdl in","mdl in (off-scale)", "downscaled","location","northwest");
%                       xlim([50,100]);
%                         yl = ylim();
%                         ylim([min(-20,yl(1)),max(20,yl(2))]);
                        ylim(yl);
                        if (range(yl) < 50)
                            yticks(-100:5:100);
                        end                        
                        grid on;
                        
                        figname = sprintf("quantile_%s.tif", underscore(siteLbl));
                        saveas(h, figname);
            
                        

        figure(fignum+1);  
        subplot(2,1,1);

        semilogx(bins, h_mdl/nyrs_mdl, "k:", "linewidth",2);
        hold on;
        semilogx(bins, h_in/nyrs_obs, "g-", "linewidth",2);
        semilogx(bins, hobs_fut/nyrs_obs_fut,"b-", "linewidth",2);
        semilogx(bins, hmdl_out/nyrs_mdl,"r--", "linewidth",2);
        hold off;
        grid on;
        title(sprintf("histogram, %s", siteLbl));
        legend("mdl_in", "obs_in", "obs_fut", "mdl_out","interpreter","none");
        xlabel("precip (log scale)");
        ylabel("wet day counts/year");

        cdf_obs_fut = cumsum(hobs_fut/sum(hobs_fut));
        cdf_mdl_out = cumsum(hmdl_out/sum(hmdl_out));
        cdf_mdl_in  = cumsum(h_mdl/sum(h_mdl));

        subplot(2,1,2);
        semilogx(bins, cdf_mdl_in, "k:", "linewidth",2);
        hold on;
        plot(bins, cdf_mdl_out, "g-", "linewidth",2);
        plot(bins, cdf_obs_fut, "b-", "linewidth",2);
        plot(bins, cdf_mdl_out, "r--", "linewidth",2);
        hold off;
        grid on;
        title(sprintf("cumulative precip, %s", siteLbl));
        legend("mdl_in", "obs_in", "obs_fut", "mdl_out","interpreter","none");
        xlabel("precip (log scale)");
        ylabel("probability");
    end
    
    tot_count = sum(obs_in > 0,"omitnan");
    tot_prcp  = sum(obs_in,"omitnan")/nyrs_obs;
    

    if(do_compare)
        probs = [.01,.05,.1,.2:.9,.95,.99]; %#ok<NASGU>
        vals = [prcp_min-1e-15, 1:20,30:10:maxval+10];
        fprintf("%s\n", siteLbl);
        fprintf("-------------------------------------------Precip per year, mm----------------------------------\n");
        fprintf("  prcp cnt_pct  obscnt mdlcnt   diff     pct       sum_pct obs_sum  mdl_sum     diff      pct   \n");  
        for i=2:length(vals)
            cnt_in   = sum(obs_in>0 & obs_in <= vals(i));
            cnt_frac =  cnt_in/tot_count*100;
            sum_in   = sum(obs_in(obs_in >0 & obs_in <= vals(i)),"omitnan")/nyrs_obs;
            sum_frac = sum_in / tot_prcp*100;
            cnt_out  = sum(mdl_out>0 & mdl_out <= vals(i));
            sum_out  = sum(mdl_out(mdl_out > 0 & mdl_out <= vals(i)),"omitnan")/nyrs_mdl;
            cnt_diff = cnt_in-cnt_out;        
            sum_diff = sum_in - sum_out;
            pct_cnt  = cnt_diff ./ cnt_in * 100;
            pct_diff = sum_diff ./ sum_in * 100;

            DSP_base.print_log("%6.1f  %6.1f  %6d %6d %6d  %6.2f %%    %8.1f %8.1f %8.1f %8.1f %8.2f %%\n", vals(i), ...
                                  cnt_frac,    cnt_in,    cnt_out,    cnt_diff,    pct_cnt,    sum_frac,    sum_in,    sum_out,    sum_diff,    pct_diff);
        end

        mdl_in  = sort(mdl_in(mdl_in>0));
        obs_fut = sort(obs_fut(obs_fut>0));
        mdl_out = sort(mdl_out(mdl_out>0));

        sum_mdl_in = cumsum(mdl_in);
        sum_obs_fut = cumsum(obs_fut);
        sum_mdl_out = cumsum(mdl_out);

        tot_mdl_in = sum_mdl_in(end);
        tot_obs_fut = sum_obs_fut(end);
        tot_mdl_out = sum_mdl_out(end);

        mdl_in_cnt = zeros(100,1);
        mdl_in_amt = zeros(100,1);
        obs_fut_cnt = zeros(100,1);
        obs_fut_amt = zeros(100,1);
        mdl_out_cnt = zeros(100,1);
        mdl_out_amt = zeros(100,1);

        mdl_diff_cnt = zeros(100,1);
        mdl_diff_cnt_pct = zeros(100,1);
        mdl_diff_amt = zeros(1-00,1);
        mdl_diff_amt_pct = zeros(100,1);

        num_mdl_in = length(mdl_in);
        num_obs_fut= length(obs_fut);
        num_mdl_out= length(mdl_out);

        DSP_base.print_log("%s\n", siteLbl);
        DSP_base.print_log("pct m_in_cnt mout_cnt diff_cnt   pct m_in_amt mout_amt diff_amt   pct  obs_cnt mout_cnt diff_cnt diff_pct   pct  obs_amt mout_amt diff_amt  diff_pct\n");  
        if (num_mdl_in ~= num_mdl_out)
            DSP_base.warn_log("NOTE: mdl in and mdl out counts differ:  %6d %6d\n", num_mdl_in, num_mdl_out); 
        end
        for i=1:100

            x = i/100*num_obs_fut;
            obs_fut_cnt(i) = my_interp(1:num_obs_fut, obs_fut, x);
            x = i/100*tot_obs_fut;
            obs_fut_amt(i) = my_interp(sum_obs_fut, obs_fut, x);

            x = i/100*num_mdl_in;
            mdl_in_cnt(i) = my_interp(1:num_mdl_in, mdl_in, x);
            x = i/100*tot_mdl_in;
            mdl_in_amt(i) = my_interp(sum_mdl_in, mdl_in, x);

            x = i/100*num_mdl_out;
            mdl_out_cnt(i) = my_interp(1:num_mdl_out, mdl_out, x);
            x = i/100*tot_mdl_out;
            mdl_out_amt(i) = my_interp(sum_mdl_out, mdl_out, x);

            mdl_diff_cnt(i) = obs_fut_cnt(i) - mdl_out_cnt(i);
            mdl_diff_cnt_pct(i) = (obs_fut_cnt(i) - mdl_out_cnt(i))/obs_fut_cnt(i)*100;

            mdl_diff_amt(i) = obs_fut_amt(i) - mdl_out_amt(i);
            mdl_diff_amt_pct(i) = (obs_fut_amt(i) - mdl_out_amt(i))/obs_fut_amt(i)*100;

            DSP_base.print_log("%3d %8.2f %8.2f %8.2f   %3d %8.2f %8.2f %8.2f   %3d %8.2f %8.2f %8.2f %6.1f %%  %3d %8.2f %8.2f %8.2f %6.2f %%  \n", i, ...
                    mdl_in_cnt(i), mdl_out_cnt(i), mdl_in_cnt(i) - mdl_out_cnt(i), ... 
                    i, ...
                    mdl_in_amt(i), mdl_out_amt(i), mdl_in_amt(i) - mdl_out_amt(i), ... 
                    i, ...
                    obs_fut_cnt(i), mdl_out_cnt(i), mdl_diff_cnt(i), mdl_diff_cnt_pct(i), ... 
                    i, ...
                    obs_fut_amt(i), mdl_out_amt(i), mdl_diff_amt(i), mdl_diff_amt_pct(i) ... 
                    );
        end
    end
    
    fid=fopen("gfdl_test.txt","a");

    worst_ix_mm  = 49+find(abs(mdl_diff_cnt_pct(50:95)) == max(abs(mdl_diff_cnt_pct(50:95))),1,"last");
    worst_ix_pct = 49+find(abs(mdl_diff_cnt_pct(50:95)) == max(abs(mdl_diff_cnt_pct(50:95))),1,"last");
    mae_50 = mean(abs(mdl_diff_cnt(1:49)));
    mae_95 = mean(abs(mdl_diff_cnt(50:95)));
    mae_99 = mean(abs(mdl_diff_cnt(96:99)));
    out_ix = [50, 75, 90, 95, 98,99];
    if (~any(out_ix == worst_ix_mm))
        out_ix = sort([out_ix, worst_ix_mm]);
    end
    if (~any(out_ix == worst_ix_pct))
        out_ix = sort([out_ix, worst_ix_pct]);
    end
    DSP_base.print_log("\n-----%s summary -----\n", siteLbl);
    DSP_base.print_log("\tMean absolute error: <50 %%  :  %6.3f mm\n\tMean absolute error: 50-95 %%:  %6.3f mm\n\tMean absolute error: 96-99 %%:  %6.3f mm\n\n", mae_50, mae_95, mae_99);
    
    fprintf(fid, "\n-----%s summary -----\n", siteLbl);
    fprintf(fid, "\tMean absolute error: <50 %%  :  %6.3f mm\n\tMean absolute error: 50-95 %%:  %6.3f mm\n\tMean absolute error: 96-99 %%:  %6.3f mm\n\n", mae_50, mae_95, mae_99);
    
    DSP_base.print_log("  pct   obs_mm  mdl_out  diff_mm  diff_pct\n");
    fprintf(fid, "  pct   obs_mm  mdl_out  diff_mm  diff_pct\n");
    for i=1:length(out_ix)
        ipct = out_ix(i);
        DSP_base.print_log("  %3d %8.2f %8.2f %8.2f %6.2f %%  ", ipct,  obs_fut_cnt(ipct), mdl_out_cnt(ipct), mdl_diff_cnt(ipct), mdl_diff_cnt_pct(ipct));
        fprintf(fid, "  %3d %8.2f %8.2f %8.2f %6.2f %%  ", ipct,  obs_fut_cnt(ipct), mdl_out_cnt(ipct), mdl_diff_cnt(ipct), mdl_diff_cnt_pct(ipct));
        if (ipct == worst_ix_mm)
            DSP_base.print_log("  <--- worst_case, mm, bottom 95 %%");
            fprintf(fid, "  <--- worst_case, mm, bottom 95 %%");
        end
        if (ipct == worst_ix_pct)
            DSP_base.print_log("  <--- worst_case, pct, 50th - 95 %%");
            fprintf(fid, "  <--- worst_case, pct, 50th - 95 %%");
        end
        DSP_base.print_log("\n");
        fprintf(fid, "\n");
    end
    fclose(fid);

%     
%     if (show_figs)
%         figure(fignum+1)
% 
%         x=1:100;
%         subplot(2,1,1)
%         plot(x, mdl_in_cnt,"k-", x, obs_fut_cnt,"b-", x, mdl_out_cnt, "r-", x, mdl_in_cnt(i) - mdl_out_cnt(i), "k:", x, obs_fut_cnt(i) - mdl_out_cnt(i), "r:", "linewidth", 2);
%         legend("mdl_in","obs_fut","mdl_out", "mdl diff", "obs diff");
%         title(sprintf("quantiles based on counts, %s", siteLbl));
%         xlabel("percent of total");
%         ylabel("precip, mm");
%         subplot(2,1,2)
%         plot(x, mdl_in_amt,"k-", x, obs_fut_amt,"b-", x, mdl_out_amt, "r-", x, mdl_in_amt(i) - mdl_out_amt(i), "k:", x, obs_fut_amt(i) - mdl_out_amt(i), "r:", "linewidth", 2);
%         legend("mdl_in","obs_fut","mdl_out", "mdl diff", "obs diff");
%         title(sprintf("quantiles based on amounts, %s", siteLbl));
%         xlabel("percent of total");
%         ylabel("precip, mm");
% 
%         figure(fignum+2);
% 
%         subplot(3,1,1);
%         x = 1:99;
%         plot(x, mdl_diff_cnt(1:99));
%         grid on;
%         xlabel("percentile");
%         ylabel("Diff, mm");
%         title(sprintf("Difference by Quantile, %s", siteLbl));
%         ylims = ylim();
%         ylims = [min(ylims(1),-2.5), max(ylims(2),2.5)];
%         ylim(ylims);
% 
%         subplot(3,1,2);
%         plot(x, mdl_diff_cnt_pct(1:99));
%         grid on;
%         xlabel("percentile");
%         ylabel("Diff, %");
%         title(sprintf("Percent Difference, %s",siteLbl));
%         ylims = ylim();
%         ylims = [min(ylims(1),-5), max(ylims(2),5)];
%         ylim(ylims);
% 
% 
%         subplot(3,1,3);
%         nobs = round(.99*num_obs_fut);
%         nmdl = round(.99*num_mdl_out);
%         ix=round(linspace(1,nobs, nmdl));
%         jx = 1:nmdl;
% 
%         ix50 = round(.50*num_obs_fut);
%     %   jx50 = find(ix>=ix50,1);
%         ix95 = round(.95*num_obs_fut);
%     %   jx95 = find(ix>=ix95,1);
%         x50 = obs_fut(ix50);
%         x95 = obs_fut(ix95);
% 
%         difs = obs_fut(ix) - mdl_out(jx);
%         scatter(obs_fut(ix), difs, 5, "filled");
%         grid on;
%         xlabel("obs precip, sorted, mm");
%         ylabel("Diff, mm");
%         title(sprintf("Difference vs Precip Amt, %s", siteLbl));
%         ylims = ylim();
%         ylims = [min(ylims(1),-2.5), max(ylims(2),2.5)];
%         ylim(ylims);
%         hold on;
%         plot([x50,x50], ylims,"r:",[x95,x95], ylims,"r--", "linewidth", 2);
%         hold off;
%         legend("diffs","50th pct","95th pct")
% 
%     end
end   

function q = my_interp(vals1, vals2, x)

    ix   = find(vals1<=x,1,"last");
    if (isempty(ix))
        v1  = 0;
        v2  = 0;
    else
        v1  = vals1(ix);
        v2  = vals2(ix);
    end
    if (ix == length(vals1))
        q = vals2(end);
    else
        dx1 = vals1(ix+1) - vals1(ix);
        dx2 = vals2(ix+1) - vals2(ix);
        frac = (x-v1)/dx1;
        q    = v2 + frac*dx2;
    end
end

function [edges, DA_obs, DA_mdl, DA_hist]  = calc_prcp_binning(DA_obs, DA_mdl, DA_hist, use_hist)
%
%   Create a non-linear set of bin edges so that the combined obs, model (and optionally hist) precip data
%   is approximately gaussian.
%
%   By using the nonlinear binning, we can keep the precip data in its
%   original units, rather than transforming it to another data space, and
%   having to transform it back later.

    prcp_distrib = DA_obs.RP.prcp_distrib; 
    
    p = strfind(DA_obs.DA_title,":");
    if (isempty(p)), p=0; end
    lbl = extractAfter(DA_obs.DA_title,p(1));
    
    nedges = 1000;
    if (use_hist)
        all_raw_data = [to_column(DA_obs.anoms); to_column(DA_mdl.anoms); to_column(DA_hist.anoms)];
    else
        all_raw_data = [to_column(DA_obs.anoms); to_column(DA_mdl.anoms)];
    end
    dmax1 = 1.1*max(all_raw_data,[],"omitnan");
    dmax = 100*ceil(dmax1/100);    
    pr_min = min(.1, min(all_raw_data));
    data = all_raw_data(all_raw_data >= pr_min);      % we assume that all_raw_data has been trimmed and thresholded.

    scaling = 1.0;
    done = false;
    pr_probthresh = .995;   % about 2.58 sigma
    while (~done)
        try
            p = fitdist(to_column(data.^scaling), prcp_distrib);
        catch
            error("fitdist failure %s  %d precip values, %s\n", prcp_distrib, length(data), lbl);
    %       error("fitdist failure processing  stn %d, %s %.4f %.4f, dist %d %s.  %d precip values\n", istn, loc, lat, lon, j, dist, length(data));
        end
        lin_edges  = linspace(pr_min,dmax,nedges);

        y1 = pdf(p, lin_edges);

%       y_pdf = y1; % [0,y1(1:end-1) + dy/2;

            % make sure we have a 0 at the beginning of our probabilities.
        y_pdf = [normpdf(-7.5),y1];
        b0 = .90*pr_min;
        lin_edges = [0, b0, linspace(pr_min,dmax,nedges)];
        diffs = diff(lin_edges);
        nedges = length(lin_edges);
        y_pdf = y_pdf./sum(y_pdf.*diffs);

        Y_cdf=[0,cumsum(y_pdf.*diffs)];                            % CDF for distribution
        Y_cdf = Y_cdf/Y_cdf(end);   % in case it doesn't sum exactly to 1...)

        pr_probmax = interp1(lin_edges,Y_cdf,max(data));
        if (pr_probmax < norminv(.999))
            done = true;
        else
            scaling = scaling*(1.002*pr_probmax/norminv(pr_probthresh));        % .999)

        end
    end

    d = diff(Y_cdf);
    zer = d==0;
    if (any(zer))
  %       fprintf("%-25s : %d zeros in Y_cdf\n", dist, sum(d==0));
        Y_cdf(zer)=[];
        lin_edges(zer)=[];

    end
    
    zedges=linspace(-7.5,9.5,nedges);            % get values for a normal CDF (mu=0, sigma=1)   from -7.5 sigmas to +9.5 sigmas
    pnorm=normcdf(zedges);                       % normal CDF, mean 0, sigma 1;

    edges = interp1(Y_cdf,lin_edges,pnorm);


    if (any(isnan(edges)))
        edges = edges(~isnan(edges));
        zedges = zedges(~isnan(edges));
    end
    
    ok = [diff(edges)>0,true];
    if (any(~ok))
        edges = edges(ok);
        zedges = zedges(ok);
    end
        
    
%   [mu, sig, skew, kurt] = pdf_stats(hcx,zbins,[],true);

        % Set the bin edges for all 3 datasets to the same edges.

    DA_obs.RP.edges  = edges;
    DA_hist.RP.edges = edges;
    DA_mdl.RP.edges  = edges;
    
    DA_obs.RP.zedges  = zedges;
    DA_hist.RP.zedges = zedges;
    DA_mdl.RP.zedges  = zedges;
    
    DA_obs.binning_checked = true;
    DA_hist.binning_checked = true;
    DA_mdl.binning_checked = true;
end   
