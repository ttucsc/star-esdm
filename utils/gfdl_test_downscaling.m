function retval = gfdl_test_downscaling(varname, e_or_c, st_yr, varargin)
%        gfdl_test_downscaling(varnames, e_or_cs, st_yr, + kwd/value pairs for dry_run, latrange, lonrange, stnnums,do_internal_plots, do_parfor, do_extended, do_continue, debug_rg)
%   setups up and runs ARRM_V2 for GFDL Perfect Model Data runs
%       By default will run tasmin & tasmax for cm3 and ecm, for the entire globe
%       
%   Inputs:
%     Required:
%       varname    tasmin, tasmax, pr or both.                         ["tasmin","tasmax"]f
%       e_or_cs     cm3, ecm or both.                                   ["cm3","ecm"]
%       st_yr       start year of data to download:  1979 (historical), or 2086 (tasmax, tasmin) or 2066 (pr)
%     Optional keyword/value pairs
%       do_parfor   true or false.  controls whether to run parallel.   [true]
%       do_continue                                                     [false]
%                   false:  run all gridboxes, regardless of previous completion status.
%                   true:   look in output folder, checks all existing grid_*.nc's for completion_flag, and runs
%                           all missing gridboxes and reruns any gridboxes where the completion_status is not 1 or 2
%                               completion_status of 1:     complete
%                               completion_status of 2:     complete, but all gridcells had insufficient data
%                                           only happens for livneh and similar inputs where water gridcells are all NAs
%       latrange    lat & lon range to run, as 2-element vector.        [ -90,  90]
%       lonrange                                                        [   0, 360]
%                       (lon can be in range -180 - 180 or 0 - 360 
%                   shortcut for CONUS:  set latrange and lonrange to "conus" to do continental US only.
%                       "conus":  lat 20 - 50, lon 230 - 305, which gets parts of mexico & canada as well.
%                       "central":
%                       "test": 
%       debug_rg    boolean.  If true, adds additional debug info to log files.     [false]
%
%       debug_info  info to turn on debugging to save ARRM_V2 output to mat files.
%                       illat(1), illat(2)
%                       illon(1), illon(2)
%                       grid lat start, grid lon start.
%
%                       defaults are set for grid_N25W100, sites 1-5 of top row.
%
%
%   For precip, use:
% gfdl_test_downscaling("pr","esm",2066,"dry_run", false,"do_parfor", true,  "region","global","pdf_map_method","clim","clim_nterms",5,"clim_sig_terms",2,"anom_nterms",2,"anom_sig_terms",[1,2], "do_compare", false, "show_figs",false, "do_clean_files", true);
%
%   For temps, use:
% gfdl_test_downscaling("tasmax","esm",2086,"dry_run", true,"do_parfor", true,  "region","central", "pdf_map_method","linear","do_compare", false, "show_figs",false, "do_clean_files", true);
%
%
%       11/18/2021 icsf added "do_extended" option, and removed "extended, false"
%                       added use_lacie option.  default: false
%                       added pdf_map_method option.  default:  "clim" for precip, "linear" for temperature variables. 
%

    if (nargin < 6)
        retval = 2;
        return;
    end
    
    retval = 1;
    try
        [varname, e_or_c, st_yr, ...
         latrange, lonrange, testsites, lats, lons, region, ...
         pdf_map_method, ...
         sigma_normalize, median_normalize, ...
         prcp_distrib, ...
         gridsize, baseout, ...
         do_clean_files, ...
         do_compare, show_figs, ...
         use_lacie, ...
         do_continue, dry_run, do_internal_plots, exit_on_error, par_info, do_extended, far_outlier_anchor_pt, far_outlier_thresh, Unmatched]  = init_params(varname, e_or_c, st_yr, varargin{:});
%         do_norm_sigmas, ...
%         far_outlier_thresh, ...
     

 
    
%------------ main program begins here------------
    
     
        interp_method = "bilinear";


        if (isfolder("/Users/iscottfl"))
            if (~use_lacie)
    %         basesrc = "/Volumes/lacie_1/data/gcm_cmip5/hires/ian_work";
    %         if (isempty(baseout))
    %             baseout = fullfile("/Volumes/lacie_1/data/downscaled/arrm_v2/GFDL_perfect_model",region);
    %         end
                basesrc = "/Users/iscottfl/gfdl_hires_ian";
                if (isempty(baseout))
                    baseout = fullfile("/Users/iscottfl/downscaled/GFDL_perfect_model",region);
                end
            else
                basesrc = "/Volumes/lacie_1/data/gfdl_hires";
                if (isempty(baseout))
                    baseout = fullfile("/Volumes/lacie_1/data/downscaled/arrm_v2/GFDL_perfect_model",region);
                end
            end            
        elseif (isfolder("/lustre/research"))
            basesrc = "/lustre/research/hayhoe/gfdl_hires_ian/joined";
            if (isempty(baseout))
                baseout = fullfile("/lustre/scratch/iscottfl/downscaled/GFDL_perfect_model",region);
            end
        end


        if (strcmp(e_or_c,"esm"))
            Xstuff="X1X2X3";
        else 
            Xstuff="X1X3X4";
        end

        ensemble="r1i1p1";
        model=sprintf("GFDL-HIRAM-C360_%s", e_or_c);
        scenario=sprintf("rcp85-%s", Xstuff);

        hist_yrs = [1979,2038];
        obs_yrs = [1979,2038];
        obsname  = sprintf("%s_day_GFDL-HIRAM-C360_amip_19790101-20381231.H2H3.hires.llt.nc",varname);
        histname = sprintf("%s_day_GFDL-HIRAM-C360_amip_19790101-20381231.H2H3.lowres.llt.nc",varname);
        if (st_yr == 2066)
            if (~strcmp(varname,"pr")), error("error:  use 2086 for future runs of %s", varname); end
            mdl_yrs = [2066,2095];
            pdf_yrstep = 30;
            yr_lbl = "future";
    %       mdldir   = basesrc; %
            mdldir   = fullfile(basesrc, "future");
            mdlname  = sprintf("c360_hiram_%s_rcp85_%s.lowres.%s.20660101-20951231.llt.nc",e_or_c, Xstuff, varname);
            if (do_compare)
                obsname_fut =  fullfile(mdldir, sprintf("c360_hiram_%s_rcp85_%s.hires.%s.20660101-20951231.llt.nc",e_or_c, Xstuff, varname));
            else
                obsname_fut = strings(0);
            end
        elseif (st_yr == 2086)
            if (strcmp(varname,"pr")), error("error:  use 2066 for future runs of %s", varname); end
            mdl_yrs = [2086,2115];
            pdf_yrstep = 30;
            yr_lbl = "future";
    %       mdldir   = basesrc; %
            mdldir   = fullfile(basesrc, "future");
            mdlname  = sprintf("c360_hiram_%s_rcp85_%s.lowres.%s.20860101-21151231.llt.nc",e_or_c, Xstuff, varname);
            if (do_compare)
                obsname_fut =  fullfile(mdldir, sprintf("c360_hiram_%s_rcp85_%s.hires.%s.20860101-21151231.llt.nc",e_or_c, Xstuff, varname));
            else
                obsname_fut = strings(0);
            end
        elseif (st_yr == 1979)
            mdl_yrs = [1979,2038];
            pdf_yrstep = 60;
            yr_lbl = "hist";
    %       mdldir   = basesrc; 
            mdldir  = fullfile(basesrc, "hist");
            mdlname  = sprintf("%s_day_GFDL-HIRAM-C360_amip_19790101-20381231.H2H3.lowres.llt.nc", varname);
            if (do_compare)
                obsname_fut  = fullfile(mdldir, sprintf("%s_day_GFDL-HIRAM-C360_amip_19790101-20381231.H2H3.hires.llt.nc", varname));
            else
                obsname_fut = strings(0);
            end
        else
            error("bad start yr: %d.  must be 1979,2066 or 2086");
        end
    %           trend_yrs = mdl_yrs;
        if (is_prcp_variable(varname))      % this is forced to false for precip runs, but we set it here for clarity.
            trend_order = 0;
        else
            trend_order = 1;
        end

        if (strcmp(region,"test_stations"))
            outname=sprintf("downscaled.%s.%s.%s.%s.gfdl_PM.%s.%4d.%4d.%%s.mat", model, ensemble, varname, scenario, region, mdl_yrs(1), mdl_yrs(2));
        else
            outname=sprintf("downscaled.%s.%s.%s.%s.gfdl_PM.%s.%4d.%4d.LLGRID.nc", model, ensemble, varname, scenario, region, mdl_yrs(1), mdl_yrs(2));
        end

        obsdir   = fullfile(basesrc, "hist");
        histdir  = fullfile(basesrc, "hist");
        if (strcmp(region,"gfdl_test_stations"))
            outsubdir = "downscaled_test";
        elseif (par_info.do_parfor)
            outsubdir = join([yr_lbl, model,ensemble,varname,scenario,"parallel"],"_");
        else
            outsubdir = join([yr_lbl, model,ensemble,varname,scenario],"_");
        end
        outdir   = fullfile(baseout, outsubdir);

                % put together a run description.
        descr = sprintf("GFDL Perfect Model Run  (%s,%s,%d) %s\nParameters:\n", varname,e_or_c,st_yr, datestr(now, "yyyy-mm-dd HH:MM:SS"));
        descr = strcat(descr, sprintf("model:         %s\n", model));
        descr = strcat(descr, sprintf("ensemble:      %s\n", ensemble));
        descr = strcat(descr, sprintf("varname:       %s\n", varname));
        descr = strcat(descr, sprintf("scenario:      %s\n", scenario));
        descr = strcat(descr, sprintf("region:        %s\n", region));
        descr = strcat(descr, sprintf("obsname:       %s\n",  vec2string(obsname, "format","%s ","brackets",'[]')));
        descr = strcat(descr, sprintf("histname:      %s\n", vec2string(histname,"format","%s ","brackets",'[]')));
        descr = strcat(descr, sprintf("mdlname:       %s\n", vec2string(mdlname, "format","%s ","brackets",'[]')));
        descr = strcat(descr, sprintf("outname:       %s\n", outname));
        descr = strcat(descr, sprintf("obsname_fut:   %s\n", obsname_fut));
        descr = strcat(descr, sprintf("obsdir:        %s\n", obsdir));
        descr = strcat(descr, sprintf("histdir:       %s\n", histdir));                
        descr = strcat(descr, sprintf("mdldir:        %s\n", mdldir));
        descr = strcat(descr, sprintf("outdir:        %s\n", outdir));
        descr = strcat(descr, sprintf("obs_yrs:       %s\n", vec2string(obs_yrs,  "format","%4d ","brackets",'[]')));
        descr = strcat(descr, sprintf("hist_yrs:      %s\n", vec2string(hist_yrs, "format","%4d ","brackets",'[]')));
        descr = strcat(descr, sprintf("mdl_yrs:       %s\n", vec2string(mdl_yrs,  "format","%4d ","brackets",'[]')));
        if (~isempty(lats))
            descr = strcat(descr, sprintf("lats:          %s\n", vec2string(lats, "format","%9.4f ", "brackets",'[]')));
            descr = strcat(descr, sprintf("lons:          %s\n", vec2string(lons, "format","%9.4f ", "brackets",'[]')));
            descr = strcat(descr, sprintf("sites:         %s\n", vec2string([testsites(1:min(5,length(testsites)));"..."])));
        else
            descr = strcat(descr, sprintf("latrange:      %s\n", vec2string(latrange, "format","%9.4f ", "brackets",'[]')));
            descr = strcat(descr, sprintf("lonrange:      %s\n", vec2string(lonrange, "format","%9.4f ", "brackets",'[]')));
            descr = strcat(descr, sprintf("gridsize:      %s\n", vec2string(gridsize, "format","%2d ",   "brackets",'[]')));
            descr = strcat(descr, sprintf("do_continue:   %d\n", do_continue));
        end
%       descr = strcat(descr, sprintf("do_norm_sigmas:%d\n", do_norm_sigmas));
        descr = strcat(descr, sprintf("sigma_normalize:%d\n", sigma_normalize));
        descr = strcat(descr, sprintf("median_normalize:%d\n", median_normalize));
        if (is_prcp_variable(varname))
            descr = strcat(descr, sprintf("prcp_distrib:  %s\n", prcp_distrib));
        end
        descr = strcat(descr, sprintf("interp_method: %s\n", interp_method));
    %   descr = strcat(descr, sprintf("trend_yrs:     %s\n", vec2string(trend_yrs,"format","%4d ","brackets",'[]')));
        descr = strcat(descr, sprintf("pdf_yrstep:    %d\n", pdf_yrstep));
        descr = strcat(descr, sprintf("pdf_yrs:       %d\n", pdf_yrstep));
        descr = strcat(descr, sprintf("trend_order:   %d\n", trend_order));
%       descr = strcat(descr, sprintf("far_outlier_thresh:     %.2f\n", far_outlier_thresh));
        descr = strcat(descr, sprintf("do_internal_plots: %s\n", vec2string(do_internal_plots, "format","%2d ",   "brackets",'[]')));
        descr = strcat(descr, sprintf("far_outlier_anchor_pt:     %d\n", far_outlier_anchor_pt));
        descr = strcat(descr, sprintf("far_outlier_thresh:     %d\n", far_outlier_thresh));
        descr = strcat(descr, sprintf("do_parfor:     %d\n", par_info.do_parfor));
        if (par_info.do_parfor)
            descr = strcat(descr, sprintf("cluster:       %s\n", par_info.cluster));
            descr = strcat(descr, sprintf("nnodes:        %d\n", par_info.nnodes));
            descr = strcat(descr, sprintf("nworkers:      %d\n", par_info.nworkers));
        end
        descr = strcat(descr, sprintf("do_extended:   %d\n", do_extended));
        descr = strcat(descr, sprintf("exit_on_error: %d\n", exit_on_error));
    %   descr = strcat(descr, sprintf("debug_rg:      %d\n", debug_rg));
        descr = strcat(descr, sprintf("do_clean_files %d\n", do_clean_files));
        descr = strcat(descr, sprintf("do_compare     %d\n", do_compare));
    %   descr = strcat(descr, sprintf("show_figs      %d\n", show_figs));
        descr = strcat(descr, sprintf("additional params:\n"));
        ufields=fields(Unmatched);
        for i=1:length(ufields)
            descr = strcat(descr, sprintf("     %-15s      %s\n", ufields{i},vec2string(Unmatched.(ufields{i}))));
        end
        run_description = strcat(descr, newline);

        if (dry_run)
            fprintf("Dry Run: %s\n", run_description);

            if (~isfile(fullfile(obsdir,  obsname)))  
                fprintf(2, "missing  obs file: %s\n", fullfile(obsdir, obsname)); 
            else
                fprintf(   "obs file exists:   %s\n", fullfile(obsdir, obsname)); 
            end
            if (~isfile(fullfile(histdir, histname)))
                fprintf(2, "missing hist file: %s\n", fullfile(histdir, histname)); 
            else
                fprintf(   "hist file exists:  %s\n", fullfile(histdir, histname)); 
            end
            if (~isfile(fullfile(mdldir, mdlname)))
                fprintf(2, "missing  mdl file: %s\n", fullfile(mdldir, mdlname)); 
            else
                fprintf(   "mdl file exists:   %s\n", fullfile(mdldir, mdlname));             
            end
            if (~isfolder(baseout))
                fprintf(2, "missing output folder: %s\n", baseout); 
            else
                fprintf(   "output folder exists:  %s\n", baseout); 
            end
            fprintf("\n\n");
        else   
            fprintf("%s\n", run_description);
            setup_dirs(baseout, outsubdir)
            if (do_clean_files)
                clean_dirs(baseout, outsubdir);
            end
            if (~isempty(lats))         % if user specified either test_stations or test_sites
                 ARRM_V2_wrapper(model,ensemble,varname,scenario,obs_yrs,hist_yrs,mdl_yrs,...
                                 "run_description", run_description, ...
                                 "obsnames",    obsname, ...
                                 "histnames",   histname, ...
                                 "mdlnames",    mdlname, ...
                                 "outname",     outname, ...
                                 "region",      region, ...
                                 "obsdir",      obsdir,...
                                 "histdir",     histdir, ...
                                 "mdldir",      mdldir,...
                                 "outdir",      outdir, ...
                                 "testsites",   testsites, ...
                                 "lats",    lats, "lons",lons, ...
                                 "llgrid_size", gridsize, ...
                                 "do_continue", false, ...
                                 "interp_method",interp_method, ...
                                 "pdf_map_method",pdf_map_method, ...
                                 "sigma_normalize",sigma_normalize, ...
                                 "median_normalize",median_normalize, ...
                                 "prcp_distrib",prcp_distrib, ...
                                 "pdf_yrstep",  pdf_yrstep, ...     % GFDL-PM-specific!
                                 "pdf_yrs",     pdf_yrstep, ...     % GFDM-PM-specific.  either doing all hist as 1 period, or all future as 1 period.
                                 "trend_order", trend_order, ...    % GFDL-PM-specific!
                                 "do_internal_plots",do_internal_plots, ...
                                 "figbase",     5000, ...
                                 "exit_on_error",exit_on_error, ...     
                                 "obsname_fut", obsname_fut, ...
                                 "do_compare",  do_compare, ...
                                 "show_figs",   show_figs, ...
                                 "do_parfor",par_info.do_parfor, "cluster",par_info.cluster, "nnodes", par_info.nnodes, "nworkers", par_info.nworkers, ...
                                 "do_extended", do_extended, ...
                                 "far_outlier_anchor_pt", far_outlier_anchor_pt, ...
                                 "far_outlier_thresh", far_outlier_thresh, ...
                                 Unmatched ...
                                 );
%                                  "do_normalize_sigmas",do_norm_sigmas, ...
                 
            else
    %                                  "debug_rg",    debug_rg, ...
    %                                      "trend_yrs",   trend_yrs, ...      % don't need to pass trend_yrs in.  ARRM_V2_run_gridded
    %                                                                     restricts the trend_yrs to the available data.
                 [~,retval] = ARRM_V2_wrapper(model,ensemble,varname,scenario,obs_yrs,hist_yrs,mdl_yrs,...
                                 "run_description", run_description, ...
                                 "obsnames",    obsname, ...
                                 "histnames",   histname, ...
                                 "mdlnames",    mdlname, ...
                                 "outname",     outname, ...
                                 "region",      region, ...
                                 "obsdir",      obsdir,...
                                 "histdir",     histdir, ...
                                 "mdldir",      mdldir,...
                                 "outdir",      outdir, ...
                                 "latrange",    latrange, "lonrange",lonrange, ...
                                 "llgrid_size", gridsize, ...
                                 "do_continue", do_continue, ...
                                 "interp_method",interp_method, ...
                                 "pdf_map_method",pdf_map_method, ...
                                 "sigma_normalize",sigma_normalize, ...
                                 "median_normalize",median_normalize, ...
                                 "prcp_distrib",prcp_distrib, ...
                                 "pdf_yrstep",  pdf_yrstep, ...     % GFDL-PM-specific!
                                 "pdf_yrs",     pdf_yrstep, ...     % GFDM-PM-specific.  either doing all hist as 1 period, or all future as 1 period.
                                 "trend_order", trend_order, ...    % GFDL-PM-specific!
                                 "do_internal_plots",do_internal_plots, ...
                                 "figbase",     5000, ...
                                 "exit_on_error",exit_on_error, ... 
                                 "obsname_fut", obsname_fut, ...
                                 "do_compare",  do_compare, ...
                                 "show_figs",   show_figs, ...
                                 "do_parfor",par_info.do_parfor, "cluster",par_info.cluster, "nnodes", par_info.nnodes, "mb_per_core", par_info.mb_per_core, "nworkers", par_info.nworkers, ...
                                 "extended", do_extended, ...
                                 "far_outlier_anchor_pt", far_outlier_anchor_pt, ...
                                 "far_outlier_thresh", far_outlier_thresh, ...
                                 Unmatched ...
                                 );
%                                  "do_normalize_sigmas",do_norm_sigmas, ...

            end
        end
    catch me
        report_me_error(me, mfilename(), 1, false);
        retval = 3;
    end
end

%     if (~exist("debug_rg","var")),                                      debug_rg = false;                 end
      

function [varname, e_or_c, st_yr, ...
          latrange, lonrange, testsites, lats, lons, region, ...
          pdf_map_method, ...
          sigma_normalize, median_normalize, ...
          prcp_distrib, ...
          gridsize, baseout, ...
          do_clean, ...
          do_compare, show_figs, ...
          use_lacie, ...
          do_continue, dry_run, do_internal_plots, exit_on_error, par_info, do_extended, far_outlier_anchor_pt, far_outlier_thresh, Unmatched]  = init_params(varname, e_or_c, st_yr, varargin)
%          normalize_sigmas, ...
%          far_outlier_thresh, ...
     
    p = inputParser;
    p.StructExpand = true;
    p.KeepUnmatched = true;
    
    testsites = strings(0);
    
    addRequired(p,"varname",                    @(s) any(strcmp(s,["tasmax","tasmin","pr"])));
    addRequired(p,"e_or_c",                     @(s) any(strcmp(s,["cm3","esm"])));
    addRequired(p,"st_yr",                      @(s) any(s==[1979, 2066, 2086]));

    addParameter(p,'latrange',      [25,30],    @(s) isnumeric(s) && (isempty(s) || length(s)==2));
    addParameter(p,'lonrange',      [-100,-95], @(s) isnumeric(s) && (isempty(s) || length(s)==2));
    addParameter(p,'lats',          [],         @(s) isnumeric(s));
    addParameter(p,'lons',          [],         @(s) isnumeric(s));
    addParameter(p,'region',        strings(0), @(s) isstring(s));
    addParameter(p,'stnnums',       [],         @(s) isnumeric(s) && (isempty(s) || (all(mod(s,1)==0) && all(s>0) && all(s<=24))));
    addParameter(p,'baseout',       strings(0), @(s) isempty(s) || isstring(s));
    addParameter(p,'use_lacie',     false,      @(s) islogical(s));     % if true and on neys or icsf_jmac, write to /Volumes/lacie_1/data/downscaled ;  else write to ~/downscaled, which is SSD and much faster.
    addParameter(p,"pdf_map_method",strings(0), @(s) isstring(s) && (isempty(s) || any(strcmp(s,["","clim","linear"]))));    
    addParameter(p,'gridsize',      [5,5],      @(s) isnumeric(s) && length(s)==2);
%   addParameter(p,"far_outlier_thresh", 0,     @(s) isnumeric(s));
    addParameter(p,"do_parfor",     true,       @(s) islogical(s) || any(s==[0,1]) || isempty(s));
    addParameter(p,"do_extended",   true,       @(s) islogical(s) || any(s==[0,1]) || isempty(s));
    addParameter(p,"cluster",       strings(0), @(s) isstring(s) && any(strcmp(s,["local","nocona","redraider R2020b"])));
    addParameter(p, "nworkers",     [],         @(s) isempty(s) || (isnumeric(s) && mod(s,1)==0 && (s>=1 && s<=512))); 
    addParameter(p, "nnodes",       1,          @(s) isinteger(s) && (s>=1 && s<=4)); 
    addParameter(p, "mb_per_core",  5200,       @(s) isempty(s) || (isnumeric(s) && s>0 && s<=515565));    % applies only to HPCC runs.  default will allow 48 notes in 250GB or 96 workers in 500G
    addParameter(p,"do_clean_files",true,       @(s) islogical(s) || any(s==[0,1]));
    addParameter(p,"do_continue",   false,      @(s) islogical(s) || any(s==[0,1]));
    addParameter(p,"dry_run",       false,      @(s) islogical(s) || isempty(s));
    addParameter(p,"do_compare",    [],         @(s) islogical(s) || isempty(s));
    addParameter(p,"show_figs",     [],         @(s) islogical(s) || isempty(s));
    addParameter(p,"exit_on_error", [],         @(s) islogical(s) || isempty(s));
    addParameter(p,"sigma_normalize", true,     @(s) islogical(s) || isempty(s));
    addParameter(p,"median_normalize", false,   @(s) islogical(s) || isempty(s));
    addParameter(p,"do_internal_plots", [],     @(s) isnumeric(s) && (isempty(s) || (all(mod(s,1)==0) && all(s>=0) && all(s<64))));
    addParameter(p,"prcp_distrib",  strings(0), @(s) isstring(s) || isempty(s));
    addParameter(p,"far_outlier_anchor_pt", .1587, @(s) isnumeric(s));
    addParameter(p,"far_outlier_thresh",    [], @(s) isnumeric(s) || isempty(s));
    
    parse(p,varname, e_or_c, st_yr, varargin{:});
    varname             = string(p.Results.varname);
    e_or_c              = string(p.Results.e_or_c);
    st_yr               = p.Results.st_yr;
    latrange            = p.Results.latrange;
    lonrange            = p.Results.lonrange;
    lats                = p.Results.lats;
    lons                = p.Results.lons;
    region              = p.Results.region;
    stnnums             = p.Results.stnnums;
    
    baseout             = p.Results.baseout;
    use_lacie           = p.Results.use_lacie;
    pdf_map_method      = p.Results.pdf_map_method;
    gridsize            = p.Results.gridsize;
    dry_run             = p.Results.dry_run;
%    far_outlier_thresh  = p.Results.far_outlier_thresh;
    do_parfor           = p.Results.do_parfor;
    do_extended         = p.Results.do_extended;
    nnodes              = p.Results.nnodes;
    mb_per_core         = p.Results.mb_per_core;
    nworkers            = p.Results.nworkers;
    pcluster            = p.Results.cluster;
    do_clean            = p.Results.do_clean_files;
    do_continue         = p.Results.do_continue;
    do_compare          = p.Results.do_compare;
    show_figs           = p.Results.show_figs;
    exit_on_error       = p.Results.exit_on_error;
    do_internal_plots   = p.Results.do_internal_plots;
%   normalize_sigmas    = p.Results.normalize_sigmas;
    sigma_normalize     = p.Results.sigma_normalize;
    median_normalize    = p.Results.median_normalize;
    prcp_distrib        = p.Results.prcp_distrib;
    far_outlier_anchor_pt= p.Results.far_outlier_anchor_pt;
    far_outlier_thresh  = p.Results.far_outlier_thresh;
    
    Unmatched = p.Unmatched;

    [hostname, on] = get_hostname();
    if (on.hpcc_system)
        onhost = "nocona";
    else
        onhost = hostname;
    end

%   if (isempty(sigma_normalize)), error("must specify sigma_normalize (to true or false)"); end
%     if(is_empty(normalize_sigmas))
%         if (is_temp_variable(varname))
%             normalize_sigmas = true;
%         else
%             normalize_sigmas = false;
%         end
%     end
     
    if (isempty(median_normalize))
        median_normalize = is_prcp_variable(varname);
    end

    if (isempty(far_outlier_thresh))
        if (is_prcp_variable(varname))
            far_outlier_thresh = .0228;
        else
            far_outlier_thresh = 1e-4;
        end
    end

    if (do_parfor)
        if (strcmp(onhost,"nocona"))
            if (isempty(pcluster)), pcluster = "redraider R2020b"; end
            if (strcmp(pcluster,"redraider R2020b"))
                maxworkers = 128*nnodes;
            else
                maxworkers = 120;
            end
        elseif (strcmp(onhost,"neys"))
            if (isempty(pcluster)), pcluster = "local"; end
            if (~strcmp(pcluster,"local"))
                error("error:  cluster %s not available on neys", pcluster);
            end
            maxworkers = 18;
        elseif (any(strcmp(onhost,["icsf-kmac","icsf-lmac"])))
            maxworkers = 8;
        else
            myCluster = parcluster('local');
            maxworkers = myCluster.NumWorkers;
        end
                 
        if (isempty(nworkers))
            nworkers = maxworkers;                
        else
            nworkers = min(nworkers, maxworkers);
        end
    end
            
    par_info = struct("do_parfor", do_parfor, "cluster", pcluster, "nworkers", nworkers, "nnodes", nnodes, "mb_per_core", mb_per_core);
    
    
    if (isempty(region))
        region = "misc_run";
        if (isempty(latrange) || isempty(lonrange))
            error("no region or lat/lon specified");
        end
        if (isempty(exit_on_error)), exit_on_error = false; end
    elseif (strcmpi(region,"conus"))
        latrange = [20,50];
        lonrange = [-130,-55];
        region = "conus";
        if (isempty(exit_on_error)), exit_on_error = false; end
    elseif (strcmpi(region,"global"))
        latrange = [-90,90];
        lonrange = [0,360];
        region = "global";
        if (isempty(exit_on_error)), exit_on_error = false; end
    elseif (strcmpi(region,"test_stations"))
        load("test_stations_25.mat","siteTbl");
        if (~isempty(stnnums))
            siteTbl = siteTbl(stnnums,:);
        end
        latrange = [];
        lonrange = [];
        lats = siteTbl.lat;
        lons = siteTbl.lon;
        testsites=string(siteTbl.stnName);
        par_info.do_parfor = false;
        if (isempty(do_compare)), do_compare = true; end
        if (isempty(show_figs)), show_figs = true; end
        if (isempty(exit_on_error)), exit_on_error = true; end
    elseif (strcmpi(region,"central"))
        latrange = [28,48];
        lonrange = [-120,-90];
        region = "central";
        if (isempty(exit_on_error)), exit_on_error = false; end
    elseif (strcmpi(region,"test_region"))
        latrange = [32,34];
        lonrange = [-102,-100];
        region = "test_region";
        if (isempty(exit_on_error)), exit_on_error = false; end
    elseif (strncmpi(region,"site",4))   % for testing an arbitrary site.  Will adjust site(s) to nearest lat/lon in GFDL hi-res data
        if (isempty(lats))
            lats = latrange;
            lons = lonrange;
        end
        nlats = length(lats);
        latrange = [];
        lonrange = [];
        testsites = strings(nlats,1);
        for i=1:nlats
            testsites(i) = sprintf("test_site %.4f %.4f",lats(i),lons(i));
        end
        if (isempty(do_compare)), do_compare = true; end
        if (isempty(show_figs)), show_figs = true; end
        par_info.do_parfor = false;
        if (isempty(exit_on_error)), exit_on_error = true; end
    end
    
    if ((isempty(region) || strcmp(region,"global") || strcmp(region,"misc_run")) && (all(lonrange==[0,360]) || all(lonrange==[-180,180])))
        lonrange(end) = lonrange(end) - .00001;   % so modulus won't screw it up...
        region = "global";
    elseif (isempty(region))
        region = "custom";
        if (isempty(exit_on_error)), exit_on_error = false; end
    end  
    
    if (isempty(do_compare)), do_compare = false; end
    if (isempty(show_figs)), show_figs = false; end
    
    if (isempty(pdf_map_method))
        if (is_prcp_variable(varname))
            pdf_map_method = "clim";
            if (isempty(prcp_distrib) || strlength(prcp_distrib)==0)
                fprintf("using root-scaling for precip\n");
                prcp_distrib = "pwr";
%               prcp_distrib = "loglogistic";       % the old default.  pwr (max .2) does slightly better.
            end
        elseif (is_temp_variable(varname))
            pdf_map_method = "linear";
            if (~isempty(prcp_distrib))
                if (strlength(prcp_distrib)==0)
                    prcp_distrib = strings(0);
                else
                    error("prcp_distrib is specified, but run is not a precip run");
                end
            end
        else
            error("unknown varname: %s\n", varname);
        end
    end

end
    



%-------------------------------------------------------------------

% From setup_test.m, my original script for running gfdl perfect model tests

%         nc=ncdf(obsname);
%         nc.loadvars([],true);
%         alllats=nc.getvardata('lat');
%         alllons=nc.getvardata('lon');
%         % dlat = abs(testlat-alllats);
        % dlon = abs(testlon-alllons);
        % [~,ixlat] = sort(dlat);
        % [~,ixlon] = sort(dlon);
        % lats = sort(alllats(ixlat(1:50)));
        % lons = sort(alllons(ixlon(1:50)));

        % large region around Denver

        % lonkeepers = alllons > 248 & alllons < 262;
        % latkeepers = alllats >  34 & alllats <  46;

        % small region around denver
        % lonkeepers = alllons > 255 & alllons < 260;
        % latkeepers = alllats >  35 & alllats <  40;
        % 

        % continental US:

        % lonkeepers = alllons > 235 & alllons < 300;
        % latkeepers = alllats >  20 & alllats <  50;
        % 
        % lats = alllats(latkeepers);
        % lons = alllons(lonkeepers);

        % lats = alllats;
        % lons = alllons;
        % 
        % lonkeepers = alllons >= 355 | alllons < 5;
        % lons = alllons(lonkeepers);

%         test_ARRM("obsYrs",obsyrs,"histYrs",histyrs,"mdlYrs",mdlyrs, ...
%                   "obsname",obsname,"histname",histname,"mdlname",mdlname, ...
%                   "latrange",latrange, "lonrange",lonrange, "llgrid_size",[5,5], ...
%                   "model",model,"ensemble",ensemble, "scenario", scenario, "varname",varname, ...
%                   "outdir",outdir, ...
%                   'cdf_append_type','sk_normal', ...
%                   'pdf_yrstep', pdf_yrstep, 'trend_order', trend_order, "do_parfor",true);
%               
