function gridded_test_downscaling(varname, model, ensemble, scenario, dry_run, latrange, lonrange, st_yr, do_parfor, do_clean_files, do_continue, stnnums, debug_rg)

WHAT IS THIS PROGRAM, IAN?
Looks like it runs the GFDL Perfect Model data for the locations given in a test site file.
%   This version tests a model for downscaling on 25 test sites.
%       By default will run tasmin & tasmax for cm3 and ecm, for the entire globe
%       
%   Inputs:
%       varnames    tasmin, tasmax, pr or both.                         ["tasmin","tasmax"]
%       e_or_cs     cm3, ecm or both.                                   ["cm3","ecm"]
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
%
%
        % % lats for lubbock, TX region
        % testlat = 33.576;
        % testlon = 258.145;

        % lat/lon for Denver region
        % testlat = 40;
        % testlon = 255;
        
        % lat/lon for central US test region
        
            % defaults should run grid_N25W100.nc , which has differences in row 1, columns 1-4.
    if (isempty(varnames)), varnames="pr"; end
    if (isempty(e_or_cs)),  e_or_cs  = "cm3";  end
	if (~exist("stnnums","var")), stnnums=[];end
    
    

%
%   Edit as needed!  Generally will need to change obsdir and obsnames below for other gridded obs data.
%
%   if (~exist("varnames","var")),                                      varnames  = ["tasmin","tasmax"];  end
%   if (~exist("e_or_cs","var")),                                       e_or_cs = ["cm3","ecm"];          end
%   if (~exist("dry_run","var")          || isempty(dry_run)),          dry_run = false;                  end
    if (~exist("latrange","var")),                                      latrange = [25,30];               end
    if (~exist("lonrange","var")),                                      lonrange = [-100,-95];            end
    if (~exist("st_yr","var")            || isempty(st_yr)),            st_yr = 2086;                     end
    if (~exist("do_parfor","var")        || isempty(do_parfor)),        do_parfor = false;                end
    if (~exist("do_continue","var")      || isempty(do_continue)),      do_continue = false;              end
    if (~exist("do_clean_files","var")),                                do_clean_files = false;           end
    if (~exist("debug_rg","var")),                                      debug_rg = false;                 end
        
    if (ischar_s(dry_run))
        if (strcmp(dry_run,"dry_run"))
            dry_run = true;
        elseif (strcmp(dry_run,"run"))
            dry_run = false;
        else
            error("dry_run (%s) should be logical or 'dry_run' or 'run'", dry_run);
        end
    end
    if (ischar_s(do_parfor))
        if (any(strcmp(do_parfor,["parfor","parallel"])))
            do_parfor = true;
        elseif (strcmp(do_parfor,"serial"))
            do_parfor = false;
        else
            error("do_parfor (%s) should be logical or 'parfor' or 'serial'", do_parfor);
        end
    end
     if (ischar_s(do_clean_files))
        if (strcmp(do_clean_files,"clean"))
            
            do_clean_files = true;
        elseif (strcmp(do_clean_files,"no_clean"))
            do_clean_files = false;
        else
            error("do_clean_files (%s) should be logical or 'clean' or 'no_clean'", do_clean_files);
        end
    end
     if (ischar_s(do_continue))
        if (strcmp(do_continue,"continue"))
            do_continue = true;
        elseif (strcmp(do_continue,"restart") || strcmp(do_continue, "start"))
            do_continue = false;
        else
            error("do_continue (%s) should be logical or 'continue' or 'restart'", do_continue);
        end
     end
    
     %-----  main work starts here
     
    interp_method = "bilinear";
    testsites=strings(0);
    
    if (ischar_s(latrange))
        if (strcmpi(latrange,"global"))
            latrange = [-90,90];
            lonrange = [0,360];
            region = "global";
        elseif (strcmpi(latrange,"conus"))
            latrange = [20,50];
            lonrange = [-130,-55];
            region = "conus";
        elseif (strcmpi(latrange,"central"))
            latrange = [28,48];
            lonrange = [-120,-90];
            region = "central";
        elseif (strcmpi(latrange,"test"))
            latrange = [32,40];
            lonrange = [-102,-94];
            region = "test_region";
        elseif (strcmpi(latrange, "test_25"))
            load("test_stations_25.mat","siteTbl");
            if (~isempty(stnnums))
                siteTbl = siteTbl(stnnums,:);
            end
            latrange = [];
            lonrange = [];
            lats = siteTbl.lat;
            lons = siteTbl.lon;
            region = "gfdl_test_stations";
            testsites=string(siteTbl.stnName);
        else
            error("error:  bad latrange: %s", latrange);
        end
    elseif ((all(latrange==[-90,90]) && (all(lonrange==[0,360]) || all(lonrange==[-180,180]))) || (isempty(latrange) && isempty(lonrange)))
        region = "global";
    else
        region = "misc_run";
    end
    
    if (isfolder("/Volumes/lacie_1/"))
        basesrc = "/Volumes/lacie_1/data/gcm_cmip5/hires/ian_work";
        baseout = fullfile("/Volumes/lacie_1/data/downscaled/arrm_v2/GFDL_perfect_model/test",region);
            
    elseif (isfolder("/lustre/research"))
        basesrc = "/lustre/research/hayhoe/gfdl_hires_ian/joined";
        baseout = fullfile("/lustre/scratch/iscottfl/downscaled/GFDL_perfect_model/test",region);
    end
        
    for vv=1:length(varnames)
        varname = varnames(vv);
        for ee=1:length(e_or_cs)
            e_or_c = e_or_cs(ee);
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
            if (st_yr == 2086)
                mdl_yrs = [2086,2115];
                pdf_yrstep = 30;
                yr_lbl = "future";
                mdldir   = basesrc; %fullfile(basesrc, "future");
                mdlname  = sprintf("c360_hiram_%s_rcp85_%s.lowres.%s.20860101-21151231.llt.nc",e_or_c, Xstuff, varname);
            elseif (st_yr == 1979)
                mdl_yrs = [1979,2038];
                pdf_yrstep = 60;
                yr_lbl = "hist";
                mdldir   = basesrc; %fullfile(basesrc, "hist");
                mdlname  = sprintf("%s_day_GFDL-HIRAM-C360_amip_19790101-20381231.H2H3.lowres.llt.nc", varname);
            else
                error("bad start yr: %d.  must be 1979 or 2086");
            end
%           trend_yrs = mdl_yrs;
            trend_order = 1;

            if (strcmp(region,"gfdl_test_stations"))
                outname=sprintf("downscaled.%s.%s.%s.%s.gfdl_PM.%s.%4d.%4d.%%s.mat", model, ensemble, varname, scenario, region, mdl_yrs(1), mdl_yrs(2));
            else
                outname=sprintf("downscaled.%s.%s.%s.%s.gfdl_PM.%s.%4d.%4d.LLGRID.nc", model, ensemble, varname, scenario, region, mdl_yrs(1), mdl_yrs(2));
            end
    
            obsdir   = fullfile(basesrc);
            histdir  = fullfile(basesrc);
            if (strcmp(region,"gfdl_test_stations"))
                outsubdir = "downscaled_test";
            elseif (do_parfor)
                outsubdir = strcat(fullfile(yr_lbl, model,ensemble,varname,scenario),"_parallel");
            else
                outsubdir = fullfile(yr_lbl, model,ensemble,varname,scenario);
            end
            outdir   = fullfile(baseout, outsubdir);

            gridsize = [5,5];
                        
            if (dry_run)
                fprintf("dry-run:  run (%d,%d) parameters:\n", vv,ee);
                fprintf("model:         %s\n", model);
                fprintf("ensemble:      %s\n", ensemble);
                fprintf("varname:       %s\n", varname);
                fprintf("scenario:      %s\n", scenario);
                fprintf("region:        %s\n", region);
                fprintf("obsnames:      %s\n",  vec2string(obsname, "format","%s ","brackets",'[]'));
                fprintf("histnames:     %s\n", vec2string(histname,"format","%s ","brackets",'[]'));
                fprintf("mdlnames:      %s\n", vec2string(mdlname, "format","%s ","brackets",'[]'));
                fprintf("outname:       %s\n", outname);
                fprintf("obsdir:        %s\n", obsdir);
                fprintf("histdir:       %s\n", histdir);                
                fprintf("mdldir:        %s\n", mdldir);
                fprintf("outdir:        %s\n", outdir);
                fprintf("obs_yrs:       %s\n", vec2string(obs_yrs,  "format","%4d ","brackets",'[]'));
                fprintf("hist_yrs:      %s\n", vec2string(hist_yrs, "format","%4d ","brackets",'[]'));
                fprintf("mdl_yrs:       %s\n", vec2string(mdl_yrs,  "format","%4d ","brackets",'[]'));
                if (~isempty(testsites))
                    fprintf("lats:          %s\n", vec2string(lats, "format","%9.4f ", "brackets",'[]'));
                    fprintf("lons:          %s\n", vec2string(lons, "format","%9.4f ", "brackets",'[]'));
                    fprintf("sites:         %s\n", vec2string([testsites(1:min(5,length(testsites)));"..."]));
                else
                    fprintf("latrange:      %s\n", vec2string(latrange, "format","%9.4f ", "brackets",'[]'));
                    fprintf("lonrange:      %s\n", vec2string(lonrange, "format","%9.4f ", "brackets",'[]'));
                    fprintf("gridsize:      %s\n", vec2string(gridsize, "format","%2d ",   "brackets",'[]'));
                    fprintf("do_continue:   %d\n", do_continue);
                end
                fprintf("interp_method: %s\n", interp_method);
%               fprintf("trend_yrs:     %s\n", vec2string(trend_yrs,"format","%4d ","brackets",'[]'));
                fprintf("pdf_yrstep:    %d\n", pdf_yrstep);
                fprintf("trend_order:   %d\n", trend_order);
                fprintf("do_parfor:     %d\n", do_parfor);
                fprintf("debug_rg:      %d\n", debug_rg);
                fprintf("do_clean_files %d\n", do_clean_files);
                fprintf("\n");
                
                if (~isfile(fullfile(obsdir,  obsname))),  fprintf("missing  obs file: %s\n", fullfile(obsdir, obsname)); end
                if (~isfile(fullfile(histdir, histname))), fprintf("missing hist file: %s\n", fullfile(histdir, histname)); end
                if (~isfile(fullfile(mdldir, mdlname))), fprintf("missing  mdl file: %s\n", fullfile(mdldir, mdlname)); end
                if (~isfolder(baseout)), fprintf("missing output folder: %s\n", baseout); end
                fprintf("\n\n");
            else   
                if (do_parfor)
                    ARRM_V2_start_workers();
                end  

                setup_dirs(baseout, outsubdir)
                if (do_clean_files)
                    clean_dirs(baseout, outsubdir);
                end
                if (strcmp(region,"gfdl_test_stations"))
                     ARRM_V2_wrapper(model,ensemble,varname,scenario,obs_yrs,hist_yrs,mdl_yrs,...
                                     "obsnames",    obsname, ...
                                     "histnames",   histname, ...
                                     "mdlnames",    mdlname, ...
                                     "outname",     outname, ...
                                     "region",      region, ...
                                     "testsites",   testsites, ...
                                     "obsdir",      obsdir,...
                                     "histdir",     histdir, ...
                                     "mdldir",      mdldir,...
                                     "outdir",      outdir, ...
                                     "lats",    lats, "lons",lons, ...
                                     "llgrid_size", gridsize, ...
                                     "do_continue", false, ...
                                     "interp_method",interp_method, ...
                                     "pdf_yrstep",  pdf_yrstep, ...     % GFDL-PM-specific!
                                     "pdf_yrs",     pdf_yrstep, ...     % GFDM-PM-specific.  either doing all hist as 1 period, or all future as 1 period.
                                     "trend_order", trend_order, ...    % GFDL-PM-specific!
                                     "do_parfor",false,"exit_on_error",true,"extended",false ...                         
                                     );
                else
%                                  "debug_rg",    debug_rg, ...
%                                      "trend_yrs",   trend_yrs, ...      % don't need to pass trend_yrs in.  ARRM_V2_run_gridded
%                                                                     restricts the trend_yrs to the available data.
                     ARRM_V2_wrapper(model,ensemble,varname,scenario,obs_yrs,hist_yrs,mdl_yrs,...
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
                                     "pdf_yrstep",  pdf_yrstep, ...     % GFDL-PM-specific!
                                     "trend_order", trend_order, ...    % GFDL-PM-specific!
                                     "do_parfor",false,"exit_on_error",true,"extended",false ...                         
                                     );
                end
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
