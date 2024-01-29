function gfdl_gridded_downscaling(varnames, e_or_cs, dry_run, latrange, lonrange, st_yr, do_parfor, do_clean_files, do_continue, pdf_map_method, nprocesses, partition)%, do_norm_sigmas , debug_rg)

This was rewritten as gfdl_test_downscaling, which also supports giving it a list of lat/lon locations for testing.

%   setups up and runs ARRM_V2 for GFDL Perfect Model Data runs
%       By default will run tasmin & tasmax for cm3 and ecm, for the entire globe
%       
%   Inputs:
%       varnames        tasmin, tasmax or both.  Or pr                      ["tasmin","tasmax"]
%       e_or_cs         cm3, ecm or both.                                   ["cm3","ecm"]
%       dry_run         true or false.  
%       latrange        lat & lon range to run, as 2-element vector.        [ -90,  90]
%       lonrange                                                        [   0, 360]
%                           (lon can be in range -180 - 180 or 0 - 360 
%                       shortcut for CONUS:  set latrange and lonrange to "conus" to do continental US only.
%                           "conus":  lat 20 - 50, lon 230 - 305, which gets parts of mexico & canada as well.
%                           "central":
%                           "test"
%       st_yr           temp:  [2086]  pr:  [2066]  No reason to change this...       
%       do_parfor       true or false.  controls whether to run parallel.   [true]
%       do_clean_files  [true]  if true, deletes any existing output files before it starts
%       do_continue                                                     [false]
%                       false:  run all gridboxes, regardless of previous completion status.
%                       true:   look in output folder, checks all existing grid_*.nc's for completion_flag, and runs
%                               all missing gridboxes and reruns any gridboxes where the completion_status is not 1 or 2
%                                   completion_status of 1:     complete
%                                   completion_status of 2:     complete, but all gridcells had insufficient data
%                                               only happens for livneh and similar inputs where water gridcells are all NAs
%       pdf_map_method  how to compress the pdf's.  default:  'linear' for tasmin/tasmax, 'clim' for precip.
%       nprocesses      # of workers to start.  default:  max possible on specific system.  [8 icsf-kmac, 18 neys, 120 HPCC nocona/redraider] 
%       partition       'local' (neys, icsf-kmac) or 'redraider R2020b' (option for HPCC) 
%       do_parfor       [true]  whether or not to run in parallel.  
%       do_continue    [false]  
% %        debug_rg    boolean.  If true, adds additional debug info to log files.     [false]
%
%
%   Currently set up to use livneh data, and downscale for CONUS.
%
        % % lats for lubbock, TX region
        % testlat = 33.576;
        % testlon = 258.145;

        % lat/lon for Denver region
        % testlat = 40;
        % testlon = 255;
        
        % lat/lon for central US test region



%
%   Edit as needed!  Generally will need to change obsdir and obsnames below for other gridded obs data.
%
%   if (~exist("varnames","var")),                            varnames  = ["tasmin","tasmax"];  end
%   if (~exist("e_or_cs","var")),                             e_or_cs = ["cm3","ecm"];          end
%   if (~exist("dry_run","var")      || isempty(dry_run)),    dry_run = false;                  end
    if (~exist("latrange","var")),                            latrange = [];                    end
    if (~exist("lonrange","var")),                            lonrange = [];                    end
    if (~exist("st_yr","var")       || isempty(st_yr)),       st_yr = [];                       end
    if (~exist("do_parfor","var")   || isempty(do_parfor)),   do_parfor = true;                 end
    if (~exist("do_clean_files","var")),                      do_clean_files = true;            end
    if (~exist("do_continue","var") || isempty(do_continue)), do_continue = false;              end
    if (~exist("pdf_map_method","var")),                      pdf_map_method = strings(0);      end
    if (~exist("nprocesses","var")),                          nprocesses = [];                  end
    if (~exist("partition","var") || isempty(partition)),     partition = "local";              end
%   if (~exist("do_norm_sigmas","var") || isempty(do_norm_sigmas)), do_norm_sigmas = false;     end
%   if (~exist("debug_rg","var")),                            debug_rg = false;                 end

    [hostname, on] = get_hostname(true);
        
    if (ischar_s(dry_run))
        if (any(strcmp(dry_run,["dry_run","dry-run"])))
            dry_run = true;
        elseif (strcmp(dry_run,"run"))
            dry_run = false;
        else
            error("dry_run (%s) should be logical or 'dry_run' or 'run'", dry_run);
        end
    end
    if (ischar_s(do_parfor))
        if (strcmp(do_parfor,"parfor"))
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
            error("do_clean_files (%s) should be logical or 'clean' or 'no_clean'", do_clean);
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
     
     if (strcmp(varnames(1),"pr"))
         doing_pr = true;
     else
         doing_pr = false;
     end
     
     if (~exist("st_yr","var") || isempty(st_yr))
         if (doing_pr)
             st_yr = 2066;
         else
             st_yr = 2086;
         end
     end
     
     if (isempty(pdf_map_method))
         if (doing_pr)
             pdf_map_method = "clim";
         else
             pdf_map_method = "linear";
         end
     end
     if (~any(strncmp(pdf_map_method,["lin","cos","cli"],3)))
         error("pdf_map_method ( %s ) should be 'linear','cosine', or 'clim'"); 
     end
     
     if (on.dev_system)
         if (~strcmp(partition,"local")), fprintf("running on dev system;  changing partition to ""local""\n"); end
         partition = "local";
     end
    
     %-----  main work starts here
     
    interp_method = "bilinear";
        
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
        else
            error("error:  bad latrange: %s", latrange);
        end
    elseif ((all(latrange==[-90,90]) && (all(lonrange==[0,360]) || all(lonrange==[-180,180]))) || (isempty(latrange) && isempty(lonrange)))
        region = "global";
    else
        region = "misc_run";
    end
    
    if (isfolder("/Volumes/lacie_1/"))
        basesrc = "/Volumes/lacie_1/data/gfdl_hires";
        baseout = fullfile("/Volumes/lacie_1/data/downscaled/arrm_v2/GFDL_perfect_model",region);
        
            
    elseif (isfolder("/lustre/research"))
        basesrc = "/lustre/research/hayhoe/gfdl_hires_ian/joined";
        baseout = fullfile("/lustre/scratch/iscottfl/downscaled/GFDL_perfect_model_test",region);
    else
        dry_run = true;
        fprintf(2,"oops.  Don't know how to run on %s.  Switching to dry_run\n", hostname)
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
                mdldir   = fullfile(basesrc, "future");
                mdlname  = sprintf("c360_hiram_%s_rcp85_%s.lowres.%s.20860101-21151231.llt.nc",e_or_c, Xstuff, varname);
            elseif (st_yr == 2066)
                mdl_yrs = [2066,2095];
                pdf_yrstep = 30;
                yr_lbl = "future";
                mdldir   = fullfile(basesrc, "future");
                mdlname  = sprintf("c360_hiram_%s_rcp85_%s.lowres.%s.20660101-20951231.llt.nc",e_or_c, Xstuff, varname);
            elseif (st_yr == 1979)
                mdl_yrs = [1979,2038];
                pdf_yrstep = 60;
                yr_lbl = "hist";
                mdldir   = fullfile(basesrc, "hist");
                mdlname  = sprintf("%s_day_GFDL-HIRAM-C360_amip_19790101-20381231.H2H3.lowres.llt.nc", varname);
            else
                error("bad start yr: %d.  must be 1979, 2066 or 2086");
            end
%           trend_yrs = mdl_yrs;
            trend_order = 1;

            outname=sprintf("downscaled.%s.%s.%s.%s.gfdl_PM.%s.%s.%4d.%4d.LLGRID.nc", model, ensemble, varname, scenario, region, pdf_map_method,mdl_yrs(1), mdl_yrs(2));
    
            obsdir   = fullfile(basesrc, "hist");
            histdir  = fullfile(basesrc, "hist");
            outsubdir= sprintf("%s_%s_%s_%s_%s_%s",yr_lbl, model,ensemble,varname,scenario, pdf_map_method);
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
                fprintf("latrange:      %s\n", vec2string(latrange, "format","%9.4f ", "brackets",'[]'));
                fprintf("lonrange:      %s\n", vec2string(lonrange, "format","%9.4f ", "brackets",'[]'));
                fprintf("gridsize:      %s\n", vec2string(gridsize, "format","%2d ",   "brackets",'[]'));
                fprintf("do_continue:   %d\n", do_continue);
                fprintf("interp_method: %s\n", interp_method);
%               fprintf("trend_yrs:     %s\n", vec2string(trend_yrs,"format","%4d ","brackets",'[]'));
                fprintf("pdf_yrstep:    %d\n", pdf_yrstep);
                fprintf("trend_order:   %d\n", trend_order);
                fprintf("do_parfor:     %d\n", do_parfor);
%               fprintf("debug_rg:      %d\n", debug_rg);
                fprintf("do_clean_files %d\n", do_clean_files);
%               fprintf("do_norm_sigmas %d\n", do_norm_sigmas);
                fprintf("pdf_map_method %s\n", pdf_map_method);
                if (do_parfor)
                    if (isempty(nprocesses))
                        fprintf("nprocesses:      (max possible)\n")
                    else
                        fprintf("nprocesses:      %d\n", nprocesses);
                    end
                    fprintf("partition:     %s\n", partition);
                end                
                fprintf("\n");
                
                if (~isfile(fullfile(obsdir,  obsname))),  fprintf("missing  obs file: %s\n", fullfile(obsdir, obsname)); end
                if (~isfile(fullfile(histdir, histname))), fprintf("missing hist file: %s\n", fullfile(histdir, histname)); end
                if (~isfile(fullfile(mdldir, mdlname))), fprintf("missing  mdl file: %s\n", fullfile(mdldir, mdlname)); end
                if (~isfolder(baseout)), fprintf("missing output folder: %s\n", baseout); end
                fprintf("\n\n");
            else   
                if (do_parfor)
                    ARRM_V2_start_workers(nprocesses, partition);
                end  

                setup_dirs(baseout, outsubdir)
                if (do_clean_files)
                    clean_dirs(baseout, outsubdir);
                end
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
                                 "pdf_map_method",  pdf_map_method, ...
                                 "pdf_yrstep",  pdf_yrstep, ...     % GFDL-PM-specific!
                                 "trend_order", trend_order, ...    % GFDL-PM-specific!
                                 "do_parfor",do_parfor,"exit_on_error",true,"extended",false ...                         
                                 );
%                                "do_normalize_sigmas", do_norm_sigmas, ...
%                                "debug_rg",debug_rg, ...
%                                "trend_yrs",   trend_yrs, ...      % don't need to pass trend_yrs in.  ARRM_V2_run_gridded
%                                                                     restricts the trend_yrs to the available data.
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
