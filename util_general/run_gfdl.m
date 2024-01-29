function retval = run_gfdl(varname, e_or_c, time_pd, region, ncores, baseout, sigma_normalize, prcp_distrib, dry_run, do_parfor, plotflag, far_outlier_thresh, cdf_append_pts, far_outlier_anchor_pt, do_continue)

    if (~any(strcmp(varname, ["pr","tasmax","tasmin"]))), error("error:  bad varname: %s\n", varname); end
    if (~any(strcmp(e_or_c,  ["esm","cm3"]           ))), error("error:  bad e_or_c:  %s\n", e_or_c);  end
%    if (~any(strcmp(region,  ["conus","global"]      ))), error("error:  bad region: %s\n",  region);  end
    if (isnumeric(region))
        if (numel(region)==2)
            region = [region(1), region(1)+1; region(2), region(2)+1];
        elseif (~all(size(region)==[2,2]))
            error("error:  numeric region must be size (2,2), [latrange; lonrange]"); 
        end
        latrange=region(1,:);
        lonrange=region(2,:);
        region=strings(0);
    else
        if (~any(strcmp(time_pd, ["hist","future"]       ))), error("error:  bad time_pd: %s\n", time_pd); end
        latrange=[];
        lonrange=[];
    end
    
            % Defaults:
    if (~exist("ncores","var")), error("error:  must specify ncores;  can be []"); end
    if (~exist("baseout","var")), error("error: must specify baseout as [] for defaults"); end
    if (~exist("sigma_normalize","var")), error("error: must specify sigma_normalize (as true or false)"); end
    if (~exist("prcp_distrib","var"))
        if (strcmp(varname,"pr")), error("error:  must specify prcp_distrib for precip run"); end
    else
        prcp_distrib = [];
    end
    if (~exist("dry_run","var") || isempty(dry_run)), dry_run = true; end
    if (~exist("do_parfor","var")), do_parfor = false; end
    if (~exist("plotflag","var") || isempty(plotflag))             % plotflag s/b 3 to actually draw the plots.
        if(~do_parfor), plotflag=2; else, plotflag=0; end
    end 
    if (~exist("far_outlier_thresh", "var"))
        if (is_prcp_variable(varname))
            error("far_outlier thresh must be set for precip");
        else
            far_outlier_thresh = 1e-4;
        end
    end
    
    if (~exist("cdf_append_pts","var") || isempty(cdf_append_pts))
        cdf_append_pts = 0;
    end
    
    if (~exist("far_outlier_anchor_pt", "var") || isempty(far_outlier_anchor_pt))
        far_outlier_anchor_pt = far_outlier_thresh;
    end
    
    if (~exist("do_continue","var") || isempty(do_continue)), do_continue = false; end
        
    if (isempty(ncores))
        if (~do_parfor)
            ncores = 1;
        else
            myhost = get_hostname(true);
            if (strcmp(myhost,"neys"))
                ncores = 18;
            else
                ncores = 96;        % user must make sure they've got 96 cores on HPCC!
            end
        end
    end
            
    
    if (~isempty(baseout))
        if (extractBetween(baseout,1,1) ~= "/")
            baseout = fullfile("/lustre/scratch/iscottfl/downscaled/", baseout);
            if (~isfolder(baseout)), mkdir(baseout); end
        end
    end
    
    [~, on] = get_hostname(true);

    if (on.hpcc_system)
        nworkers = ncores;
        use_lacie = false;
    elseif (on.neys || on.icsf_kmac)
        use_lacie = true;
        nworkers = [];
    end

    if (strncmp(time_pd, "hist",4))
        start_yr = 1979;
    elseif (strncmp(time_pd, "fut",3))
        if (strcmp(varname,"pr"))
            start_yr = 2066;
        else
            start_yr = 2086;
        end
    end
    retval = gfdl_test_downscaling(varname,e_or_c,start_yr,"dry_run", dry_run,"do_parfor", do_parfor,  "baseout", baseout, ...
                                   "latrange",latrange,"lonrange", lonrange,"region",region, ...
                                   "sigma_normalize", sigma_normalize,"median_normalize",false, "prcp_distrib", prcp_distrib, ...
                                   "far_outlier_thresh", far_outlier_thresh , "cdf_append_pts", cdf_append_pts, "far_outlier_anchor_pt", far_outlier_anchor_pt, ...
                                   "do_compare", false, "show_figs",false, "do_continue", do_continue, "do_clean_files", ~do_continue, ...
                                   "use_lacie", use_lacie, "cluster", "local", "nworkers", nworkers, "internal_plotflag",plotflag ...
                                   );
end

% gfdl_test_downscaling("tasmax","esm",1979,"dry_run", false,"do_parfor", true,  "latrange",[],"lonrange", [],"region","conus", "do_compare", false, "show_figs",false, "do_clean_files", true, "use_lacie", true, "cluster", "local", "nworkers", nworkers);
% gfdl_test_downscaling("tasmax","cm3",1979,"dry_run", false,"do_parfor", true,  "latrange",[],"lonrange", [],"region","conus", "do_compare", false, "show_figs",false, "do_clean_files", true, "use_lacie", true, "cluster", "local", "nworkers", nworkers);
% gfdl_test_downscaling("tasmin","esm",1979,"dry_run", false,"do_parfor", true,  "latrange",[],"lonrange", [],"region","conus", "do_compare", false, "show_figs",false, "do_clean_files", true, "use_lacie", true, "cluster", "local", "nworkers", nworkers);
% gfdl_test_downscaling("tasmin","cm3",1979,"dry_run", false,"do_parfor", true,  "latrange",[],"lonrange", [],"region","conus", "do_compare", false, "show_figs",false, "do_clean_files", true, "use_lacie", true, "cluster", "local", "nworkers", nworkers);
% gfdl_test_downscaling("pr","esm",1979,"dry_run", false,"do_parfor", true,  "latrange",[],"lonrange", [],"region","conus", "do_compare", false, "show_figs",false, "do_clean_files", true, "use_lacie", true, "cluster", "local", "nworkers", nworkers);
% gfdl_test_downscaling("pr","cm3",1979,"dry_run", false,"do_parfor", true,  "latrange",[],"lonrange", [],"region","conus", "do_compare", false, "show_figs",false, "do_clean_files", true, "use_lacie", true, "cluster", "local", "nworkers", nworkers);
% 
% gfdl_test_downscaling("tasmax","esm",2086,"dry_run", false,"do_parfor", true,  "latrange",[],"lonrange", [],"region","conus", "do_compare", false, "show_figs",false, "do_clean_files", true, "use_lacie", true, "cluster", "local", "nworkers", nworkers);
% gfdl_test_downscaling("tasmax","cm3",2086,"dry_run", false,"do_parfor", true,  "latrange",[],"lonrange", [],"region","conus", "do_compare", false, "show_figs",false, "do_clean_files", true, "use_lacie", true, "cluster", "local", "nworkers", nworkers);
% gfdl_test_downscaling("tasmin","esm",2086,"dry_run", false,"do_parfor", true,  "latrange",[],"lonrange", [],"region","conus", "do_compare", false, "show_figs",false, "do_clean_files", true, "use_lacie", true, "cluster", "local", "nworkers", nworkers);
% gfdl_test_downscaling("tasmin","cm3",2086,"dry_run", false,"do_parfor", true,  "latrange",[],"lonrange", [],"region","conus", "do_compare", false, "show_figs",false, "do_clean_files", true, "use_lacie", true, "cluster", "local", "nworkers", nworkers);
% gfdl_test_downscaling("pr","esm",2066,"dry_run", false,"do_parfor", true,  "latrange",[],"lonrange", [],"region","conus", "do_compare", false, "show_figs",false, "do_clean_files", true, "use_lacie", true, "cluster", "local", "nworkers", nworkers);
% gfdl_test_downscaling("pr","cm3",2066,"dry_run", false,"do_parfor", true,  "latrange",[],"lonrange", [],"region","conus", "do_compare", false, "show_figs",false, "do_clean_files", true, "use_lacie", true, "cluster", "local", "nworkers", nworkers);

