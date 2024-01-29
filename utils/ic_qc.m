function ic_qc(obs_src, mdl, scenario, grgn, latrng, lonrng, daterange, figdir, show_figs, part_figs, use_tllc_time)
%
%   create output map arrays
%       QC5:  tmax > max(obs_tmax)
%       QC6:  tmax < min(obs_tmax)
%       QC7:  tmin > max(obs_tmin)
%       QC8:  tmin < min(obs_tmin)
%       QC18: mean(mdl) - mean(obs)
%
%   fn_obs s/b LLT order
%   read lats,lons from fn_obs
%   parfor each lon:
%       read all obs lats, all days 
%           (convert to 365-day calendar)  not an issue for livneh, which is already 365-day for our data
%       read all mdl lats, all days
%       for each lat:
%           find values for QC maps, put in output arrays.
%       end
%   end
%
%   draw maps for R-version's QC's 1-17, plus avg diff & # NA's in Obs.
%
%   Need to customize the colormaps.
%   should add map of std_dev, skewness & kurtosis.
%
%   Ian Scott-Fleming, 10/2021
%
%--------------------------

    start_tic = tic;
    [hostname, on] = get_hostname();
    
    fprintf("running on host %s, obs_src %s, mdl %s, scenario %s\n", hostname, obs_src, mdl, scenario);
    
    if (on.hpcc_system)
        if any(strcmp(obs_src,["livneh","conus_livneh"]))
            obs_src_dir = "/lustre/scratch/iscottfl/obs/livneh";
            if (use_tllc_time)
                mdl_src_dir = fullfile("/lustre/scratch/iscottfl/downscaled/cmip6/conus_livneh",scenario,"tllc_time");
            else
                mdl_src_dir = fullfile("/lustre/scratch/iscottfl/downscaled/cmip6/conus_livneh",scenario);
            end                
        elseif (strcmp(obs_src,"nclimgrid"))
            obs_src_dir = "/lustre/scratch/iscottfl/obs/nclimgrid";
            if (use_tllc_time)
                mdl_src_dir = fullfile("/lustre/scratch/iscottfl/downscaled/cmip6/nclimgrid",scenario,"tllc_time");
            else
                mdl_src_dir = fullfile("/lustre/scratch/iscottfl/downscaled/cmip6/nclimgrid",scenario);
            end                
        else
            error("unknown src.  Must be livneh or nclimgrid");
        end
    else
        if any(strcmp(obs_src,["livneh","conus_livneh"]))
            obs_src_dir = "/Volumes/lacie_1/data/obs/livneh";
        elseif (strcmp(obs_src,"nclimgrid"))
            obs_src_dir = "/Volumes/lacie_1/data/obs/nclimgrid";
        else
            error("unknown src.  Must be livneh or nclimgrid");
        end
        mdl_src_dir = "/Volumes/lacie_1/data/work";
    end
    
    if any(strcmp(obs_src,["livneh","conus_livneh"]))
        fn_obs_max = "livneh_1_16th_Tmax.1950.2015.llt.nc";
        fn_obs_min = "livneh_1_16th_Tmin.1950.2015.llt.nc";
        fn_mdl_max = sprintf("downscaled.%s.r1i1p1f1.tasmax.%s.%s.livneh_1_16th.conus_livneh.1950.2100.tllc_time.nc", mdl, scenario, grgn);
        fn_mdl_min = sprintf("downscaled.%s.r1i1p1f1.tasmin.%s.%s.livneh_1_16th.conus_livneh.1950.2100.tllc_time.nc", mdl, scenario, grgn);
        obs_max_vname = "Tmax";
        obs_min_vname = "Tmin";
    else
        fn_obs_max = "nclimgrid.day.tmax.1951.2019.llt.nc";
        fn_obs_min = "nclimgrid.day.tmin.1951.2019.llt.nc";
        if (use_tllc_time)
            fn_mdl_max = sprintf("downscaled.%s.r1i1p1f1.tasmax.%s.%s.nclimgrid.nclimgrid.1950.2100.tllc_time.nc", mdl, scenario, grgn);
            fn_mdl_min = sprintf("downscaled.%s.r1i1p1f1.tasmin.%s.%s.nclimgrid.nclimgrid.1950.2100.tllc_time.nc", mdl, scenario, grgn);
        else
            fn_mdl_max = sprintf("downscaled.%s.r1i1p1f1.tasmax.%s.%s.nclimgrid.nclimgrid.1950.2100.nc", mdl, scenario, grgn);
            fn_mdl_min = sprintf("downscaled.%s.r1i1p1f1.tasmin.%s.%s.nclimgrid.nclimgrid.1950.2100.nc", mdl, scenario, grgn);
        end
        obs_max_vname = "tmax";
        obs_min_vname = "tmin";
    end
    
    if (~exist("daterange","var") || isempty(daterange))
        daterange = [1950,1,1;2000,12,31];
    end
    
    if (~exist("show_figs","var") || isempty(show_figs))
        show_figs = true;
    end
    
    if (~exist("part_figs","var") || isempty(part_figs))
        part_figs = show_figs;
    end
    
    mdl_max_vname = "tasmax";
    mdl_min_vname = "tasmin";

    fprintf("obs_src_dir: %s\n", obs_src_dir);
    fprintf("mdl_src_dir: %s\n", mdl_src_dir);
    fprintf("fn_obs_max:  %s\n", fn_obs_max);
    fprintf("fn_obs_min:  %s\n", fn_obs_min);
    fprintf("fn_mdl_max:  %s\n", fn_mdl_max);
    fprintf("fn_mdl_min:  %s\n", fn_mdl_min);
    ls("-l",fullfile(obs_src_dir,fn_obs_max),fullfile(obs_src_dir,fn_obs_min), fullfile(mdl_src_dir,fn_mdl_max),fullfile(mdl_src_dir,fn_mdl_min));
        
    nc_obs_max = ncdf(fullfile(obs_src_dir, fn_obs_max));
    nc_obs_min = ncdf(fullfile(obs_src_dir, fn_obs_min));
    nc_mdl_max = ncdf(fullfile(mdl_src_dir, fn_mdl_max));
    nc_mdl_min = ncdf(fullfile(mdl_src_dir, fn_mdl_min));
    
    [olats,olons] = ncdf_get_latlons(nc_obs_max);
    [mlats,mlons] = ncdf_get_latlons(nc_mdl_max);
        % make sure lats & lons are doubles.  nclimgrid uses single-precision, but geoRefCells in make_refmat needs doubles.
    olats = double(olats);
    olons = double(olons);
    
        % make sure lons are in range [-180,180)
    olons = mod(olons+180,360)-180;
    mlons = mod(mlons+180,360)-180;
    
    if (exist('latrng','var') && ~isempty(latrng))
        olats = olats((olats>= latrng(1)-.0001) & (olats<= latrng(2)+.0001));
        mlats = mlats((mlats>= latrng(1)-.0001) & (mlats<= latrng(2)+.0001));
    end
    if (exist('lonrng','var') && ~isempty(lonrng))
        lonrng = mod(lonrng+180,360)-180;
        olons = olons((olons>= lonrng(1)-.0001) & (olons<= lonrng(2)+.0001));
        mlons = mlons((mlons>= lonrng(1)-.0001) & (mlons<= lonrng(2)+.0001));
    end
    if (~exist("figdir","var") || isempty(figdir))
        figdir = "./qc_figs";
    end
    matdir = fullfile(figdir,"mats");
    
        % # of lats & lons
    nolats = length(olats);
    nolons = length(olons);
    nmlats = length(mlats);
    nmlons = length(mlons);
    
    R = make_refmat(olats,olons, nolats,nolons, olats(8)-olats(7), olons(8)-olons(7));      % for drawing maps.
        
    
%         % # of days in model data.
%     mdl_vix = find(strcmp(nc_mdl_max.varlist(), mdl_max_vname),1);
%     v = nc_mdl_max.Variables(mdl_vix);
%     ndays_mdl = v.Size(3);
    
    
    fprintf("latrange:  %9.4f - %9.4f (%d lats)\n", olats(1),olats(end), nolats);
    fprintf("lonrange:  %9.4f - %9.4f (%d lons)\n", olons(1),olons(end), nolons);
    
    if (any(abs(olats-mlats)>.0001) || any(abs(olons-mlons)>.0001) || nolats ~= nmlats || nolons ~= nmlons)
        error("lat/lon mismatch");
    end
    
    olatrange = [olats(1), olats(end)];
%   olonrange = [olons(1), olons(end)];
    mlatrange = [mlats(1), mlats(end)];
%   mlonrange = [mlons(1), mlons(end)];
    
    nqcs = 26;
    qc = nan(nqcs,nolats,nolons);

        % setup dirs for figures.
    if (~isfolder(figdir)), error("error:  output folder %s does not exist", figdir); end
    sdir = sprintf("%s_%s", obs_src, mdl);
    setup_dirs(figdir,["mats", fullfile(sdir,"pngs"),fullfile(sdir,"figs")]);
    figdir = fullfile(figdir,sdir);
    if (~isfolder(figdir)), error("error:  output folder %s does not exist", figdir); end
    
    reads_time = zeros(1,nolons);
    calcs_time = zeros(1,nolons); 
    draw_time = 0;
    save_time = 0;
    yr_range = [ 10*365,  40*365-1; ...  % 1960-1989
                 70*365, 100*365-1; ...  % 2020-2049
                100*365, 130*365-1; ... % 2050-2079
                130*365, 150*365-1; ... % 2080-2099
                ];
    doyr = [1,59; 182, 243]; % jan/feb & jul/aug
    
    if (nolons <= 100)
        jmod = ceil(nolons/20);
        figmod = 5*jmod;
    elseif (nolons <= 400)
        jmod = ceil(nolons/100);
        figmod = 5*jmod;
    else
        jmod = 3;
        figmod = 3*jmod;
    end
    loop_tic = tic;
    t_end = 51*365; % 51 years, 18615 days
    fid = fopen("ic_qc.log","w");
    qc_ttl = strings(1,nqcs);
    qc_ttl(1)  = sprintf("QC1,  # of days Tmax(mdl) >  60 C 1950-2100 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(2)  = sprintf("QC2,  # of days Tmax(mdl) < -40 C 1950-2100 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(3)  = sprintf("QC3,  # of days Tmin(mdl) >  50 C 1950-2100 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(4)  = sprintf("QC4,  # of days Tmin(mdl) > -60 C 1950-2100 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(5)  = sprintf("QC5,  # of days Tmax(mdl) > max(Tmax(obs))+1 1950-2000 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(6)  = sprintf("QC6,  # of days Tmax(mdl) < min(Tmax(obs))-1 1950-2000 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(7)  = sprintf("QC7,  # of days Tmin(mdl) > max(Tmin(obs))+1 1950-2000 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(8)  = sprintf("QC8,  # of days Tmin(mdl) < min(Tmin(obs))-1 1950-2000 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(9)  = sprintf("QC9,  # of days Tmin(mdl) > Tmax(mdl) 1950-2100 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(10) = sprintf("QC10, # of NA's Tmax(mdl) 1950-2100 (thresh=16) %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(11) = sprintf("QC11, # of NA's Tmin(mdl) 1950-2100 (thresh=16) %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(12) = sprintf("QC12, Avg Tmax (future) < Tmax(hist) or Tmax increase > 5C %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(13) = sprintf("QC13, Avg Tmin (future) < Tmin(hist) or Tmin increase > 5C  %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(14) = sprintf("QC14, # of days abs(Tmax(i)-Tmax(i-1)) < 1e-5 1950-2100 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(15) = sprintf("QC15, # of days abs(Tmin(i)-Tmin(i-1)) < 1e-5 1950-2100 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(16) = sprintf("QC16, # of years mean(Tmax(jan/feb)) > mean(Tmax(jul/aug)) 1950-2100 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(17) = sprintf("QC17, # of years mean(Tmin(jan/feb)) > mean(Tmin(jul/aug)) 1950-2100 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(18) = sprintf("QC18, mean(Tmax(mdl)) - mean(Tmax(obs)) 1950-2000 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(19) = sprintf("QC19, mean(tmin(mdl)) - mean(tmin(obs)) 1950-2000 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(20) = sprintf("QC20, # of NA's Tmax(obs) 1950-2000 (thresh=1500) %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(21) = sprintf("QC21, # of NA's Tmin(obs) 1950-2000 (thresh=1500) %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(22) = sprintf("QC22, # of days Tmin(obs) > Tmax(obs) 1950-2015  %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(23)  = sprintf("QC9a, # of days Tmin(mdl) > Tmax(mdl) 1950-2100 (thresh=55) %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(24) = sprintf("QC22a, # of days Tmin(obs) > Tmax(obs) 1950-2015 (thresh=18) %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(25) = sprintf("QC25, # of NA's Tmax(obs) 1950-2015 (thresh=1500) %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(26) = sprintf("QC26, # of NA's Tmin(obs) 1950-2015 (thresh=1500) %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    
    
    for j=1:nolons
        
        reads_tic = tic;
        nco_max = ncdf_read_file(nc_obs_max, obs_max_vname, olatrange, olons(j), [], "365-day");
        ncm_max = ncdf_read_file(nc_mdl_max, mdl_max_vname, mlatrange, mlons(j), [], "365-day");
        nco_min = ncdf_read_file(nc_obs_min, obs_min_vname, olatrange, olons(j), [], "365-day");
        ncm_min = ncdf_read_file(nc_mdl_min, mdl_min_vname, mlatrange, mlons(j), [], "365-day");
        
        mdl_vix = find(strcmp(ncm_max.varlist(), mdl_max_vname),1);
        obs_vix = find(strcmp(nco_max.varlist(), obs_max_vname),1);
        
        reads_time(j) = toc(reads_tic);
        calcs_tic = tic;
        obs_tmax = nco_max.Variables(obs_vix).vdata;
        obs_tmin = nco_min.Variables(obs_vix).vdata;
        mdl_tmax = ncm_max.Variables(mdl_vix).vdata;
        mdl_tmin = ncm_min.Variables(mdl_vix).vdata;
        
        for i=1:nolats
            
             if (abs(olons(j) - (-113.7812)) < .001 && i == 30)
                fprintf("lat,lon = %.5f %.5f\n",olats(i), olons(j));
                fprintf("here\n");
             end
        
           mtmax = mdl_tmax(:,i);
            mtmin = mdl_tmin(:,i);
            otmax = obs_tmax(:,i);
            otmin = obs_tmin(:,i);
                % get the data for just 1950-2000
            mtmax2k = mdl_tmax(1:t_end,i);
            mtmin2k = mdl_tmin(1:t_end,i);
            otmax2k = obs_tmax(1:t_end,i);
            otmin2k = obs_tmin(1:t_end,i);
            
                % QC 12-13:  checking for increasing tmax & tmin over decades, and that mean value doesn't increase by > 5 C
            means_max = zeros(1,4);
            means_min = zeros(1,4);
            for k=1:4
                means_max(i) = mean(mdl_tmax(yr_range(k,1):yr_range(k,2)));
                means_min(i) = mean(mdl_tmin(yr_range(k,1):yr_range(k,2)));
            end
            qc(12,i,j) = means_max(2)>means_max(1) && means_max(3)>means_max(2) && means_max(4)>means_max(3) && means_max(4)-means_max(1) <= 5;
            qc(13,i,j) = means_min(2)>means_min(1) && means_min(3)>means_min(2) && means_min(4)>means_min(3) && means_min(4)-means_min(1) <= 5;
            
%             if (~any(isnan(means_max)) && ~any(isnan(means_min)))
%                 zz = true; %#ok<NASGU>
%             else
%                 zz = false; %#ok<NASGU>
%             end
%             
                % QC 14-15.  repetitive values
                % R version only counts if repeats >= 3x.  I'm counting any repeats w/in .0001 C
            difs_max = diff(mtmax);
            difs_min = diff(mtmin);
            rep_max = sum(abs(difs_max)<1e-5);
            rep_min = sum(abs(difs_min)<1e-5);
            qc(14,i,j) = rep_max;
            qc(15,i,j) = rep_min;
           
                % QC 16-17.  #years where mean Jan/Feb Tmax >  mean Jul/Aug Tmax 
            nyrs = size(mdl_tmax,1)/365;
            if (mod(nyrs,1) ~= 0)
                fprintf("oops\n");
            end
            max_3d = reshape(mtmax, 365,nyrs);
            min_3d = reshape(mtmin, 365,nyrs);
            qc(16:17,i,j)=0;
            for k=1:nyrs
                janfeb = mean(max_3d(doyr(1,1):doyr(1,2)),"all");
                julaug = mean(max_3d(doyr(2,1):doyr(2,2)),"all");
                qc(16,i,j) = qc(16,i,j) + (janfeb>julaug);
                janfeb = mean(min_3d(doyr(1,1):doyr(1,2)),"all");
                julaug = mean(min_3d(doyr(2,1):doyr(2,2)),"all");
                qc(17,i,j) = qc(17,i,j) + (janfeb>julaug);
            end
            
            qc(1,i,j)  = sum(mtmax  >  60);              % tmax >  60 C
            qc(2,i,j)  = sum(mtmax  < -40);              % tmax < -40 C
            qc(3,i,j)  = sum(mtmin  >  50);              % tmin >  50 C
            qc(4,i,j)  = sum(mtmin  < -60);              % tmin < -60 C
            qc(5,i,j)  =  sum(mtmax2k  > max(otmax2k)+1);  % # of days where mdl tmax > max observed tmax
            qc(6,i,j)  =  sum(mtmax2k  < min(otmax2k)-1); % # of days where mdl tmax < min observed tmax
            qc(7,i,j)  =  sum(mtmin2k  > max(otmin2k)+1); % # of days where mdl tmin > max observed tmin
            qc(8,i,j)  =  sum(mtmin2k  < min(otmin2k)-1); % # of days where mdl tmin < min observed tmin
            qc(9,i,j)  =  sum(mtmin  > mtmax);          % # of days where mdl tmin > mdl tmax
            qc(10,i,j) = min(16,sum(isnan(mtmax)));     % # of tmax NAs.  But using log scale so we can see if > 0.
            qc(11,i,j) = min(16,sum(isnan(mtmin)));     % # of tmin NAs.  This is different from R version.
            % qc12-17 are up above
            qc(18,i,j) = mean(mtmax2k) - mean(otmax2k);
            qc(19,i,j) = mean(mtmin2k) - mean(otmin2k);
            qc(20,i,j) = min(1500,sum(isnan(otmax2k)));
            qc(21,i,j) = min(1500,sum(isnan(otmin2k)));     
            qc(22,i,j) =  sum(otmin  > otmax);        % # of days where mdl tmin > mdl tmax
            qc(23,i,j)  =  min(55,sum(mtmin  > mtmax));        % # of days where mdl tmin > mdl tmax, truncated @ 20
            qc(24,i,j)  =  min(18,sum(otmin  > otmax));        % # of days where obs tmin > obs tmax, truncated @ 20
            qc(25,i,j) = min(1500,sum(isnan(otmax)));
            qc(26,i,j) = min(1500,sum(isnan(otmin)));            
            
                % set counts to NA if any NA's in mdl data
            if (qc(10,i,j) > 0)
                qc([1,2,5,6,9,12,14,16,18,23],i,j) = nan;
            end
            if (qc(11,i,j) > 0)
                qc([3,4,7,8,9,13,15,17,19,23],i,j) = nan;
            end
                
                % set counts to NA if any NA's in obs data
            if (qc(20,i,j) > 0)
                qc([5,6,22,24],i,j) = nan;
            end
            if (qc(21,i,j) > 0)
                qc([7,8,22,24],i,j) = nan;
            end
                
            
            
        end
        
        calcs_time(j) = toc(calcs_tic);
        
        fprintf(fid,"%d\n", j);
        if (mod(j,jmod)==0)
            [~,done] = system(" wc -l ic_qc.log | awk '{print $1}'");
            done = str2double(done);
            elapsed = toc(loop_tic);
            remain = elapsed/j*(nolons-j);
            fprintf("%4d of %4d complete  read_time %6.1f mins (%6.1f secs),  calc_time %8.3f secs.  elapsed: %6.1f mins remaining: %6.1f mins\n", done, nolons, sum(reads_time(1:j))/60, reads_time(j), sum(calcs_time(1:j)), elapsed/60, remain/60);
        end
        
%         if (j==floor(nolons/8))
%             fprintf("here\n");
%             fprintf("hello\n");
%         end
        
        cmap = parula(55);
        cmap(end,:) = [1,.25,.25];
        if (part_figs && mod(j, figmod)==0)
            draw_tic = tic;
            kqc=[5,6,23,24];
            for kk=1:4
                k=kqc(kk);
                min_cnt = 0;
                max_cnt = max(qc(k,:,:), [], "all","omitnan");
                if (kk==3)
                    display_map(qc(k,:,:), qc_ttl(k), R, 99, [2,2,kk], min_cnt, max_cnt, 2, cmap, [], [], [100 140 2000 1200], 10, true, ~part_figs);
                else
                    display_map(qc(k,:,:), qc_ttl(k), R, 99, [2,2,kk], min_cnt, max_cnt, 2, [], [], [], [100 140 2000 1200], 10, true, ~part_figs);
                end
            end
            draw_time = draw_time + toc(draw_tic);
                
        end
            
    end
            
%   draw maps

    h_nm   = strings(1,nqcs);
    max_cnt = zeros(1,nqcs);
    min_cnt = zeros(1,nqcs);
    for i=1:nqcs
        h_nm(i)  = sprintf("%s_%s_%s_QC%02d",  obs_src,mdl,scenario,i);
    end
    
%     h_nm(5)  = sprintf("%s_%s_%s_QC5.png",  obs_src,mdl,scenario);
%     h_nm(6)  = sprintf("%s_%s_%s_QC6.png",  obs_src,mdl,scenario);
%     h_nm(7)  = sprintf("%s_%s_%s_QC7.png",  obs_src,mdl,scenario);
%     h_nm(8)  = sprintf("%s_%s_%s_QC8.png",  obs_src,mdl,scenario);
%     h_nm(18) = sprintf("%s_%s_%s_QC18.png", obs_src,mdl,scenario);
%     h_nm(19) = sprintf("%s_%s_%s_QC19.png", obs_src,mdl,scenario);
    
    for i=1:nqcs
        max_cnt(i) = max(qc(i,:,:), [], "all","omitnan");
        if (i==18 || i==19)
            min_cnt(i) = min(qc(i,:,:), [], "all", "omitnan");
        else
            min_cnt(i) = 0;
        end
    end
        
    h = zeros(1,22);
    
    for i = 1:nqcs
        err=false;
        try
            draw_tic = tic;
            if any(i==[10,11])
                cmap = parula(16);
                cmap(end,:,:) = [1, .25,.25];
            elseif (i==23)
                cmap = parula(55);
                cmap(end,:,:) = [1,.25,.25];
            elseif (i==24)
                cmap=parula(18);
                cmap(end,:,:) = [1,.25,.25];
            else
                cmap=parula(256);
            end
            h(i) = display_map(qc(i,:,:), qc_ttl(i), R, i, [], min_cnt(i), max_cnt(i), 2, cmap, [], [], [50+3*i,70-3*i,2000, 1200], 12, true, ~show_figs);
            draw_time = draw_time + toc(draw_tic);
        catch
            fprintf(2,"error drawing  figure %d\n", i);
            err=true;
        end
        if (~err)
            save_tic = tic;
            fprintf("saving figure %s\n",fullfile(figdir,"pngs",strcat(h_nm(i),".png")));
            saveas(h(i),  fullfile(figdir,"pngs",strcat(h_nm(i),".png")));
%             fprintf(" %s\n", fullfile(figdir,"figs",strcat(h_nm(i),".fig")));
%             savefig(h(i), fullfile(figdir,"figs",strcat(h_nm(i),".fig")));
            if (~show_figs), close(i); end
            save_time = save_time + toc(save_tic);
        end
    end
    
    qc = qc([5,6,7,8,9],:,:);
    qc_ttl = qc_ttl([5,6,7,8,9]);
    matname = fullfile(matdir,sprintf("%s_%s_%s_QC.mat",  obs_src,mdl,scenario));

    fprintf("data saved to %s\n", matname);
    save(matname, "obs_src","mdl","scenario", "qc","qc_ttl");
    
        
    elapsed = toc(start_tic);
    
    fprintf("total time:  %6.1f mins\n", elapsed/60);
    fprintf("read  time:  %6.1f mins\n", sum(reads_time)/60);
    fprintf("calc  time:  %6.1f mins\n", sum(calcs_time)/60);
    fprintf("draw  time:  %6.1f mins\n", sum(draw_time)/60);
    fprintf("save  time:  %6.1f mins\n", sum(save_time)/60);
    
    clear qc;
end

