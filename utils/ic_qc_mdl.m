function ic_qc_mdl(mdl, scenario, grgn, latrng, lonrng, daterange, show_figs, part_figs)
%
%   program to run some basic QC on CMIP6 model data.
%
%   
%
%   Ian Scott-Fleming, 10/2021
%
%--------------------------

    start_tic = tic;
    [hostname, on] = get_hostname();
    
    fprintf("running on host %s, obs_src %s, mdl %s, scenario %s\n", hostname, obs_src, mdl, scenario);
    
    if (on.hpcc_system)
            mdl_src_dir = fullfile("/lustre/research/hayhoe/cmip6_tllc",scenario);
    else
        mdl_src_dir = "/Volumes/lacie_1/data/work";
    end
    
    fn_mdl_max_fut  = sprintf("tasmax.%s.r1i1p1f1.%s.%s.2015.2100.nc", mdl, scenario, grgn);
    fn_mdl_max_hist = sprintf("tasmax.%s.r1i1p1f1.historical.%s.%s.1950.2100.nc", mdl, scenario, grgn);
    fn_mdl_min_fut  = sprintf("tasmin.%s.r1i1p1f1.%s.%s.1950.2014.nc", mdl, scenario, grgn);
    fn_mdl_min_hist = sprintf("tasmin.%s.r1i1p1f1.historical.%s.%s.1950.2100.nc", mdl, scenario, grgn);

    
    if (~exist("daterange","var") || isempty(daterange))
        daterange = [1950,1,1;2000,12,31];
    end
    
    if (~exist("show_figs","var") || isempty(show_figs))
        show_figs = true;
    end
    
    if (~exist("par_figs","var") || isempty(part_figs))     % to display a few figures while processing.
        part_figs = show_figs;
    end
    
    mdl_max_vname = "tasmax";
    mdl_min_vname = "tasmin";

    fprintf("fn_mdl_max:  %s %s\n", fn_mdl_max_hist, fn_mdl_max_fut);
    fprintf("fn_mdl_min:  %s %s\n", fn_mdl_min_hist, fn_mdl_min_fut);    
        
    nc_mdl_max_fut = ncdf(fullfile(mdl_src_dir, scenario, fn_mdl_max_fut));
    nc_mdl_min_fut = ncdf(fullfile(mdl_src_dir, scenario, fn_mdl_min_fut));
    nc_mdl_max_hist = ncdf(fullfile(mdl_src_dir, "historical", fn_mdl_max_hist));
    nc_mdl_min_hist = ncdf(fullfile(mdl_src_dir, "historical", fn_mdl_min_hist));
    
    nc_mdl_max_names = [nc_mdl_max_hist, nc_mdl_max_fut];
    nc_mdl_min_names = [nc_mdl_min_hist, nc_mdl_min_fut];
    
    [mlats,mlons] = ncdf_get_latlons(nc_mdl_max_fut);
    
        % make sure lons are in range [-180,180)
    mlons = mod(mlons+180,360)-180;
    
    if (exist('latrng','var') && ~isempty(latrng))
        mlats = mlats((mlats>= latrng(1)-.0001) & (mlats<= latrng(2)+.0001));
    end
    if (exist('lonrng','var') && ~isempty(lonrng))
        mlons = mlons((mlons>= lonrng(1)-.0001) & (mlons<= lonrng(2)+.0001));
    end
    
        % # of lats & lons
    nmlats = length(mlats);
    nmlons = length(mlons);
    
    R = make_refmat(mlats,mlons, nmlats,nmlons, mlats(22)-mlats(21), mlons(22)-mlons(21));      % for drawing maps.
            
   
    fprintf("latrange:  %9.4f - %9.4f (%d lats)\n", mlats(1),mlats(end), nmlats);
    fprintf("lonrange:  %9.4f - %9.4f (%d lons)\n", mlons(1),mlons(end), nmlons);
    
    if (any(abs(mlats-mlats)>.0001) || any(abs(mlats-mlats)>.0001) || nmlats ~= nmlats || nmlons ~= nmlons)
        error("lat/lon mismatch");
    end
    
    mlatrange = [mlats(1), mlats(end)];
%   mlonrange = [mlons(1), mlons(end)];
    
    nqcs = 24;
    qc = nan(nqcs,nmlats,nmlons);

    
    reads_time = zeros(1,nmlons);
    calcs_time = zeros(1,nmlons); 
    yr_range = [ 10*365,  40*365-1; ...  % 1960-1989
                 70*365, 100*365-1; ...  % 2020-2049
                100*365, 130*365-1; ... % 2050-2079
                130*365, 150*365-1; ... % 2080-2099
                ];
    doyr = [1,59; 182, 243]; % jan/feb & jul/aug
    
    if (nmlons <= 100)
        jmod = ceil(nmlons/20);
        figmod = 5*jmod;
    elseif (nmlons <= 400)
        jmod = ceil(nmlons/100);
        figmod = 5*jmod;
    else
        jmod = 5;
        figmod = 25;
    end
    loop_tic = tic;
    t_end = 51*365; % 51 years, 18615 days
    fid = fopen("ic_qc.log","w");
    qc_ttl = strings(1,nqcs);
    qc_ttl(1)  = sprintf("QC1,  # of days Tmax(mdl) >  60 C 1950-2100 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(2)  = sprintf("QC2,  # of days Tmax(mdl) < -40 C 1950-2100 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(3)  = sprintf("QC3,  # of days Tmin(mdl) >  50 C 1950-2100 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(4)  = sprintf("QC4,  # of days Tmin(mdl) > -60 C 1950-2100 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(9)  = sprintf("QC9,  # of days Tmin(mdl) > Tmax(mdl) 1950-2100 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(10) = sprintf("QC10, # of NA's Tmax(mdl) 1950-2100 (16max) %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(11) = sprintf("QC11, # of NA's Tmin(mdl) 1950-2100 (16max) %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(12) = sprintf("QC12, Avg Tmax increasing && increase < 5C %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(13) = sprintf("QC13, Avg Tmin increasing && increase < 5C %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(14) = sprintf("QC14, Repetition: # of days abs(Tmax(i)-Tmax(i-1)) < 1e-5 1950-2100 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(15) = sprintf("QC15, Repetition: # of days abs(Tmin(i)-Tmin(i-1)) < 1e-5 1950-2100 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(16) = sprintf("QC16, # of years mean(Tmax(jan/feb)) > mean(Tmax(jul/aug)) 1950-2100 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(17) = sprintf("QC17, # of years mean(Tmin(jan/feb)) > mean(Tmin(jul/aug)) 1950-2100 %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    qc_ttl(23)  = sprintf("QC9a, # of days Tmin(mdl) > Tmax(mdl) 1950-2100 (20max) %s, %s, %s %d - %d", obs_src, mdl, scenario, daterange(1,1), daterange(2,1));
    
    
    for j=1:nmlons
        
        reads_tic = tic;
        ncm_max = ncdf_read_file(nc_mdl_max_names, mdl_max_vname, mlatrange, mlons(j), [], "365-day");
        ncm_min = ncdf_read_file(nc_mdl_min_names, mdl_min_vname, mlatrange, mlons(j), [], "365-day");
        
        mdl_vix = find(strcmp(ncm_max.varlist(), mdl_max_vname),1);
        
        reads_time(j) = toc(reads_tic);
        calcs_tic = tic;
        mdl_tmax = ncm_max.Variables(mdl_vix).vdata;
        mdl_tmin = ncm_min.Variables(mdl_vix).vdata;
        
        for i=1:nmlats
            
             if (abs(mlons(j) - (-113.7812)) < .001 && i == 30)
                fprintf("lat,lon = %.5f %.5f\n",mlats(i), mlons(j));
                fprintf("here\n");
             end
        
           mtmax = mdl_tmax(:,i);
            mtmin = mdl_tmin(:,i);
%                 % get the data for just 1950-2000
%             mtmax2k = mdl_tmax(1:t_end,i);
%             mtmin2k = mdl_tmin(1:t_end,i);
%             
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
            qc(9,i,j)  =  sum(mtmin  > mtmax);          % # of days where mdl tmin > mdl tmax
            qc(10,i,j) = min(20,sum(isnan(mtmax)));     % # of tmax NAs.  But using log scale so we can see if > 0.
            qc(11,i,j) = min(20,sum(isnan(mtmin)));     % # of tmin NAs.  This is different from R version.
            % qc12-17 are up above
            qc(23,i,j)  =  min(20,sum(mtmin  > mtmax));        % # of days where mdl tmin > mdl tmax, truncated @ 20
            
            
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
            remain = elapsed/j*(nmlons-j);
            fprintf("%4d of %4d complete  read_time %6.1f mins (%6.1f secs),  calc_time %8.3f secs.  remaining: %6.1f mins\n", done, nmlons, sum(reads_time(1:j))/60, reads_time(j), sum(calcs_time(1:j)), remain/60);
        end
        
        if (j==floor(nmlons/8))
            fprintf("here\n");
            fprintf("hello\n");
        end
        
        if (part_figs && mod(j, figmod)==0)
           kqc=[5,6,23,24];
           for kk=1:4
                k=kqc(kk);
                min_cnt = 0;
                max_cnt = max(qc(k,:,:), [], "all","omitnan");
                display_map(qc(k,:,:), qc_ttl(k), R, 99, [2,2,kk], min_cnt, max_cnt, 2, [], [], [], [100 140 2000 1200], 12, false, ~show_figs);
            end
                
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
    draw_time = 0;
    save_time = 0;
    
    for i = 1:nqcs
        err=false;
        try
            draw_tic = tic;
            h(i) = display_map(qc(i,:,:), qc_ttl(i), R, i, [], min_cnt(i), max_cnt(i), 2, [], [], [], [50+3*i,70-3*i,2000, 1200], 16, false, ~show_figs);
            draw_time = draw_time + toc(draw_tic);
        catch
            fprintf(2,"error drawing  figure %d\n", i);
            err=true;
        end
        if (~err)
            save_tic = tic;
            saveas(h(i),  strcat(h_nm(i),".png"));
            savefig(h(i), strcat(h_nm(i),".fig"));
            save_time = save_time + toc(save_tic);
        end
    end
        
    elapsed = toc(start_tic);
    
    fprintf("total time:  %5.1f\n", elapsed/60);
    fprintf("read  time:  %5.1f\n", sum(reads_time));
    fprintf("calc  time:  %5.1f\n", sum(calcs_time));
    fprintf("draw  time:  %5.1f\n", sum(draw_time));
    fprintf("save  time:  %5.1f\n", sum(save_time));
end

