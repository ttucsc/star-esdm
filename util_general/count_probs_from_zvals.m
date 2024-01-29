function count_probs_from_zvals(mdl, scenario, grgn, basedir, use_tllc_time, latrange, lonrange, skip_read, fignum, matname)
% count_probs_from_zvals(mdl, scenario, grgn, basedir, use_tllc_time, latrange, lonrange, skipit)
%
% gets histograms of probability values for locations with tmax<tmin (bad) and where no points have tmax<tmin (good).
% Prints fractional counts for several 'standard' quantiles to see if there are biases in the distributions.
% Breaks counts down by months, by seasons (DJF, MAM, etc.) and total counts.
%
%   Program is designed to verify that ARRM_V2's distributions are in aggregate unbiased, and to see if there is a 
%   difference in the distributions between locations with tmax<tmin vs good locations.  
%
%   In general, tmax<tmin occurs at high elevations;  there is a bimodal distribution of these points.  See program
%   qc_tmin_problem.m

    if (~exist("latrange","var")), latrange = []; end
    if (~exist("lonrange","var")), lonrange = []; end
    if (~exist("fignum","var") || isempty(fignum)), fignum = 10; end
    if (~exist("matname","var")), matname = strings(0); end

    if (isempty(matname) || strlength(matname)==0)
        if (~isempty(latrange))
            latstr=sprintf("_lat%.4f_%.4f", latrange);
        else
            latstr="";
        end
        if (~isempty(lonrange))
            lonstr=sprintf("_lon%.4f_%.4f",lonrange);
        else
            lonstr="";
        end
        matname=sprintf("zval_counts/zval_counts_%s_%s%s%s.mat",mdl, scenario,latstr,lonstr);
    end
    [~,on] = get_hostname();
    if (~isempty(basedir))
        src_dir = basedir;
    elseif (on.hpcc_system)
        if (use_tllc_time)
            src_dir = fullfile("/lustre/scratch/iscottfl/downscaled/cmip6/nclimgrid",scenario,"tllc_time");
        else
            src_dir = fullfile("/lustre/scratch/iscottfl/downscaled/cmip6/nclimgrid",scenario);                
        end
    else
        error("not supported yet");
    end

    if (use_tllc_time)
        fn_tmax = sprintf("downscaled.%s.r1i1p1f1.tasmax.%s.%s.nclimgrid.nclimgrid.1950.2100.tllc_time.nc", mdl, scenario, grgn);
        fn_tmin = sprintf("downscaled.%s.r1i1p1f1.tasmin.%s.%s.nclimgrid.nclimgrid.1950.2100.tllc_time.nc", mdl, scenario, grgn);
    else
        fn_tmax = sprintf("downscaled.%s.r1i1p1f1.tasmax.%s.%s.nclimgrid.nclimgrid.1950.2100.nc", mdl, scenario, grgn);
        fn_tmin = sprintf("downscaled.%s.r1i1p1f1.tasmin.%s.%s.nclimgrid.nclimgrid.1950.2100.nc", mdl, scenario, grgn);
    end

    nc_tmax = ncdf(fullfile(src_dir, fn_tmax));
    nc_tmin = ncdf(fullfile(src_dir, fn_tmin));

    [lats,lons] = ncdf_get_latlons(nc_tmax);
    [~,~,~,~,~,tstamps] = ncdf_get_time_info(nc_tmax, "time",true);
    ndays = length(tstamps);
    nyrs = ndays/365;
    
    if (~isempty(latrange))
        latix = find(lats>latrange(1)-.0001 & lats<latrange(end)+.0001);
    else
        latix = 1:length(lats);
    end
    nlats = length(latix);
    latrng = lats([latix(1),latix(end)]);
    if (~isempty(lonrange))
        if (lonrange(1)<0), lonrange = lonrange+360; end
        lonix = find(lons>lonrange(1)-.0001 & lons<lonrange(end)+.0001);
    else
        lonix = 1:length(lons);
    end
    nlons = length(lonix);
    lonrng = lons([lonix(1),lonix(end)]);

    fprintf("tmax:  %s\n", nc_tmax.Filename);
    fprintf("tmin:  %s\n", nc_tmin.Filename);
    fprintf("lats:  %4d range %9.4f - %9.4f\n", nlats, latrng());
    fprintf("lons:  %4d range %9.4f - %9.4f\n", nlons, lonrng());

    zedges = -6.0:.01:6.0;
    nedges=length(zedges);
    nbins = nedges-1;
    bins  = (zedges(1:end-1)+zedges(2:end))/2;
    
        % monthly counts for z-vals
    z_hmax_good = zeros(13,nbins);  % counts    tmax z-vals on good points by month (13=total)
    z_hmin_good = zeros(13,nbins);  %           tmin z-vals on good points
    z_hmax_bad  = zeros(13,nbins);
    z_hmin_bad  = zeros(13,nbins);
    z_hmax_good_cum = zeros(13,nbins);          % cumulative distributions of z-vals for tmax on good points
    z_hmin_good_cum = zeros(13,nbins);
    z_hmax_bad_cum  = zeros(13,nbins);
    z_hmin_bad_cum  = zeros(13,nbins);
    nz_hmax_good = zeros(13,1);                 % # of tmax good readings
    nz_hmin_good = zeros(13,1);
    nz_hmax_bad  = zeros(13,1);
    nz_hmin_bad  = zeros(13,1);
    
        % season counts for z-vals
    z_hmax_good_ssn = zeros(5,nbins);   % z-vals by season. 
    z_hmin_good_ssn = zeros(5,nbins);
    z_hmax_bad_ssn  = zeros(5,nbins);
    z_hmin_bad_ssn  = zeros(5,nbins);
    z_hmax_good_cum_ssn = zeros(5,nbins);
    z_hmin_good_cum_ssn = zeros(5,nbins);
    z_hmax_bad_cum_ssn  = zeros(5,nbins);
    z_hmin_bad_cum_ssn  = zeros(5,nbins);
    nz_hmax_good_ssn = zeros(5,1);
    nz_hmin_good_ssn = zeros(5,1);
    nz_hmax_bad_ssn  = zeros(5,1);
    nz_hmin_bad_ssn  = zeros(5,1);
    
       
        % days of year for each month 
    mon_days={1:31;32:59;60:90;91:120;121:151;152:181;182:212;213:243;244:273;274:304;305:334;335:365;1:365};
    mon_lbls=["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Year"];
    ssn     =[1,1,2,2,2,3,3,3,4,4,4,1];  % which months go in to which seasons
    ssn_lbls=["DJF","MAM","JJA","SON","Year"];
    
  if (skip_read)   % if data has been run already, we saved the counts in test_counts.mat.
      my_fignum = fignum;
      load(matname); %#ok<LOAD>
      fignum = my_fignum;
  else
% [ ncdf_out, vdata, nc_in ] = ncdf_read_file(fname, varName, latrange, lonrange, daterange, out_calendar, do_random, out_varName, is_indexes, closure_flags)
    for ilon = 1:nlons
        
        careful!  this needs variable zval_name extracted from input files.  could be named packed_zvals (older files, prior to 12/30/21) or zvals (newer, as of 12/30/21)
                % read entire column of latitudes for current longitude
        [~,tmax_vals]  = ncdf_read_file(nc_tmax,"tasmax",      latix,lonix(ilon),[],"365-day",false,"tasmax",      true);
        [~,tmin_vals]  = ncdf_read_file(nc_tmin,"tasmin",      latix,lonix(ilon),[],"365-day",false,"tasmin",      true);
        [~,tmax_zvals] = ncdf_read_file(nc_tmax,"packed_zvals",latix,lonix(ilon),[],"365-day",false,"packed_zvals",true);
        [~,tmin_zvals] = ncdf_read_file(nc_tmin,"packed_zvals",latix,lonix(ilon),[],"365-day",false,"packed_zvals",true);
        
        bad_days = tmax_vals < tmin_vals;
        bad_sites = sum(bad_days,1,"omitnan")>0;
        good_sites = ~bad_sites;
%         tmax_vals  = reshape(tmax_vals,  365,nyrs,inf);
%         tmin_vals  = reshape(tmin_vals,  365,nyrs,inf);
        tmax_zvals = reshape(tmax_zvals, 365,nyrs,nlats);   % reshape by year.
        tmin_zvals = reshape(tmin_zvals, 365,nyrs,nlats);
        
            % for each month, 
       for jmon=1:13
                % separate good and bad sites for current month
            good_tmax = tmax_zvals(mon_days{jmon},:,good_sites);
            good_tmin = tmin_zvals(mon_days{jmon},:,good_sites);
            bad_tmax  = tmax_zvals(mon_days{jmon},:,bad_sites);
            bad_tmin  = tmin_zvals(mon_days{jmon},:,bad_sites);
                % histogram data for all good and all bad sites
            tmhg = histcounts(good_tmax(:),zedges);
            tnhg = histcounts(good_tmin(:),zedges);
            tmhb = histcounts(bad_tmax(:), zedges);
            tnhb = histcounts(bad_tmin(:), zedges);
                    % add to monthly histograms
            z_hmax_good(jmon,:) = z_hmax_good(jmon,:) + tmhg;
            z_hmin_good(jmon,:) = z_hmin_good(jmon,:) + tnhg;
            z_hmax_bad(jmon,:)  = z_hmax_bad(jmon,:)  + tmhb;
            z_hmin_bad(jmon,:)  = z_hmin_bad(jmon,:)  + tnhb;
                    % and add to seasons histograms
            if (jmon ~= 13)
                z_hmax_good_ssn(ssn(jmon),:) = z_hmax_good_ssn(ssn(jmon),:) + tmhg;
                z_hmin_good_ssn(ssn(jmon),:) = z_hmin_good_ssn(ssn(jmon),:) + tnhg;
                z_hmax_bad_ssn(ssn(jmon),:)  = z_hmax_bad_ssn(ssn(jmon),:)  + tmhb;
                z_hmin_bad_ssn(ssn(jmon),:)  = z_hmin_bad_ssn(ssn(jmon),:)  + tnhb;
            else
                z_hmax_good_ssn(5,:) = z_hmax_good_ssn(5,:) + tmhg;
                z_hmin_good_ssn(5,:) = z_hmin_good_ssn(5,:) + tnhg;
                z_hmax_bad_ssn(5,:)  = z_hmax_bad_ssn(5,:)  + tmhb;
                z_hmin_bad_ssn(5,:)  = z_hmin_bad_ssn(5,:)  + tnhb;
            end
        end
        show_progress(ilon, nlons, min(100,nlons));
        if (min(tmax_zvals,[],"all","omitnan")<-6 || max(tmax_zvals,[],"all","omitnan")>6 || min(tmin_zvals,[],"all","omitnan")<-6 || max(tmin_zvals,[],"all","omitnan")>6)
            fprintf("zvals out of range at %.4f\n", lons(lonix));
        end
    end
       clear tmax_vals tmin_vals tmax_zvals tmin_zvals;
    save(matname)
  end
  
    figs=fignum:fignum+3;
    fignum = fignum+4;
    for i=1:4
        h=figure(figs(i));
        clf;
        h.Position=[2600+40*(i-1),250-40*(i-1),2000,1100];
        h.Units="pixels";
    end

%   jpos = [1,2,4,5,3]; % position on plot.
    normvals = normpdf(bins)*.01;
    % now calculate, plot and print out probability distributions.
        % season CDFs
    for j=1:5
        
        z_hmax_good_cum_ssn(j,:) = cumsum(z_hmax_good_ssn(j,:))/sum(z_hmax_good_ssn(j,:));
        z_hmin_good_cum_ssn(j,:) = cumsum(z_hmin_good_ssn(j,:))/sum(z_hmin_good_ssn(j,:));
        z_hmax_bad_cum_ssn(j,:)  = cumsum(z_hmax_bad_ssn(j,:))/sum(z_hmax_bad_ssn(j,:));
        z_hmin_bad_cum_ssn(j,:)  = cumsum(z_hmin_bad_ssn(j,:))/sum(z_hmin_bad_ssn(j,:));
        nz_hmax_good_ssn(j) = sum(z_hmax_good_ssn(j,:));
        nz_hmin_good_ssn(j) = sum(z_hmin_good_ssn(j,:));
        nz_hmax_bad_ssn(j)  = sum(z_hmax_bad_ssn(j,:));
        nz_hmin_bad_ssn(j)  = sum(z_hmin_bad_ssn(j,:));
        figure(figs(1));
        subplot(3,2,j);
        plot(bins, normvals,"c-","linewidth",5);
        hold on;
        plot(bins,z_hmax_good_ssn(j,:)/nz_hmax_good_ssn(j),"b-",bins,z_hmax_bad_ssn(j,:)/nz_hmax_bad_ssn(j),"r-");
        hold off;
        title(sprintf("tasmax by season, %s", ssn_lbls(j)));
        grid on;
        legend("normal","good","bad");
        figure(figs(2));
        subplot(3,2,j);
        plot(bins, normvals,"c-","linewidth",5);
        hold on;
        plot(bins,z_hmin_good_ssn(j,:)/nz_hmin_good_ssn(j),"b-",bins,z_hmin_bad_ssn(j,:)/nz_hmin_bad_ssn(j),"r-");
        hold off;
        title(sprintf("tasmin by season, %s", ssn_lbls(j)));
        grid on;
        legend("normal","good","bad");
        figure(figs(3));
        subplot(3,2,j);
        plot(bins,bins,"c-","linewidth",5);
        hold on;
        z1 = norminv(z_hmax_good_cum_ssn(j,:));
        z2 = norminv(z_hmax_bad_cum_ssn(j,:));
        plot(bins,z1,"b-",bins,z2,"r-");
        hold off;
        title(sprintf("tasmax by season, %s", ssn_lbls(j)));
        grid on;
        legend("normal","good","bad","location","southeast");
        figure(figs(4));
        subplot(3,2,j);
        plot(bins,bins,"c-","linewidth",5);
        hold on;
        z1 = norminv(z_hmin_good_cum_ssn(j,:));
        z2 = norminv(z_hmin_bad_cum_ssn(j,:));
        plot(bins,z1,"b-",bins,z2,"r-");
        hold off;
        title(sprintf("tasmin by season, %s", ssn_lbls(j)));
        grid on;
        legend("normal","good","bad","location","southeast");
    end

    figs=fignum:fignum+3;
%   fignum = fignum+4;
    for i=1:4
        h=figure(figs(i));
        clf;
        h.Position=[2600+40*(i+3),250-40*(i+3),2000,1100];
        h.Units="pixels";
    end

%   jpos = 1:13; % position on plot.
            % Monthly CDFs
    for j=1:13
        z_hmax_good_cum(j,:) = cumsum(z_hmax_good(j,:))/sum(z_hmax_good(j,:));
        z_hmin_good_cum(j,:) = cumsum(z_hmin_good(j,:))/sum(z_hmin_good(j,:));
        z_hmax_bad_cum(j,:)  = cumsum(z_hmax_bad(j,:))/sum(z_hmax_bad(j,:));
        z_hmin_bad_cum(j,:)  = cumsum(z_hmin_bad(j,:))/sum(z_hmin_bad(j,:));
        nz_hmax_good(j) = sum(z_hmax_good(j,:));
        nz_hmin_good(j) = sum(z_hmin_good(j,:));
        nz_hmax_bad(j)  = sum(z_hmax_bad(j,:));
        nz_hmin_bad(j)  = sum(z_hmin_bad(j,:));
        figure(figs(1));
        subplot(4,4,j);
        plot(bins, normvals,"c-","linewidth",5);
        hold on;
        plot(bins,z_hmax_good(j,:)/nz_hmax_good(j),"b-",bins,z_hmax_bad(j,:)/nz_hmax_bad(j),"r-");
        hold off;
        title(sprintf("tasmax by month, %s", mon_lbls(j)));
        grid on;
        legend("normal","good","bad");
        figure(figs(2));
        subplot(4,4,j);
        plot(bins, normvals,"c-","linewidth",5);
        hold on;
        plot(bins,z_hmin_good(j,:)/nz_hmin_good(j),"b-",bins,z_hmin_bad(j,:)/nz_hmin_bad(j),"r-");
        hold off;
        title(sprintf("tasmin by month, %s", mon_lbls(j)));
        grid on;
        legend("normal","good","bad");
        figure(figs(3));
        subplot(4,4,j);
        plot(bins,bins,"c-","linewidth",5);
        hold on;
        z1 = norminv(z_hmax_good_cum(j,:));
        z2 = norminv(z_hmax_bad_cum(j,:));
        plot(bins,z1,"b-",bins,z2,"r-");
        hold off;
        title(sprintf("tasmax by month, %s", mon_lbls(j)));
        grid on;
        legend("normal","good","bad","location","southeast");
        figure(figs(4));
        subplot(4,4,j);
        plot(bins,bins,"c-","linewidth",5);
        hold on;
        z1 = norminv(z_hmin_good_cum(j,:));
        z2 = norminv(z_hmin_bad_cum(j,:));
        plot(bins,z1,"b-",bins,z2,"r-");
        hold off;
        title(sprintf("tasmin by month, %s", mon_lbls(j)));
        grid on;
        legend("normal","good","bad","location","southeast");
    end
    probs = [.0001,.001,.01,normcdf(-2),.05,.1,normcdf(-1),.5,normcdf(1),.9,.95,normcdf(2),.99,.999,.9999,1];
    np = length(probs);
    pzvals = norminv(probs);
    pzvals_ix = zeros(1,np);
    for i=1:(np-1)
        pzvals_ix(i) = find(zedges<=pzvals(i),1,"last");
    end
    pzvals_ix(np) = nbins;

        % seasonal counts
    fprintf("%s %s Seasonal\n", mdl, scenario)
    for j=1:5
        fprintf("         prob: ");
        fprintf("%7.4f  ", probs);
        fprintf("\n");
        print_fracs(z_hmax_good_cum_ssn(j,:), pzvals_ix, nz_hmax_good_ssn(j), "tmax_good", ssn_lbls(j));
        print_fracs(z_hmin_good_cum_ssn(j,:), pzvals_ix, nz_hmin_good_ssn(j), "tmin_good", ssn_lbls(j));
        print_fracs(z_hmax_bad_cum_ssn(j,:),  pzvals_ix, nz_hmax_bad_ssn(j),  "tmax_bad ", ssn_lbls(j));
        print_fracs(z_hmin_bad_cum_ssn(j,:),  pzvals_ix, nz_hmin_bad_ssn(j),  "tmin_bad ", ssn_lbls(j));
        fprintf("\n");
    end
    fprintf("\n");
    
    fprintf("         prob: ");
    fprintf("%7.4f  ", probs);
    fprintf("\n");
    for j=1:5
        print_fracs(z_hmax_good_cum_ssn(j,:), pzvals_ix, nz_hmax_good_ssn(j), "tmax_good", ssn_lbls(j));
    end
    fprintf("\n");
    fprintf("         prob: ");
    fprintf("%7.4f  ", probs);
    fprintf("\n");
    for j=1:5
        print_fracs(z_hmin_good_cum_ssn(j,:), pzvals_ix, nz_hmin_good_ssn(j), "tmin_good", ssn_lbls(j));
    end
    fprintf("\n");
    fprintf("         prob: ");
    fprintf("%7.4f  ", probs);
    fprintf("\n");
    for j=1:5
        print_fracs(z_hmax_bad_cum_ssn(j,:),  pzvals_ix, nz_hmax_bad_ssn(j),  "tmax_bad ", ssn_lbls(j));
    end
    fprintf("\n");
    fprintf("         prob: ");
    fprintf("%7.4f  ", probs);
    fprintf("\n");
    for j=1:5
        print_fracs(z_hmin_bad_cum_ssn(j,:),  pzvals_ix, nz_hmin_bad_ssn(j),  "tmin_bad ", ssn_lbls(j));
    end

            % Monthly counts
    fprintf("\n%s %s Monthly\n", mdl, scenario)
    for j=1:13
        fprintf("         prob: ");
        fprintf("%7.4f  ", probs);
        fprintf("\n");
        print_fracs(z_hmax_good_cum(j,:), pzvals_ix, nz_hmax_good(j), "tmax_good", mon_lbls(j));
        print_fracs(z_hmin_good_cum(j,:), pzvals_ix, nz_hmin_good(j), "tmin_good", mon_lbls(j));
        print_fracs(z_hmax_bad_cum(j,:),  pzvals_ix, nz_hmax_bad(j),  "tmax_bad ", mon_lbls(j));
        print_fracs(z_hmin_bad_cum(j,:),  pzvals_ix, nz_hmin_bad(j),  "tmin_bad ", mon_lbls(j));
        fprintf("\n");
    end
    fprintf("         prob: ");
    fprintf("%7.4f  ", probs);
    fprintf("\n");
    for j=1:13
        print_fracs(z_hmax_good_cum(j,:), pzvals_ix, nz_hmax_good(j), "tmax_good", mon_lbls(j));
    end
    fprintf("\n");
    fprintf("         prob: ");
    fprintf("%7.4f  ", probs);
    fprintf("\n");
    for j=1:13
        print_fracs(z_hmin_good_cum(j,:), pzvals_ix, nz_hmin_good(j), "tmin_good", mon_lbls(j));
    end
    fprintf("\n");
    fprintf("         prob: ");
    fprintf("%7.4f  ", probs);
    fprintf("\n");
    for j=1:13
        print_fracs(z_hmax_bad_cum(j,:),  pzvals_ix, nz_hmax_bad(j),  "tmax_bad ", mon_lbls(j));
    end
    fprintf("\n");
    fprintf("         prob: ");
    fprintf("%7.4f  ", probs);
    fprintf("\n");
    for j=1:13
        print_fracs(z_hmin_bad_cum(j,:),  pzvals_ix, nz_hmin_bad(j),  "tmin_bad ", mon_lbls(j));
    end
end

function print_fracs(zh_cum, pix, nz, lbl, mon_lbl)

    fprintf("%-9s %-4s:", lbl, mon_lbl);
    np = length(pix);
    for i=1:np
        fprintf("%8.5f ", zh_cum(pix(i)));
    end
    fprintf(" count: %12d\n", nz);
end
                







