function [elevs, avg_elevs, diff_elevs] = conus_elevations(fignum, elevs_in, avg_elevs_in, dmax)


    nc=ncdf("/Volumes/lacie_1/data/downscaled/arrm_v2/conus/downscaled.GFDL-HIRAM-C360_cm3.r1i1p1.tasmax.rcp85.X1X3X4.2086.2115.nc");
    nc.loadvars([],true);
    lats=nc.getvardata("latitude");
    lons=nc.getvardata("longitude");
    [grlons,grlats] = meshgrid(lons,lats);
    grlons=grlons(:);
    grlats=grlats(:);
    
    nclow = ncdf("/Volumes/lacie_1/projects/ian/work/hist/tasmax_day_GFDL-HIRAM-C360_amip_19790101-20381231.H2H3.lowres.llt.nc");
    nclow.loadvars([], true);
    latslow = nclow.getvardata("lat");
    lonslow = nclow.getvardata("lon");
    
    if (~exist('elevs_in','var') || isempty(elevs_in))
        elevs=gdem2_elevation(grlats,grlons, true);
        elevs=reshape(elevs,120,208);
    elseif (ischar_s(elevs_in))
        e = load(elevs_in,"elevs");
        elevs = e.elevs;
    else
        elevs = elevs_in;
    end
    
    elevs(isnan(elevs))=0;

        % get 1-degree average elevations
        
    if (~exist('avg_elevs_in','var') || isempty(avg_elevs_in))

        avg_elevs = nan(size(elevs));

        for x=1:208
            lonlow1 = find(lonslow <= lons(x), 1, 'last');   
            lonlow2 = find(lonslow > lons(x), 1);

            if (isempty(lonlow1)), lonlow1 = lons(x); end
            if (isempty(lonlow2)), lonlow2 = lons(x); end
            mylons = lons >= lonslow(lonlow1) & lons < lonslow(lonlow2);
            for y=1:120
                latlow1 = find(latslow <= lats(y),1,'last');
                latlow2 = find(latslow > lats(y),1);
                if (isempty(latlow1)), latlow1 = lats(y); end
                if (isempty(latlow2)), latlow2 = lats(y); end
                mylats = lats >= latslow(latlow1) & lats < latslow(latlow2);

                gridelevs = elevs(mylats,mylons);

                avg_elevs(mylats, mylons) = nanmean(gridelevs(:));
            end
        end
    elseif (ischar_s(avg_elevs_in))
        e = load(avg_elevs_in);
        avg_elevs = e.avg_elevs;
    else
        avg_elevs = avg_elevs_in;
    end
    
    avg_elevs(isnan(avg_elevs))=0;
    
    diff_elevs = elevs - avg_elevs;

    minhi = min(elevs(:));
    maxhi = max(elevs(:));
    minlow = min(avg_elevs(:));
    maxlow = max(avg_elevs(:));
    mindiff = min(diff_elevs(:));
    maxdiff = max(diff_elevs(:));
    
    if (~exist('dmax','var') || isempty(dmax))
        dmax = max(abs(diff_elevs(:)));
    end
    
%   save("conus_elevs.mat", "elevs", "avg_elevs", "diff_elevs", "lats","lons", "latslow", "lonslow");
    
    
    if (~exist('fignum','var') || isempty(fignum)), return; end
    
    elevs11 = elevs(1,1);
    avg_elevs11 = avg_elevs(1,1);
    avg_elevsee = avg_elevs(end,end);
    
    elevs(1,1) = -max(elevs(:));
    avg_elevs(1,1) = -max(elevs(:));
    avg_elevs(end,end) = max(elevs(:));
    
    diff_elevs(1,1) = -dmax;
    diff_elevs(end,end) = dmax;
    diff_elevs(diff_elevs > dmax) = dmax;
    diff_elevs(diff_elevs < -dmax) = dmax;
    
    R=georefcells([20,50],[-125,-60],[120,208]);

    figure(fignum);
    subplot(3,1,1);
    worldmap(elevs, R);
    meshm(elevs,R);
    ncamerica_and_state_boundaries([],true);
    lowgrid(latslow, lonslow, [20,50],[-125,-60]);
    colorbar("location","eastoutside");
    title(sprintf("Elevations (m)  range: %.0f-%.0f", minhi, maxhi));
    hold off;
    

    subplot(3,1,2);
    worldmap(elevs, R);
    meshm(avg_elevs,R);
    ncamerica_and_state_boundaries([],true);
    lowgrid(latslow, lonslow, [20,50],[-125,-60]);
    colorbar("location","eastoutside");
    title(sprintf("Mean Elevation for gridcells (m)  range:  %.0f-%.0f", minlow, maxlow));
    hold off;
    
    subplot(3,1,3);
    worldmap(elevs, R);
    meshm(diff_elevs,R);
    ncamerica_and_state_boundaries([],true);
    lowgrid(latslow, lonslow, [20,50],[-125,-60]);
    colorbar("location","eastoutside");
    title(sprintf("Elevation difference (m)  range:  %.0f-%.0f", mindiff, maxdiff));
    hold off;
    
    colormap(jet(16));
    
    elevs(1,1) = elevs11;
    avg_elevs(1,1) = avg_elevs11;
    avg_elevs(end,end) = avg_elevsee;
    
end

function lowgrid(latslow, lonslow, latrange, lonrange)

    lonslow = mod(lonslow+180, 360) - 180;
    keepers = latslow >= latrange(1) & latslow <= latrange(2);
    latslow = latslow(keepers);
    keepers = lonslow >= lonrange(1) & lonslow <= lonrange(2);
    lonslow = lonslow(keepers);
    
    for i=1:length(latslow)
        plotm(latslow(i)*ones(size(lonslow)), lonslow, '-','color',[.75,.75,.75]);
    end
    for i=1:length(lonslow)
        plotm(latslow, lonslow(i)*ones(size(latslow)), '-','color',[.75,.75,.75]);
    end
end