function compare_ncdata(fname1, fname2, varname, lats, lons)
%
%   Plots data and difference of data between two netcdf files for specified lat/lon locations.
%   See "compare_ncfile_data" for similar plots of difference of station data.
%
%   if lats,lons not specified, does for lat/lons closest to 1st 10 test stations.  
%   if lats given, but lons not specified, takes lats to be # of stations to display from test stations.

    if (~exist('lons','var'))
        tbl=readtable("test_stations_25.csv");
        if (exist('lats','var'))
            nsites=lats;            
        else
            nsites=10;
        end
        nsites = min(nsites, length(tbl.lat));
        lats = tbl.lat(1:nsites);
        lons = tbl.lon(1:nsites);
        sitenames = string(tbl.stnName(1:nsites));
    else
        nsites = length(lons);
        sitenames = strings(nsites,1);
        for i=1:length(lons)
            sitenames(i) = sprintf("site %d",i);
        end
    end

%     nc1 = ncdf(fname1,'do_create',false);
%     nc1.loadvars([],true);
%     nc2 = ncdf(fname2,'do_create',false);
%     nc2.loadvars([],true);
%     lats1 = nc1.getvardata('lat');
%     lons1 = nc1.getvardata('lon');
%     lats2 = nc2.getvardata('lat');
%     lons2 = nc2.getvardata('lon');
%     
    for i=1:nsites
%         [latix1, lonix1, latout1, lonout1] = closets_lat_lon_range(lats(i), lons(i), lats1, lons1, 1, 1);
%         [latix2, lonix2, latout2, lonout2] = closets_lat_lon_range(lats(i), lons(i), lats2, lons2, 1, 1);
        
        if (contains(fname1,"gfdl"))       % old hires files, which begin w/ GFDL, were of on longitude by 1/2 gridcell
            nc1 = ncdf_read_file(fname1, varname, lats(i), lons(i)+.15625,[],"365-day");
        else
            nc1 = ncdf_read_file(fname1, varname, lats(i), lons(i),[],"365-day");
        end
        if (contains(fname2,"gfdl")) 
            nc2 = ncdf_read_file(fname2, varname, lats(i), lons(i)+.15625,[],"365_day");
        else
            nc2 = ncdf_read_file(fname2, varname, lats(i), lons(i),[],"365_day");
        end
        
        d1 = nc1.getvardata(varname);
        d1=d1(:,1,1);
        d2 = nc2.getvardata(varname);
        d2=d2(:,1,1);
        
        difs = d2 - d1;
        nonz = sum(difs ~= 0);
        
        x1 = (1:length(d1))/365;
        x2 = (1:length(d2))/365;
        
        figure(10+i);
        if (nonz > 0)
            subplot(2,1,1);
        end
        plot(x1,d1,x2,d2);
        title(sprintf("%s (%.4f %.4f)", sitenames(i), lats(1), lons(i)));
        
        if (nonz > 0)
            subplot(2,1,2);
            plot(x1, difs);
            title(sprintf("difs, mean:  %.5g RMS:  %.5g, max %.5g  min %.5g", mean(difs), std(difs), max(difs), min(difs)));
        end
    end
end
