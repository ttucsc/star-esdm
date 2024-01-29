function QC_add_info_to_latlon_file(latlon_fname, latlon_outname, varname, inventory_name)
%
%   Adds elevation columns to all_stations, using GHCN stnID's elevation when available, and gdem2_elevation(...) when
%   not.
%
%   Should modify this to add timezone as well.
%
%   Really want a version of this to add metadata columns to station netcdf files, rather than to the lat/lon file...
%   This is really just an attempt to add the missing GHCN station-data  info to Ranjini's Lat/Lon station file.
%
%   latlon_fname columns:
%       stnID, stnName, lat, lon, state       Ranjini's lat/lon file  all_ncamerica_latlon.csv
%
%   latlon_outname:  something like:   all_ncamerica_latlon_elev.csv
%
%   varname:  Tmax, Tmin or Prec
%
%   inventory_name:  "ghcnd-inventory_TMAX.csv", etc.

        % be sure to update ghcn_tbl location for future uses!!!

    latlon_fname = string(latlon_fname);
    latlon_outname = string(latlon_outname); 
    varname = string(varname);
    inventory_name = string(inventory_name);    % ghcnd-inventory_TMAX.csv, etc.
    
    start_tic = tic();
    
    maxdist_deg = .001  ;   % if within .001 degrees of a ghcn station, use GHCN elevation..  This is ~111 m, between 3 & 4 pixels in GDEM2 elevation database.
                            % when stn doesn't match a GHCN stnID, use nearby GHCN stnID's lat/lon if within 2 pixels.
                            % this will catch case where location has been rounded to 3 decimals (.001 degrees), ~55 m.
    
    tbl = readtable(latlon_fname);
    
%   tbl = tbl(1:50:end,:);             % for debugging!
    
    tbl.stnID = string(tbl.stnID);
    tbl.stnName = string(tbl.stnName);
    tbl.state = string(tbl.state);
    nr = size(tbl,1);
    tbl.elev = nan(nr,1);
    
        % flag the non-ghcn stations
        
    tbl.is_ghcn = true(nr,1);
    ghcn_tbl = ghcn_station_info("/Volumes/jcsf_data/ghcn/2019-02-13/ghcnd-stations.txt",false, true);    
    [~,nonghcn] = setdiff(tbl.stnID, ghcn_tbl.stn_id);    
    tbl.is_ghcn(nonghcn) = false;
        
        % now remove any  stations that don't have data for this variable.

    varbase = '/Volumes/jcsf_data/data/obs/QC_csv';     % for icsf-jmac.
    % varbase = '/data/obs/organized_raw_stations/all_stations/ncamerica/qc_final';     % for killarney
    dirname = fullfile(varbase,varname);
    vardir = dir(dirname);
    fnames=string({vardir(:).name})';
    [~, ix] = setdiff(tbl.stnID, fnames);
        
    keepers = true(nr,1);
    keepers(ix) = false;
        
        % remove any ghcn stations we don't have data for for this variable.
    tbl = tbl(keepers,:);
    
        % compare directory with inventory, see if we're missing any GHCN stations.
        
    invbase = '/Volumes/jcsf_data/ghcn/2019-02-13';
    inv_tbl = readtable(fullfile(invbase,inventory_name));
    ghcn_ncamerica_stns = region_stations(inv_tbl, 'ncamerica');
    invmissing = setdiff(ghcn_ncamerica_stns, fnames);
    ninvmissing = length(invmissing);
    fprintf("%d stations in GHCN inventory (%s) which are missing QC output\n", ninvmissing, inventory_name);
    if (ninvmissing < 500)
        fprintf("missing stations: \n");
        fprintf("%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s \n", invmissing);
        fprintf("\n");
    end
    
    
%     tbl = tbl(1:100:end,:);    % for debugging;
    
%   [elev, ghcn_dist,ghcn_az, is_ghcn, ghcnName, ghcnLat, ghcnLon, gsn, hcn, crn, wmoID] = get_ghcn_info(tbl.stnID, tbl.lat, tbl.lon, maxdist_deg, true);   
%   [elev, ghcn_dist,ghcn_az, ~,                ghcnName, ghcnLat, ghcnLon, gsn, hcn, crn, wmoID] = get_ghcn_info(tbl.stnID, tbl.lat, tbl.lon, maxdist_deg, true);   
    [elev, ghcn_dist, ghcn_az, ~,       ghcnID, ghcnName, ghcnLat, ghcnLon, ghcnElev, gsn, hcn, crn, wmoID] = get_ghcn_info(tbl.stnID, tbl.lat, tbl.lon, maxdist_deg, true, varname) ;   
    tbl.elev = elev;
    tbl.ghcnID = ghcnID;
    tbl.ghcn_dist = ghcn_dist;
    tbl.ghcn_az = ghcn_az;
    tbl.ghcn_elev = ghcnElev;
%   tbl.is_ghcn = is_ghcn;
    tbl.gsn = gsn;
    tbl.hcn = hcn;
    tbl.crn = crn;
    tbl.wmoID = wmoID;    
    
    
    missing = isnan(tbl.elev);
    fprintf("\n%d missing elevations.  Using gdem2_elevations to replace missing elevations\n", sum(missing));
    tbl.elev(missing) = gdem2_elevation(tbl.lat(missing),tbl.lon(missing),true);

    

    ix = find(tbl.is_ghcn);
    nghcn = length(ix);
    fprintf('ghcn matches:  %d\n', nghcn);
%     ghcnName(ix(1))="hello world";            % for debugging, to force some mismatches.
%     ghcnLat(ix(2))=-99;
%     ghcnLon(ix(3))= -999;
    
        % for ghcn stations, report if lat, lon or name is different.
    compare_ghcn_info(tbl, ghcnName, ghcnLat, ghcnLon, 100);

    writetable(tbl, latlon_outname, 'QuoteStrings',true);
    
    elapsed = toc(start_tic);
    mins = floor(elapsed/60);
    secs = elapsed - 60*mins;
    fprintf('runtime:  %d:%.1f elapsed;  results written to %s\n', mins, secs, latlon_outname);

end

function compare_ghcn_info(tbl, ghcnName, ghcnLat, ghcnLon, maxdiffs)

    npts = size(tbl,1);
    
    name_mismatches = tbl.is_ghcn & ~(strcmp(tbl.stnName, ghcnName));
    lat_mismatches = tbl.is_ghcn & (tbl.lat ~= ghcnLat);
    lon_mismatches = tbl.is_ghcn & (tbl.lon ~= ghcnLon);
    fprintf("ghcn comparison:  %d name mismatches, %d lat mismatches, %d lon mismatches\n", sum(name_mismatches), sum(lat_mismatches), sum(lon_mismatches));
    
    mismatches = name_mismatches | lat_mismatches | lon_mismatches;
    
    nmismatches = sum(mismatches);
    if (nmismatches == 0 || nmismatches > maxdiffs), return; end

    nameflags = repmat(' ',npts,1);
    latflags  = repmat(' ', npts,1);
    lonflags  = repmat(' ', npts,1);
    nameflags(name_mismatches) = '*';
    latflags(lat_mismatches) = '*';
    lonflags(lon_mismatches) = '*';
    
    radius = 6371;  % radius of earth in km
	km_per_degree = 2*pi*radius/360;

fprintf('\n index  stnID          latitude    ghcn_lat longitude    ghcn_lon dist(km) stnName                          ghcnName\n\n');
    for i=1:npts
        if (mismatches(i))
            my_dist = distance(tbl.lat(i), tbl.lon(i), ghcnLat(i), ghcnLon(i)) * km_per_degree;
            fprintf('%6d %-12s %9.4f %c %9.4f %9.4f %c %9.4f %8.3f %-30s %c %-30s\n', i, tbl.stnID(i), ...
                     tbl.lat(i), latflags(i), ghcnLat(i), tbl.lon(i), lonflags(i), ghcnLon(i), my_dist, ...
                     tbl.stnName(i), nameflags(i), ghcnName(i));
        end
    end
end     

function region_stnIDs = region_stations(inv_tbl, region)

    cc_codes = GHCN_country_codes(region);
    keepers = false(size(inv_tbl,1),1);
    for i=1:length(cc_codes)
        code = cc_codes(i);
        if (strlength(code) ~= 2), error('error:  bad country code: %s', code); end
        keepers = keepers | strncmpi(inv_tbl.stnID, code,2);
    end
    region_stnIDs = inv_tbl.stnID(keepers);
end
    
    