function [elevs,ghcn_dist, ghcn_az, is_ghcn, ghcnID, ghcnName, ghcnLat, ghcnLon, ghcnElev, gsn, hcn, crn, wmoID] = get_ghcn_info(stnID, stnLat, stnLon, maxdist_deg, verbose, varname)
% returns GHCN station info about a site, from an extension of the ghcnd-stations and ghcn-inventory files.
%   tries to match first using stnID against ghcn table's stnID column.
%   if no match, looks for the closest station within maxdist-deg and uses that.
%   if still no match, elev for that station is set to nan.
%
%   Inputs:
%       stnID       single station ID, or matrix of strings of stnIDs
%       stnLat      latitudes for each stnID
%       stnLon      longitudes for each stnID
%       maxdist_deg max distance, in degrees, to use ghcn elevation.  (.001 deg, (~111m,  is reasonable.  This is 3
%                       pixels in the GDEM2 database.)
%       verbose     boolean.  If true, reports progress (can take a while to run...)
%       ghcn_tbl    table read from ghcnd-stations.txt file.  Optional.  Will read from location below if not present.
%
%   Outputs:
%       elevs       column vector of elevations, 1 per station.  elev(i) set to nan if no match.
%       ghcn_dist   distance in km to nearest ghcn station.
%       ghcn_az     az direction to nearest ghcn station.
%       is_ghcn     boolean, true if stnID is in ghcn_tbl (from ghcnd-stations.txt)
%       ghcnID      stnID of nearest GHCN station
%       ghcnName    name of nearest GHCN station
%       ghcnLat     latitude of nearest GHCN station
%       ghcnLon     longitude of nearest GHCN station
%       gsn         gsn, hcn, crn & wmo_id fields for any station where is_ghcn is true.
%       hcn
%       crn
%       wmoID       
%
%-------------------------------

    if (~exist('verbose','var') || isempty(verbose)), verbose = false; end

    basedir = "/Volumes/jcsf_data/ghcn/2019-02-13";
    ghcn_stn_fname = fullfile(basedir, "ghcnd-stations.txt");
    ghcn_stns = ghcn_station_info(ghcn_stn_fname, false, true);
    ghcn_inv_fname = fullfile(basedir,sprintf("ghcnd-inventory_%s.csv", varname));
    ghcn_inv  = readtable(ghcn_inv_fname);

    stnID = string(stnID);     % in case we get cell array of chars.

    gpts = size(ghcn_inv,1);
    npts = length(stnID);
    elevs       = nan(npts,1);
    ghcn_dist   = zeros(npts,1);
    ghcn_az     = zeros(npts,1);
    is_ghcn     = false(npts,1);
    ghcnID      = strings(npts,1);  % stnID of nearest ghcn station
    ghcnName    = strings(npts,1);
    ghcnLat     = zeros(npts,1);
    ghcnLon     = zeros(npts,1);
    ghcnElev    = nan(npts,1);
    gsn         = false(npts,1);
    hcn         = false(npts,1);
    crn         = false(npts,1);
    wmoID       = strings(npts,1);
    nmismatches = 0;
    
    radius = 6371;  % radius of earth in km
    km_per_degree = 2*pi*radius/360;
    
    nmissing = 0;
    missing_stns = strings(0,0);
    nmissing_elevs = 0;
    for i=1:npts
            % make sure it
        jx=find(strcmp(stnID(i), ghcn_inv.stnID), 1);
        ix=find(strcmp(stnID(i), ghcn_stns.stn_id), 1);
        if (~isempty(ix))       % if it's in the inventory for the variable, 
                                % i.e., if we have data for this station for this variable, then use it directly
            if (isnan(ghcn_stns.elev(ix)))
                ghcn_stns.elev(ix) = gdem2_elevation(ghcn_stns.lat(ix), ghcn_stns.lon(ix), false); 
                nmissing_elevs = nmissing_elevs+1;
%                 fprintf("\nmissing ghcn elevation %d retrieved for %d of %d\n", nmissing_elevs, i, npts);
%                 disp(ghcn_stns(ix,:));
                if (verbose && mod(nmissing_elevs,10)==0), fprintf('!'); end
            end
            elevs(i) = ghcn_stns.elev(ix);
            is_ghcn(i)=true;
            ghcn_dist(i) = 0;
            ghcn_az(i)   = 0;
            ghcnID(i)    = ghcn_stns.stn_id(ix);
            ghcnName(i)  = ghcn_stns.stn_name(ix);
            ghcnLat(i)   = ghcn_stns.lat(ix);
            ghcnLon(i)   = ghcn_stns.lon(ix);
            ghcnElev(i)  = ghcn_stns.elev(ix);
            gsn(i)       = ghcn_stns.gsn(ix);
            hcn(i)       = ghcn_stns.hcn(ix);
            crn(i)       = ghcn_stns.crn(ix);
            wmoID(i)     = ghcn_stns.wmo_id(ix);
            if (stnLat(i) ~= ghcn_stns.lat(ix) || stnLon(i) ~= ghcn_stns.lon(ix))
                nmismatches=nmismatches+1;
                if (verbose && nmismatches <= 100)
                    fprintf('\nwarning:  stn %6d:  %s lat/lon (%.5f, %.5f) mismatch with ghcn stn %6d: %s (%s) lat/lon (%.5f, %.5f)\n',...
                            i, stnID(i), stnLat(i),stnLon(i), ix, ghcn.stnID(ix), ghcn.stnName(i), ghcn.lat(ix),ghcn.lon(ix));
                end
            end
            if (isempty(jx))
                if (verbose)
                    fprintf("\nwarning:  stn %6d: %s : found in QC'd data and ghcn station file %s, but not found in inventory file %s for var %s .  '%s'\n", i, ghcnID(i), ghcn_stn_fname, ghcn_inv_fname, varname, ghcnName(i));
                end
                nmissing = nmissing + 1;
                missing_stns(end+1,1) = ghcnID(i); %#ok<AGROW>
            end
        else
                % find the closest station in the inventory table
            stnlat = repmat(stnLat(i), gpts,1);
            stnlon = repmat(stnLon(i), gpts, 1);
            [dists,az] = distance(stnlat, stnlon, ghcn_inv.lat, ghcn_inv.lon);
            jx = find(dists == min(dists),1);
                % get distance and direction from inventory table
            ghcn_dist(i) = dists(jx) * km_per_degree;
            ghcn_az(i) = az(jx);
                % find it in the main ghcn-stations table
            ix=find(strcmp(ghcn_inv.stnID(jx), ghcn_stns.stn_id), 1);
            ghcnID(i)    = ghcn_stns.stn_id(ix);
            ghcnName(i)  = ghcn_stns.stn_name(ix);
            ghcnLat(i)   = ghcn_stns.lat(ix);
            ghcnLon(i)   = ghcn_stns.lon(ix);
            if (isnan(ghcn_stns.elev(ix)))
                ghcn_stns.elev(ix) = gdem2_elevation(ghcn_stns.lat(ix), ghcn_stns.lon(ix), false); 
                nmissing_elevs = nmissing_elevs+1;
%                 fprintf("\nmissing ghcn elevation %d retrieved for %d of %d\n", nmissing_elevs, i, npts);
%                 disp(ghcn_stns(ix,:));
%                 fprintf("for non-ghcn station at lat/lon: %.4f %.4f\n", stnLat(i), stnLon(i));
                if (verbose && mod(nmissing_elevs,10)==0), fprintf('!'); end
            end
            ghcnElev(i)  = ghcn_stns.elev(ix);
                % if we're very close, use GHCN station's elevation.  
                % Otherwise, we'll look for all missing elevations later.
            if (dists(jx) < maxdist_deg)
                elevs(i) = ghcn_stns.elev(ix);
            end
        end        
            % some GHCN stations are missing their elevations, so get elevation if it's still missing.
        if (verbose)   % show some progress
            show_progress(i, npts);
        end
    end
    if (verbose)
        fprintf("%3d GHCN stations missing elevations\n", nmissing_elevs);
        fprintf("%3d stations missing in inventory\n", nmissing);
        fprintf("%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s \n", missing_stns);
        fprintf("\n");
    end    
    ghcn_dist = round(ghcn_dist, 3);    % to nearest meter
    ghcn_az = round(ghcn_az,4);         % 1 arc-second = 2.78e-4.
     
end

