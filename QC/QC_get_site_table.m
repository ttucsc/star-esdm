function [siteTbl, nsites, varNames, npts, allTbl] = QC_get_site_table(ncName, varargin)
%function [siteTbl, nsites, varName, npts] = QC_get_site_table(ncName, [key/value pairs])
%
% THIS NEEDS: some sort of date/timestamp and/or hash for verifying table is uptodate if ncName is a QC_stntbl or .mat file
%
%   originally:  function [siteTbl, nsites, varName, npts] = QC_get_site_table(ncName, dateRange, full_or_minyrs, boundRange, show_dates, csvname, loadTblVars)
%returns site table info, with all sites in ncName and their data.
%
%   Inputs:
%       ncName          name of station netcdf file, or .mat file w/ contents of QC_stntbl, or  QC_stntbl from a
%                           previous call to QC_get_site_table
%     optional inputs:%                           
%       start_date      year, or [year,month,day[,hr,min,sec]] of starting date                                 [empty] 
%       end_date        year, or [year,month,day[,hr,min,sec]] of ending date                                   [empty]
%                           If specified, include only stations with data in date range specified
%                           (see also fullYears and minYears option below)
%                           If start_date empty or missing, use file's starting date as start_date
%                           If end_date empty or missing, use file's ending date as end_date.
%
%     optional keyword/value pairs
%       "outcal"            calendar type for output:  "day-365","julian", etc.  see function calendar_length() for list 
%                               of calendar strings accepted.
%       "removeLeaps",t/f   (obsolete;  to support older calls.)default:  [].  true:  sets outcal to "365-day".  false: sets outcal to input file's calendar.    
%       "do_random",t/f     if true and calendar_length(outcal) > calendar_length(calendar), then inserts extra days at
%                               random locations, so inserted nans don't fall on the same day of the year.
%       "fullYears", t/f    if true, include only stations that cover entire date range specified.              [false]
%                           if false, include any station with data anywhere in the date range specified.
%       "minYears", #       if > 0, include only stations with at least minYears of data                        [ 0 ] 
%                                   in station's start/end date range, or in date range, if specified
%       "minPct", #         if > 0, include only stations with at least minPct percentage of valid data         [ 0 ]
%                                   in station's start/end date range, or in date range, if specified
%                                   value is percentage (0.0 - 100.0), not 0-1.  values < 1 will be multiplied by 100.
%       "showDates", t/f    if true, dates returned as datestr(...,'yyyy-mm-dd');  else dates returned as datenums [false] 
%       loadTblVars         if true, loads all netcdf variables (except actual data) which use station as a dimension.  
%                                   i.e., all variables associated with station, except for data.               [false]
%                                   if false, only loads the standard QC station file columns.
%       "stnID",stnids      array of stn identifiers.  stnids can be string array or cell array of chars with stnIDs [empty] 
%                                   or stnNames (or partial names), or array of indexes from original netcdf file.
%                                   if missing, loads info for all stations in file.  (which can be subsetted with a
%                                   future call to QC_get_site_table)
%                                   Can also be an array of [lat,lon] locations to specify locations to keep. (see "closest" for search_type below) 
%       "searchType", typ   "stnID","stnName","any"                                                             [stnID]
%                           or "like_stnID", "like_stnName", "like_any"
%                               add "like" to match in searchType on any partial match.
%                               otherwise, match is exact match only.
%                           or "closest" for when lats and lons are specified.  If [lat,lon] points given and not using closest, then will only keep exact matches.   
%       "tryMatfile", 0,1 or 2, if >0, looks for mat file in same location as netcdf file and loads it if present. [false] 
%                               if == 1, will use .nc file     if timestamp info doesn't match original .nc file
%                               if == 2, will use matfile even if timestamp info doesn't match original .nc file
%       "latrange", [lat1,lat2] keep only stations in half-open latitude  range [lat1,lat2);
%       "lonrange", [lon1,lon2] keep only stations in half-open longitude range [lon1,lon2);
%       "loadData", t/f     If true, loads data for climate variable(s) as well as table information
%                               caution:  can use a LOT of memory unless you limit the latrange & lonrange or provide a
%                               list if stnIDs to include.
%       "maxmem",nbytes     Set maximum memory to use.  If loading data will use more than this, QC_get_site_table will
%                               abort. default:  1 GByte.
%
%   Outputs
%       siteTbl         site table, one entry for each stnID
%       nsites          # of sites in table
%       varName         name of data variable found in netcdf file
%       npts            number of datapoints (one for each time step) in file. 
%                           
%---------------------QC Table output----------------
%   QC table (output) columns:
%           stnID           station ID
%           lat, lon        location  (decimal degrees, usually to 4 decimals;  some only to 2 or 3)
%           stnName         station Name
%           elev            elevation, in meters
%           startDate       start timestamp of first valid record for site, as Matlab datenum (for data's calendar type)
%           endDate         end   timestamp of last  valid record for site, as Matlab datenum (for data's calendar type)
%           pctValid        percent of valid points between startDate and endDate
%           index           index of record in original input file
%       possible additional fields
%           time_zone       (currently set to NA)
%           is_ghcn         boolean.  true if site is a GHCN station
%           nearest_ghcn    stnID of nearest GHCN station, if is_ghcn is false
%           ghcn_elev       elevation of nearest GHCN stations
%           ghcn_dist       distance to nearest GHCN station
%           ghcn_az         direction (degrees AZ from true North) to nearest GHCN station.
%
%   QC table Properties.UserData fields & example values:
%       Properties.UserData contains useful table-wide data extracted from netcdf file's global attributes.
%
%        varName: 'TMAX'                                variable name
%       longName: 'daily max temp'                      variable long name
%          units: 'degC'                                variable units
%         nsites: 70                                    # of sites being returned
%           npts: 24106                                 # of data points (dates) in file.  All stations have the same
%                                                       npts;  data before first recorded date and after last recorded
%                                                       date for any station are set to NAs in the netcdf file.
%         NAFlag: 1.0000e+20                            file's original NA value ('missing or invalid data' flag)
%         ncName: "/Volumes/jcsf_data/data/obs/stations_netcdf/stations.Tmax.ncamerica.1850.2018.20_yrs_min.nc"
%                                                       original input file, full pathname
%  fileDateRange: "1850-01-01 2018-12-31"               date range of file (or of extracted data for some QC programs)
%      dateRange: "1940-01-01 2015-12-31"               date range of data read from file 
%                                                           (== fileDateRange if start_date or end_date not specified)
%       calendar: 'standard'                            calendar type.  365-day, 360-day, standard, etc.
%      timeunits: 'days since 1850-01-01 00:00:00 UTC'  netcdf time units
%          isUTC: 1                                     boolean. flags whether timestamps are UTC or local time.
%          dates: [24106,1 double]                      datenums_cal for each data value (i.e., matlab datenums based on calendar).
%           day1: 675700                                matlab datenum of 'days since' day of timeunits)
%  modifiedDates: false                                 boolean.  flags whether date information has been modified from
%                                                           original source.  I.e., table has subset of original dates
%                                                           true:   date info is subset of date info in original nc file 
%                                                           false:  date info is identical to date info in orig. nc file
%     isQCstntbl: 1                                     boolean.  true if table is a QC station table.
%       allSites: false                                 boolean.  flags whether all stations from netcdf file included
%                                                           true:   all stations included;  
%                                                           false:  only some stations from original file included.
%      fullYears: 0                                     if true, only keeping stations that span entire date range
%                                                           specified by start_date and end_date.
%                                                           if false, keeping stations with data at any point in
%                                                           start/end date range
%       minYears: 0                                     if>0 only keeping stations with at least this # of years of data
%     nctblNames: ["stnID"    "lat"    "lon"    "stnName"    "elevation"    "start_date"    "end_date"    "pct_valid"    "time_zone"    "is_ghcn"    "nearest_ghcn"    "ghcn_distance"    "ghcn_az"    "ghcn_elev"]
%                                                       netCDF file's variable names corresponding to the table column
%                                                           names.  (some are slightly different, e.g. 
%                                                               startDate in table, vs start_date in the file.  
%   table_source: [1?1 struct]                          struct with directory info about the underlying original netcdf
%                                                           file  (name, folder, date, size, file timestamp)
%
%       Percent Valid is percent valid over dates between each station's startDate and endDate as reported.
%       CAUTION:  if loadData is false, percentValid is based on all data in original file, 
%                 even if start_date and end_date were specified.
%
%   2022-09-19      added search by lat/lon, and keeping closest point withing some tolerance.  This lets lats/lons not
%                       be exact.
%
%
%-----------------------------------
    [ncName, start_date, end_date, fullYears, minYears, minPct, allSites, showDates, ...
      loadTblVars,stnID, searchType, loadData, maxmem, tryMatfile, latrange, lonrange, latlon_keepers, varNames, use_varNames, outcal, do_random] = initParams(ncName, nargout, varargin{:});
        % make sure filename has path if in local folder, so we can put full path info into UserData
        
        % NOTE:  start_date, end_date (as specified by user) are datevecs, not datenums.
        
    if (isQCstntbl(ncName) || isNCtbl(ncName))
        siteTbl = ncName;
        nc = [];
    else
        [siteTbl, nc] = read_netcdf_QC_table(ncName, tryMatfile);        % will also read .mat file if ncName is matfile. 
    end
    
%   incal = siteTbl.Properties.UserData.calendar;
    
    if (loadTblVars)
         [siteTbl, nc] = load_table_vars(siteTbl, nc, loadTblVars);
%       siteTbl = load_table_vars(siteTbl, nc, loadTblVars, siteTbl.Properties.UserData.nctblNames, siteTbl.Properties.UserData.varName);
    end       
    
    allTbl = siteTbl;
    
    if (isempty_s(varNames))
        varNames  = siteTbl.Properties.UserData.varName;
    elseif (length(varNames)==1 && strcmpi(varNames,"all"))
%       vnames=nc.varlist();
        keepers = is_climate_variable(nc.varlist());
        varNames = nc.varlist(keepers);
        
%         varNames = strings(0);
%         for i=1:length(vnames)
%             is_clim = is_climate_variables(vnames(i));
%             if (is_clim)
%                 varNames(end+1) = vnames(i); %#ok<AGROW>
%             end
%         end
    end
    calendar = siteTbl.Properties.UserData.calendar;
    dates    = siteTbl.Properties.UserData.dates;
    
        % we need these as datenums, not datevecs, but needed to know the calendar before we could convert them.
    dates_modified = false;
    if (~isempty(start_date))
        start_dnum = datenum_cal(start_date,calendar);
        if (start_dnum ~= dates(1))
            dates_modified = true;
        end
    else
        start_dnum = dates(1);
    end
    if (~isempty(end_date))
        end_dnum = datenum_cal(  end_date,calendar); 
        if (end_dnum ~= dates(end))
            dates_modified = true;
        end
    else
        end_dnum = dates(end);
    end
    
        % in case start, end dates are in text strings instead of datenums.
    if (ischars(siteTbl.startDate)), siteTbl.startDate = datenum_cal(datevec(siteTbl.startDate), calendar); end
    if (ischars(siteTbl.endDate)),   siteTbl.endDate   = datenum_cal(datevec(siteTbl.endDate),   calendar); end
    
            % if stnID is not empty, limit list of stations to those in the list of stnIDs
            %       stnID could be:  indexes to keep, stnIDs or stnNames.
    if (ischars(stnID)), stnID = string(stnID); end
    
    if (~isempty(stnID) && ~(length(stnID)==1 && strcmp(stnID,"all")))
        nsites = size(siteTbl,1);
        keepers = false(nsites,1);
        
        if (isnumeric(stnID))
            for i=1:length(stnID)
                ix = find(siteTbl.index == stnID(i));
                for j=1:length(ix)
                    keepers(ix(j)) = true;
                end
            end
        else
            for i=1:length(stnID)        
                ix = find_site(stnID(i), siteTbl, string(searchType));
                for j=1:length(ix)
                    keepers(ix(j)) = true;
                end
            end
        end
        
        if (sum(keepers)==0), error("error:  no sites matching stnID criteria"); end 
        siteTbl = siteTbl(keepers,:);
    end
    
                % limit stations to specified latlons or to lat/lon region if lats or lons specified.
                
    if (~isempty(latlon_keepers))
        siteTbl = siteTbl(latlon_keepers, :);
    else

        if (~isempty(latrange))
            keepers = siteTbl.lat >=latrange(1) & siteTbl.lat <= latrange(2);
            if (sum(keepers)==0), error("error:  no sites in lat range %.4f - %.4f)", latrange); end 
            siteTbl = siteTbl(keepers,:);
        end        

        if (~isempty(lonrange))
            if (lonrange(1) > 0)
                lons = longitude_360(siteTbl.lon);
            else
                lons = longitude_180(siteTbl.lon);
            end
            keepers = lons >=lonrange(1) & lons <= lonrange(2);
            if (sum(keepers)==0), error("error:  no sites in lon range %.4f - %.4f)", lonrange); end 
            siteTbl = siteTbl(keepers,:);
        end        
    end        
            % if start, end date given, limit list of stations by date 

%     if (isempty(start_date)), my_start_date = dates(1);   else, my_start_date = start_date; end
%     if (isempty(  end_date)),   my_end_date = dates(end); else,   my_end_date = end_date;   end
%     siteTbl.Properties.UserData.dateRangeLimit = sprintf("%s %s", string(datestr(datevec_cal([my_start_date; my_end_date],calendar),'yyyy-mm-dd')));
%     npts = my_end_date - my_start_date + 1;

    if ((start_dnum > dates(1) || end_dnum < dates(end)) || fullYears)
    
        if (~fullYears)
            keepers = (siteTbl.startDate <= end_dnum & siteTbl.endDate >= start_dnum);
        else
            keepers = (siteTbl.startDate <= start_dnum & siteTbl.endDate >= end_dnum);
        end
        if (sum(keepers) == 0)
            error('no valid data in date range %s to %s', datestr_cal(start_dnum,calendar,'yyyy-mm-dd'),datestr_cal(end_dnum,calendar,'yyyy-mm-dd')); 
        end

        siteTbl = siteTbl(keepers,:);
    end
        % if minPct given, limit list of stations to those with at least minPct valid data.
    
    if (minPct > 0)
        keepers = siteTbl.pctValid >= minPct;
        if (sum(keepers) == 0)
            error('no stations left after applying minPct > %.3f %% criteria', minPct); 
        end    
        siteTbl = siteTbl(keepers,:);
    end
    
    
        % if minYears given, limit list of stations to those with at least minYears with valid data.
        
    if (minYears > 0)
        bounded_startDate = max(siteTbl.startDate, start_dnum);
        bounded_endDate   = min(siteTbl.endDate,   end_dnum);

            % years rounded so leap-year issues don't bite us.  If divide by 365.25, but we're short 1 leap-day,
            % then test of >= nyears might fail incorrectly otherwise.

        stn_years = round((bounded_endDate - bounded_startDate+1)/calendar_length(calendar),2);

        keepers = stn_years >= minYears;       % exclude sites with less than minYears data
        if (sum(keepers) == 0)
            error('no stations left after applying minYears > %.1f criteria', minYears); 
        end    
        siteTbl = siteTbl(keepers,:);
    end

   
    nsites = size(siteTbl,1);
        
%     npts = end_date - start_date + 1;    
%     siteTbl.Properties.UserData.npts = npts;
    siteTbl.Properties.UserData.allSites = allSites;
    siteTbl.Properties.UserData.fullYears = fullYears;
    siteTbl.Properties.UserData.minYears  = minYears;
    siteTbl.Properties.UserData.nsites = nsites;
    siteTbl.Properties.UserData.modifiedDates = dates_modified;
    siteTbl.Properties.UserData.dateRange = sprintf("%s %s", string(datestr(datevec_cal([start_dnum; end_dnum],calendar),'yyyy-mm-dd')));
    
    if (showDates)
        try
            siteTbl.startDate = string(datestr_cal(siteTbl.startDate,calendar, "yyyy-mm-dd"));
            siteTbl.endDate   = string(datestr_cal(siteTbl.endDate,  calendar, "yyyy-mm-dd"));
        catch
        end
    end
    
    if (~isempty(stnID) && ~contains(searchType,"like"))
        siteTbl = reorder(siteTbl, stnID, searchType);
    end
    
    if (loadData)
%  This needs to load all data variables.  input "variable" s/b vector of strings, or "all", and loads all variables named or all data variables.  Also needs to adjust calendar.
        siteTbl = QC_get_data(siteTbl, nc, datevec_cal(start_dnum,calendar), datevec_cal(end_dnum,calendar), outcal, maxmem, varNames, use_varNames, do_random, false, latlon_keepers);        
    end
    
    npts = siteTbl.Properties.UserData.npts;
    
    siteTbl.Properties.UserData = fix_field_order(siteTbl.Properties.UserData);
end

function siteTbl = reorder(siteTbl, stnID, searchType)

        % find the location in the output of each stnID requested.
    nsites = size(siteTbl,1);
    locs = nan(nsites,1);
    siteTbl.neworder = locs;    % add a temporary new column to siteTbl for sorting on.
    jx=0;
    if (isnumeric(stnID))
        for i=1:length(stnID)
            ix = find(siteTbl.index == stnID(i));
            for j=1:length(ix)
                if (~ismember(ix, locs(1:jx)))  % if it's not already in the reordering list, add it
                    jx=jx+1;
                    locs(jx)=ix(j);
                end
            end
        end
    else
        for i=1:length(stnID)        
            ix = find_site(stnID(i), siteTbl, string(searchType));
            for j=1:length(ix)
                if (~ismember(ix, locs(1:jx)))   % if it's not already in the reordering list, add it
                    jx=jx+1;
                    locs(jx)=ix(j);
                end
            end
        end
    end
    
    for i=1:nsites
        siteTbl.neworder(locs(i)) = i;
    end
    
    siteTbl = sortrows(siteTbl, {'neworder'});
    siteTbl.neworder = [];
    
end

function userdata = fix_field_order(userdata)
%   re-orders the UserData fields so they're easy to scan

    neworder = { 'varName', 'longName', 'units', 'nsites', 'npts', 'NAFlag', 'ncName', ...
                'dateRange', 'dateRangeLimit', 'calendar', 'timeunits',  'isUTC', 'dates', 'day1', 'modifiedDates', ...
                'isQCstntbl', 'allSites', 'fullYears','minYears', 'nctblNames', 'table_source'};
            
    fldnames = fieldnames(userdata);
    
    nf=length(neworder);
    for i=nf:-1:1
        ix=find(strcmp(neworder(i), fldnames),1);
        if (isempty(ix))
            neworder(i)=[];
        else
            fldnames(ix)=[];
        end
    end
    
    neworder = [neworder, fldnames{:}];
    userdata = orderfields(userdata, neworder);

end

function ix = find_site(ID, tbl, typ)

    if (contains(lower(typ),"id"))
        if (contains(typ,'like'))
            ix = find(contains(lower(tbl.stnID), lower(ID)));
        else
            ix = find(tbl.stnID == ID);
        end
    elseif (contains(lower(typ),"name"))
        if (contains(typ,'like'))
            ix = find(contains(lower(tbl.stnName), lower(ID)));
        else
            ix = find(tbl.stnName == ID);
        end
    elseif (contains(lower(typ),"like"))
        ix1 = find(contains(lower(tbl.stnID), lower(ID)));
        ix2 = find(contains(lower(tbl.stnName), lower(ID)));
        ix = unique(sort([ix1;ix2]));
    else
        ix1 = find(tbl.stnID == ID);
        ix2 = find(tbl.stnName == ID);
        ix = unique(sort([ix1;ix2]));
    end
end

function [ncName, start_date, end_date, fullYears, minYears, minPct, allSites, showDates, ...
          loadTblVars, stnID, searchType, loadData, maxmem, tryMatfile, latrange, lonrange, latlon_keepers_ix, varnames, use_varnames, outcal, do_random] = initParams(ncName, nargsout, varargin)

    search_types = ["stnID","stnName","any", "like_stnID", "like_stnName", "like_any","closest"];
    
    p = inputParser;
    p.StructExpand = true;
    p.KeepUnmatched = true;        

    flag = nargsout==0;
    
    addRequired(p,"ncName",                     @(s) ischar_s(s) || isQCstntbl(s) || isNCtbl(s));
    addOptional(p,"start_date",   [],           @(s) isempty(s) || isnumeric(s) && isrow(s) && any(length(s)==[1,3,6]));
    addOptional(p,"end_date",     [],           @(s) isempty(s) || isnumeric(s) && isrow(s) && any(length(s)==[1,3,6]));
    addParameter(p,"fullYears",   false,        @(s) islogical(s) || (isnumeric(s) && any(s==[0,1])));
    addParameter(p,"minYears",    0,            @(s) isnumeric(s) && s>=0);
    addParameter(p,"minPct",      0,            @(s) isnumeric(s) && s>=0 && s <= 100);
    addParameter(p,"showDates",  flag,          @(s) islogical(s) || (isnumeric(s) && any(s==[0,1])));
    addParameter(p,"loadTblVars", false,        @(s) islogical(s) || (isnumeric(s) && any(s==[0,1])));
    addParameter(p,"stnID",       strings(0),   @(s) isnumeric(s) || ischars(s));
    addParameter(p,"stations",    strings(0),   @(s) isnumeric(s) || ischars(s));
    addParameter(p,"latlons",     [],           @(s) isempty(s) || size(s,2)==2);
    addParameter(p,"searchType", "stnID",       @(s) ischar_s(s) && any(strcmpi(s,search_types)));       % numeric, for gridded
    addParameter(p,"loadData",    false,        @(s) islogical(s) || (isnumeric(s) && any(s==[0,1])));
    addParameter(p,"varnames",     strings(0),   @(s) isempty_s(s) || ischar_s(s));
    addParameter(p,"use_varnames",false,        @(s) islogical(s) || (isnumeric(s) && any(s==[0,1])));
    addParameter(p,"removeLeaps", [],           @(s) isempty(s) || islogical(s) || (isnumeric(s) && any(s==[0,1])));
    addParameter(p,"maxmem",      1e9,          @(s) isnumeric(s) && (s<=500 || s>=10e6));
    addParameter(p,"tryMatfile",  1,            @(s) islogical(s) || (isnumeric(s) && s>=0 && s<=2));
    addParameter(p,"latrange",    [],           @(s) isempty(s) || (isnumeric(s) && numel(s)==2 && s(1)>= -90 && s(1)<  90 && s(2)>s(1) && s(2)<=90));
    addParameter(p,"lonrange",    [],           @(s) isempty(s) || (isnumeric(s) && numel(s)==2 && s(1)>=-180 && s(1)< 180 && s(2)>s(1) && s(2)<=360));
    addParameter(p,"outcal",      strings(0),   @(s) isstring(s));       % NOT USED YET!
    addParameter(p,"calendar",    strings(0),   @(s) isstring(s));       % NOT USED YET!
    addParameter(p,"do_random",   false,        @(s) isempty(s) || islogical(s) || (isnumeric(s) && mod(s,1)==0)); 

    parse(p, ncName, varargin{:});
    Parms = p.Results;
    Unmatched = p.Unmatched;        % save rest of input params for later
    
    if (~isempty(fields(Unmatched)))
        flds=string(fields(Unmatched));
        msg = strcat(sprintf("error:  unexpected input keywords:  \n\t"), sprintf("%s ", flds));
        error(msg);
    end
    
    ncName       = Parms.ncName;
    start_date   = Parms.start_date;
    end_date     = Parms.end_date;
    fullYears    = Parms.fullYears;
    minYears     = Parms.minYears;
    minPct       = Parms.minPct;
    showDates    = Parms.showDates;
    loadTblVars  = Parms.loadTblVars;
    stnID        = string(Parms.stnID);
    stations     = string(Parms.stations);
    searchType   = Parms.searchType;
    loadData     = Parms.loadData;
    removeLeaps  = Parms.removeLeaps;
    maxmem       = Parms.maxmem;
    tryMatfile   = Parms.tryMatfile;
    latrange     = Parms.latrange;
    lonrange     = Parms.lonrange;
    latlons      = Parms.latlons;
    varnames     = Parms.varnames;
    use_varnames = Parms.use_varnames;
    outcal       = Parms.outcal;
    do_random    = Parms.do_random;
    latlon_keepers_ix = [];    % will find these later.
    
    if (~isempty(removeLeaps) && removeLeaps)
        outcal = "365-day";
    end
    if (isempty(outcal)), outcal = Parms.calendar; end
    
    if (maxmem < 500), maxmem = maxmem*1024*1024*1024; end  % if small, treat as GB spec.
    
        % make sure latrange and lonrange are doubles
    latrange = double(latrange);
    lonrange = double(lonrange);
    
    allSites    = ~(minYears>0 || ~isempty(start_date) || ~isempty(end_date) || minPct > 0 || ~isempty(stnID) || ~isempty(stations) || ~isempty(latlons));
     
    if (~isQCstntbl(ncName) && ~isNCtbl(ncName))
                % make sure we have full pathname to file
        [dir,~,~] = fileparts(char(ncName));
        if (isempty(dir)), ncName=fullfile(pwd(),ncName); end

        if (~isfile(ncName)), error("error:  file does not exist:  %s", ncName); end
    end
    
    if (strcmp(start_date,"all"))
        start_date = [];
        end_date = [];
    end
    if (length(start_date)==1), start_date = [start_date, 1, 1]; end
    if (length(  end_date)==1),   end_date = [  end_date,12,31]; end        % Q:  problem here if calendar is 360-day?
    
    if (minPct > 0 && minPct <=1)
        minPct = 100.0*minPct;
        fprintf(2,"warning:  minPct %.3f < 1.  Using %.1f for minPct", minPct/100, minPct);
    end
    
    if (length(stations) > 1 && isrow(stations)), stations = stations'; end
    if (length(stnID) > 1    && isrow(stnID)),    stnID    = stnID'; end
    if (~isempty(stations))
        if (isempty_s(stnID))
            stnID = stations;
        else
            stnID = [stnID(:);stations(:)];
        end
    end
    
        % get the list of station IDs
    if (~isempty_s(stnID))
                % is stnID the name of a file?
        if (ischar_s(stnID) && isfile(stnID)) 
            fname = stnID; 
            stnID = readtable(fname);
            if(isfield(stnID.Properties.VariableNames,"stnID"))
                stnID = stnID.stnID;
                searchType = "stnID";
            elseif (isfield(stnID.Properties.VariableNames,"stnName"))
                stnID = stnID.stnName;
                searchType = "stnName";       
            else    % no column stnID or stnName.  Reread file with no headers and use the first column.
                tbl = readtable(fname,'ReadVariableNames',false);
                stnID = tbl{:,1};
                if (~ischars(stnID)), error("error:  No stnID or stnName column, and 1st column is not strings");end
                fprintf("using 1st column of file %s for stnIDs", fname);
            end
        end

        if (istable(stnID))
            if(isfield(stnID.Properties.VariableNames,"stnID"))
                stnID = stnID.stnID;
                searchType = "stnID";
            elseif (isfield(stnID.Properties.VariableNames,"stnName"))
                stnID = stnID.stnName;
                searchType = "stnName";       
            else
                stnID = stnID{:,1};     % No stnID or stnName in table.  use 1st column.  
                if (~ischars(stnID)), error("error:  No stnID or stnName column, and 1st column is not strings");end
                fprintf("using 1st column (labeled %s) of %s for stnIDs", fname);            
            end
        end
    elseif (~isempty(latlons))
        [all_lats,all_lons] = ncdf_get_latlons(ncName);
        nlatlons = size(latlons,1);
        if (strcmp(searchType,"closest"))
            tol = .5;   % accept closest point within .5 degrees, about 55 km
        else                                   
            tol = .01;  % accept closest point within .01 degrees, about 1.1 km.  Assumes latlons are essentially correct.
        end
        [latlon_keepers_ix, llkeepers_in_tolerance] = closest_latlon(latlons, all_lats, all_lons, tol);
        if (sum(llkeepers_in_tolerance) ~= nlatlons), error("error:  cannot match latlons to lats & lons in file with tolerance of %.3 degrees", tol); end       
    end   
end

function [siteTbl, matname, ncName] = read_QC_matfile(ncName,  tryMatfile)
% returns siteTbl if ncName is a QC_site_table matfile.  
%   returns [] if not.
%       tryMatfile:  1      use matfile's siteTbl only if date/time stamp and size match
%                    2      use matfile's siteTbl even if timestamp, etc. doesn't match.

    [p,d,ext] = fileparts(ncName);
    siteTbl = [];
    if (~exist('tryMatfile', 'var') || isempty(tryMatfile)), tryMatfile=1; end  % (1= must 
    
    matname = fullfile(p,sprintf('%s%s',d,".mat"));
    if (strcmpi(ext,'.mat'))
        ncName  = fullfile(p, sprintf('%s%s',d,'.nc'));
    end
    if (~isfile(matname))       % if matfile doesn't exist, check canonical path as well.
        ncNameC = canonical_path(ncName);
        if (~strcmp(ncNameC, ncName))
            [p2,d2,~] = fileparts(ncNameC);
            matname  = fullfile(p2, sprintf('%s%s',d2,'.mat'));
        end
    end

    if (isfile(matname))
        try
                % look for siteTbl in matfile.
            minfo=whos(matfile(matname),'siteTbl');
            if (~isempty(minfo) && strcmp(minfo.class,'table'))
                load(matname,"siteTbl");
                finfo = dir_canonical(ncName);
                if (isfield(siteTbl.Properties.UserData,'table_source') && isequal(finfo, siteTbl.Properties.UserData.table_source))
                    return;
                end
                if (tryMatfile==1 && strcmpi(ext,".nc"))
                    fprintf("NOTE:  timestamp/size mismatch between matfile and netcdf file\n\t%s\n\t%s\nUsing nc file\n", matname, ncName);
                    siteTbl = [];
                end
            end
        catch
            siteTbl = [];
            if (strcmpi(ext,".nc"))
                fprintf("problem reading matfile %s.  Trying nc file %s instead", matname, ncName);
            else
                error("error reading %s", matname);
            end
        end
    elseif (strcmpi(ext,".mat"))        % if told read a matfile and it doesn't exist, error out.
        error("error reading %s", matname);
    end     
end

function [siteTbl,nc] = read_netcdf_QC_table(ncName, tryMatfile)
%   Reads full list of stations from netcdf file and returns a QC_stntbl.
%   Does not load any data variables.
%   This could probably be sped up considerably, Ian, with single call to ncdf once list of variables needed is
%   determined.
%   


    ncName = canonical_path(ncName);
    [~,~,ext] = fileparts(ncName);
    
    if ((tryMatfile) || strcmpi(ext,".mat"))  % reading matfile is much faster, if it's available.
        
        [siteTbl, ~, ncName] = read_QC_matfile(ncName, tryMatfile);
        nc=[];
        
        if (~isempty(siteTbl)), return; end  % if unsuccessful, then try reading netcdf file, which is slower.
    end
    
    [yn, nc] = isQCnetcdf(ncName);
    if (~yn), error("error:  %s is not a QC_site_table", ncName); end
%   nc = ncdf(ncName);
%   vars = ["Tmax","Tmin","Prec","tasmax","tasmin","pr", "temp_F","rh_F","temperature","rhsmax","rhsmin","relhum"];
%   vars = [  "tmax","tasmax",  "tmin","tasmin","tavg","temp_F","temp_C", "temperature", "Prec","prcp","precipitation",  "pr","precip","rh","rh_F","rhsmax","rhsmin","relhum"];

        % get first climate variable name from the file (excluding zvals variable)
    is_climvar=is_climate_variable(nc.varlist()) & ~is_zvals_variable(nc.varlist);
    vix = find(is_climvar,1);
    if (~isempty(vix))
        varName = nc.Variables(vix).Name;
    else
        warning("warning:  cannot find a climate variable in %s\n", nc.FileName);
    end
    
%     varName="";
%     for i=1:length(vars)
%         vix = find(strcmpi(vars(i), nc.varlist()),1);
% %         if (isempty(vix))
% %             vix = find(contains(nc.varlist(),"zvals"),1);     % we don't want to return the zvals variable
% %         end
%         if (~isempty(vix))
%             varName = nc.Variables(vix).Name;
%             break;
%         end
%     end   
% 
    if (strlength(varName)~=0)

        ix = find(strcmp({nc.Dimensions.Name},'time'),1);
        npts = nc.Dimensions(ix).Length;
%         timeunits = ncreadatt(ncName, 'time','units');   % should be:  "days since 1850-01-01"
        timeunits = nc.getattvalue('time/units');
        [day1vec, timescale, isUTC] = nc_parse_date_str(timeunits);
        calendar = nc.getattvalue('time/calendar');
        day1 = datenum_cal(day1vec, calendar);
%         dates = day1 + ncread(ncName,'time')/timescale;    
%         dates = day1 + nc.loadvar('time')/timescale;  
        tstamps = nc.getvardata('time');
        if (isempty(tstamps))
            tstamps = nc.loadvar('time');
        end
        dates = day1 + tstamps/timescale;
        var   = nc.get(varName);
        try
            longName= var.getattvalue('long_name');
        catch
            longName = '';
        end
        try
            NAFlag  = var.getattvalue('_FillValue');
        catch
            NAFlag=[];
        end
        try
            varUnits= var.getattvalue('units');
        catch
            varUnits='';
        end
    else
                % not a QC file with variables.  
        error("error:  cannot find a valid variable name in the netcdf file"); 
    end

      tblNames = ["stnID","lat","lon","stnName","elev",     "startDate", "endDate", "pctValid", "time_zone","is_ghcn","nearest_ghcn","ghcn_distance","ghcn_az","ghcn_elev","index"];    % names in the table
    nctblNames = ["stnID","lat","lon","stnName","elevation","start_date","end_date","pct_valid","time_zone","is_ghcn","nearest_ghcn","ghcn_distance","ghcn_az","ghcn_elev","index"];    % equivalent names in the netcdf file.
    ntbl = length(tblNames);
    stringNames = ["stnID","stnName","nearest_ghcn"];         % shouldn't have tz_name anymore, but some old files may still have that columns.

            % read QC table variables from nc file
    vals = cell(ntbl,1);
    varkeepers = true(ntbl,1);         % flags which variables were in the netcdf file.
    for i=1:ntbl
        vname = char(nctblNames(i));
        try
%             vals{i} = ncread(ncName, vname);
            vals{i} = nc.readvar(vname);
            if (ismember(vname, stringNames)) 
                vals{i} = strtrim(string(vals{i}')); 
            end
            if (i==1), nsites = length(vals{i}); end
        catch
            vals{i} = [];
            varkeepers(i) = false;
        end
    end

            % put variables into siteTbl
    siteTbl = table(vals{1}, 'VariableNames',{'stnID'});
    for i=2:ntbl
        vname = char(tblNames(i));
        if (varkeepers(i))
            try
            siteTbl.(vname) = vals{i};
            catch
                fprintf(2, "warning:  problem inserting variable %s into table", vname);
            end
        end
    end

        % startDate and endDate are in netcdf file as days-since-day1.  Convert them to matlab datenums.
    siteTbl.index     = (1:nsites)';
    vnames=["startDate","endDate"];
    for i=1:length(vnames)
        vname=char(vnames(i));
        if (ismember(vname, siteTbl.Properties.VariableNames))
            siteTbl.(vname) = siteTbl.(vname) + day1;
        end
    end

    siteTbl.Properties.UserData.timeunits = timeunits;
    siteTbl.Properties.UserData.calendar = calendar;
    siteTbl.Properties.UserData.day1 = day1;
    siteTbl.Properties.UserData.isUTC = isUTC;
    siteTbl.Properties.UserData.isQCstntbl=true;       % flags this as a QC stn table.
    siteTbl.Properties.UserData.dates = dates;    
    siteTbl.Properties.UserData.fileDateRange=sprintf("%s ", string(datestr_cal([dates(1); dates(end)],calendar,"yyyy-mm-dd")));
    siteTbl.Properties.UserData.dateRange = siteTbl.Properties.UserData.fileDateRange;
    nsites = size(siteTbl,1);
    siteTbl.Properties.UserData.ncName = ncName;                    % capturing output.  add Userdata info to table
    siteTbl.Properties.UserData.varName = varName;
    siteTbl.Properties.UserData.npts = npts;
    siteTbl.Properties.UserData.NAFlag = NAFlag;
    siteTbl.Properties.UserData.longName = longName;
    siteTbl.Properties.UserData.units = varUnits;
    siteTbl.Properties.UserData.nsites = nsites;
    
    table_source = dir_canonical(ncName);
    siteTbl.Properties.UserData.table_source=table_source;
    
    nctblNames = nctblNames(varkeepers);
    siteTbl.Properties.UserData.nctblNames = nctblNames;
    
    siteTbl.Properties.UserData.globalAttributes = ncdf_get_attributes(nc);
    siteTbl.Properties.UserData.RunParams = ncdf_get_attributes(nc,"RunParams");
    siteTbl.Properties.UserData.DownscalingParams = ncdf_get_attributes(nc,"DownscalingParams");
    
end  

function [siteTbl, nc] = load_table_vars(siteTbl, nc, loadTblVars)
% loads any non-standard 1-D variable from the table which uses stnnum as the dimension. 
%   or any character variable with stnnum as the 2nd dimension and 1st dimension < 256. 

    if (islogical(loadTblVars) && ~loadTblVars), return; end
    
    if (isempty(nc)), nc=ncdf(siteTbl.Properties.UserData.ncName); end
    
    nctblNames = siteTbl.Properties.UserData.nctblNames;
    varName = siteTbl.Properties.UserData.varName;
    
    nvars = length(nc.Variables);
    for i=1:nvars
        v=nc.Variables(i);
        if (strcmpi(v.Name, varName) || any(strcmpi(v.Name,nctblNames))), continue; end    % we already loaded these.
        [isStnVar, isStnText] = is_station_variable(nc.Variables(i));
        if (isStnVar)
            vname=v.Name;
            vals=nc.readvar(v.Name);
            if (isStnText)
                vals = strtrim(string(vals'));
            end
            try
                siteTbl.(vname) = vals;
            catch
                fprintf("oops.  problem adding variable %s to table\n", vname);
            end
        end
    end    
end

function  [isStnVar, isStnText] = is_station_variable(ncVariable)
% flags if ncVariable is a data variable tied to the stn_num dimension.
% returns true if stn_num is ncVariable's single dimensions || stn_num is the 2nd dimension and the variable is a char
% variable.  (but false if length of strings > 256, in case it's a char data variable)
    stn_num_dim = 0;
    ndims = length(ncVariable.Dimensions);
    dimsize=ncVariable.Size;
    istext=strcmp(ncVariable.Datatype,'char');
    for i=1:ndims
        if (strcmp(ncVariable.Dimensions(i).Name,'stn_num'))
            stn_num_dim = i;
            break;
        end
    end
    isStnText = (istext && ndims==2 && stn_num_dim == 2 && dimsize(1) <= 256);
    isStnVar  = ( (stn_num_dim == 1 && ndims == 1) || isStnText);
end

