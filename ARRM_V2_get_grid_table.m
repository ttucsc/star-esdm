function [ncTbl, nsites, varName, npts] = ARRM_V2_get_grid_table(ncNames, varargin)
%function [siteTbl, nsites, varName, npts] = QC_get_site_table(ncName, startDate, endDate, [key/value pairs])
%
% THIS NEEDS: some sort of date/timestamp and/or hash for verifying table is uptodate if ncName is an ARRM_V2_grid_table or .mat file
%
%   Inputs:
%       ncName          name of station netcdf file, or .mat file w/ contents of QC_stntbl, or  QC_stntbl from a
%                           previous call to QC_get_site_table
%     optional inputs:%                           
%       startDate      year, or [year,month,day[,hr,min,sec]] of starting date                                 [empty] 
%       endDate        year, or [year,month,day[,hr,min,sec]] of ending date                                   [empty]
%                           If specified, include only stations with data in date range specified
%                           (see also fullYears and minYears option below)
%                           If startDate empty or missing, use file's starting date as startDate
%                           If endDate empty or missing, use file's ending date as endDate.
%
%     optional keyword/value pairs
%       "showDates", t/f    if true, dates returned as datestr(...,'yyyy-mm-dd');  else dates returned as datenums [false]                                                  [empty]
%       loadTblVars         if true, loads all netcdf variables (except actual data) which use station as a dimension.  
%                                   i.e., all variables associated with station, except for data.               [false]
%                                   if false, only loads the standard QC station file columns.
%       "calendar", outcalendar  changes dates from original calendar to new calendar, inserting or deleting days
%                                   to switch between 360-, 365- and 365.25-day year lengths.
%       "latrange", [lat1,lat2] keep only locations in half-open latitude  range [lat1,lat2);
%       "lonrange", [lon1,lon2] keep only locations in half-open longitude range [lon1,lon2);
%
%
%   Outputs
%       gridTbl         grid table, one line for each lat/lon location
%       nsites          # of sites in table
%       varName         name of data variable found in netcdf file
%       npts            number of datapoints (one for each time step) in file.
%                           
%---------------------QC Table output----------------
%   ARRM_V2_grid_table (output) columns:
%           lat, lon        location  (decimal degrees, usually to 4 decimals;  some only to 2 or 3)
%           startDate       start timestamp of first valid record for location, as Matlab datenum (for data's calendar type)
%           endDate         end   timestamp of last  valid record for location, as Matlab datenum (for data's calendar type)
%           pctValid        percent of valid points between startDate and endDate
%           index           [i,j] index of record in original input file
%           stnID           string, "(%8.4f,%9.4f)" , text string identifying lat & lon of center of gridcell
%           stnName         string, "[%4d,%4d]", text string identifying the (i,j) index of gridcell in original file.
%
%       possible additional fields
%
%   QC table Properties.UserData fields & example values:
%       Properties.UserData contains useful table-wide data extracted from netcdf file's global attributes.
%
%         ncName: "/Volumes/jcsf_data/data/obs/stations_netcdf/stations.Tmax.ncamerica.1850.2018.20_yrs_min.nc"
%                                                       original input file, full pathname
%        varName: 'TMAX'                                variable name
%           npts: 24106                                 # of data points (dates) in file.  All sites have same npts.
%                                                       date for any station are set to NAs in the netcdf file.
%   allLocations: false                                 boolean.  flags whether all sites from netcdf file included
%                                                           true:   all sites included;  
%                                                           false:  only some sites from original file included.
%  modifiedDates: false                                 boolean.  flags whether date information has been modified from
%                                                           original source.  I.e., table has subset of original dates
%                                                           true:   date info is subset of date info in original nc file 
%                                                           false:  date info is identical to date info in orig. nc file
%   missingDates: false                                 boolean.  true if any dates are missing in original file.
%         NAFlag: 1.0000e+20                            file's original NA value ('missing or invalid data' flag)
%       longName: 'daily max temp'                      variable long name
%          units: 'degC'                                variable units
%         nsites: 70                                    # of grid locations being returned (called sites to match
%                                                           QC_netcdf talbes
%      timeunits: 'days since 1850-01-01 00:00:00 UTC'  netcdf time units
%       calendar: '365-day'                             calendar type.  365-day, 360-day, standard, etc.
%   filecalendar: '365-day'                             original calendar type of file
%           day1: 675700                                matlab datenum_cal of  'days since' day of timeunits
%       fileday1: 675700                                matlab datenum_cal of  'days since' day of file's original timeunits
%          isUTC: 1                                     boolean. flags whether timestamps are UTC or local time.
%        isNCtbl: 1                                     boolean.  true if table is a grid table.
%          dates: [24090x1 double]                      datenums_cal for each data value (i.e., matlab datenums based on calendar).
%      filedates: [24090x1 double]                      datenums_cal for original dates in file.
%  fileDateRange: "1950-01-01 2015-12-31"               date range of data in original file
%      dateRange: "1950-01-01 2015-12-31"               date range of requested data (as specified by start_ & endDate)
%           lats: [180x1 double]                        1-D vector of latitudes in table
%           lons: [360x1 double]                        1-D vector of longitudes in table
%       fileLats: [180x1 double]                        1-D vector of original lats in file
%       fileLons: [360x1 double]                        1-D vector of original lons in file
%       
%    tableSource: [dir struct]                          directory info struct of original file
%     nctblNames: [1x1 string]                          string array of data variables
% globalAttributes: [5x1 cell array]                    cell array of global attributes
%      RunParams: [20x1 cell array]                     cell array of ARRM_V2 RunParams object used to downscale data
%                                                           (empty if not ARRM_V2 downscaled output)
% DownscalingParams" 30x1 cell array]                   cell array of ARRM_V2 downscaling params object used to
%                                                           downscale data
%                                                           (empty if not ARRM_V2 downscaled output)
%-----------------------------------

% steps for gridded:

%   1.  Read file, make basic table  (add outcal)
%           if tryMatFile, 
%               read from matfile if timestamp matches.
%           else
%               check it's QCnetcdf (reads nc in as well)
%               get varname from nc
%               get time variable
%                   timeunits, calendar
%                   day1
%                   read dates
%               get var attributes:  name, longname, NAflag, units, variable-order (LaLoT, TLaLo  TLoLa, etc.)
%               read table vars
%                   lat, lon, generate (stnID), (stnName) (elev) (index as (start,step,stride)  (maybe others?  elevation, time_zone)
%               Create table
%                   add variables
%                   add index
%           Fill in UserData
%               timeunits, calendar, day1, isUTC, isARRMtbl, isQCstntbl/false,dates,daterange
%               ncname,varname,npts,NAflag,longname,varUnits,nsites
%               nctblNames(varkeepers)
%               globalAttributes
%               RunParams (if present)
%               DownscalingParams (if present)
%               

%   2.  load additional table vars if present
%           open as ncdf
%           get nctblNames from UserData
%           for each variable in ncdf:
%               if in nctblNames, skip it (we already read it)
%               if it uses lat/lon as dimensions AND
%                   no other dimension, or ischar and lat/lon are 2nd, 3rd dimensions
%                   read data
%                   add variable as column

%   3.  Get metadata:
%           varname
%           calendar, dates, start & end dates for sites.  (will be same for all gridded locations.
%           limit file to locations specified
%               either by indexes or by lat/lon boundaries.
%           no neeed to limit by 
%               start/end date.
%               pct valid
%               minyears
%               fullyears
%
%   4.  Set/update UserData
%           allsites, nsites, modifiedDates (rename to limitedDates
%
%   5.  If showDates, convert dates to text.
%   6.  If loadData, read data from file.
%   7.  re-order UserData fields.

%   startDate, endDate as datenums(outcal);
%   dates as datenums(outcal);

Not sure if this has been completed or tested.  Ian 8/5/22

    ncNames = string(ncNames);
    nfiles = length(ncNames);
    
    startDate = [];
    endDate = [];
    if (length(varargin)>=1)
        if (isnumeric(varargin{1}))
            startDate = varargin{1};
            varargin=varargin(2:end);
            if (length(varargin)>=1 && isnumeric(varargin{1}))
                endDate = varargin{1}; 
                varargin=varargin(2:end);
            end
        end
    end
    
    if (nfiles == 1)
        [ncTbl, nsites, varName, npts] = ARRM_V2_get_grid_table_sub(ncNames(1), startDate, endDate, varargin{:});        
        
    else
        [ncTbl, startDate, endDate, varName, loadData, outcalendar] = scan_files(ncNames, startDate, endDate, varargin{:}); 
    
        if (loadData)
            [ncTbl, nsites, varName, npts]  = join_data(ncNames, ncTbl, startDate, endDate, varargin{:});     
        else
            nsites = size(ncTbl,1);
            npts = datenum_cal(endDate,outcalendar) - datenum_cal(startDate, calendar) + 1;
        end
    end
end

function [ncTbl, startDate, endDate, varName, loadData, outcalendar] = scan_files(ncNames, startDate, endDate, varargin)
            
    nfiles = length(ncNames);
        % get info from all files.
        
%   [ncName, startDvec, endDvec, showDates, loadTblVars, loadData, ...
%    varName, outcalendar, do_random, maxmem, tryMatfile, latrange, lonrange] = initParams(ncNames(1), startDate, endDate, varargin{:});
    [~,      ~,         ~,       ~,         ~,            loadData, ...
     varName, outcalendar, ~,         ~,      tryMatfile, latrange, lonrange] = initParams(ncNames(1), startDate, endDate, varargin{:});
 
        
    for n=1:nfiles
        ncName = ncNames(n);

        tbl = ARRM_V2_get_grid_table_sub(ncNames(n), startDate, endDate, 'showDates', false, ...
                                        'loadTblVars', false, 'loadData', false, ...
                                         'tryMatfile', tryMatfile, ...
                                         'latrange', latrange, 'lonrange', lonrange);
        if (n==1)
                % get info from first file
            varName     = tbl.Properties.UserData.varName;
            incalendar  = tbl.Properties.UserData.calendar;
            units       = tbl.Properties.UserData.units; 
            start_date  = tbl.Properties.UserData.dates(1);
            end_date    = tbl.Properties.UserData.dates(end);
            timeunits   = tbl.Properties.UserData.timeunits;
            day1        = tbl.Properties.UserData.day1;
            ncTbl = tbl;
        else
            start_date = min(start_date, tbl.Properties.UserData.dates(1));
            end_date   = max(end_date,   tbl.Properties.UserData.dates(end));
            ncTbl = join_tables(ncTbl, tbl);
                % get info from file
            varName2     = tbl.Properties.UserData.varName;
            calendar2    = tbl.Properties.UserData.calendar;
            units2       = tbl.Properties.UserData.units; 
            timeunits2   = tbl.Properties.UserData.timeunits;
            day12        = tbl.Properties.UserData.day1;
            [~,fn,ext]   = fileparts(ncName);        
                % check for errors
            if (calendar_length(incalendar) ~= calendar_length(calendar2))
                error("calendar  mismatch %s %s %s", incalendar, calendar2,  sprintf("%s%s", fn, ext)); 
            end
            if (~strcmp(varName,   varName2)),    error("varName    mismatch %s %s %s", varName,   varName2,   sprintf("%s%s", fn, ext)); end 
            if (~strcmp(units,     units2)),      error("units     mismatch %s %s %s", units,      units2,     sprintf("%s%s", fn, ext)); end
            if (~strcmp(timeunits, timeunits2)),  error("timeunits mismatch %s %s %s", timeunits,  timeunits2, sprintf("%s%s", fn, ext)); end
            if (day1 ~= day1),                    error("day1      mismatch %d %d %s", day1,       day12,      sprintf("%s%s", fn, ext)); end
        end
    end
    if (isempty(startDate))
        startDate = datevec_cal(start_date, incalendar); 
    end
    if (isempty(  endDate))
        endDate   = datevec_cal(end_date, incalendar);
    end

end

function [ncTbl] = join_tables(ncTbl, tbl)

    [~, new_ix] = setdiff(tbl(:,{'lon','lat'}), ncTbl(:,{'lon','lat'}));
    
    ncTbl.Properties.UserData.ncName = [string(ncTbl.Properties.UserData.ncName); string(tbl.Properties.UserData.ncName)];
    if (~isempty(new_ix))
        ncTbl = [ncTbl; tbl(new_ix,:)];    
    end
end

function [ncTbl, nsites, varName, npts] = join_data(ncNames, ncTbl, startDate, endDate, varargin)

    [~, startDvec, endDvec, showDates, loadTblVars, loadData, ...
     varName, outcalendar, do_random, maxmem, tryMatfile, latrange, lonrange] = initParams(ncNames(1), startDate, endDate, varargin{:});
 
    sdnum = datenum_cal(startDate, outcalendar);
    ednum = datenum_cal(  endDate, outcalendar);
    
    dates = (sdnum:ednum)';
    nsites = size(ncTbl,1);
    npts = length(dates);
    ncTbl.Properties.UserData.dates = dates;
    ncTbl.Properties.UserData.calendar = outcalendar;
    ncTbl.Properties.UserData.modifiedDates = true;
    ncTbl.Properties.UserData.npts = npts;
    ncTbl.Properties.UserData.nsites = nsites;
    ncTbl.Properties.UserData.dateRange = sprintf("%s %s", datestr_cal(sdnum, outcalendar), datestr_cal(ednum, outcalendar));
    ncTbl.data = nan(nsites,npts);
    
    for n=1:length(ncNames)
        tbl = ARRM_V2_get_grid_table_sub(ncNames(n), startDvec, endDvec, 'showDates', showDates, ...
                                        'loadTblVars', loadTblVars, 'loadData', loadData, ...
                                         'calendar',outcalendar, 'do_random',do_random, ...
                                         'maxmem',maxmem, 'tryMatfile', tryMatfile, ...
                                         'latrange', latrange, 'lonrange', lonrange);
%         UD = tbl.Properties.UserData;
%         ud1 = find(UD.dates >= sdnum, 1);
%         ud2 = find(UD.dates <= ednum, 1, 'last');

        mydata = nan(nsites,npts);
%         ix1 = UD.dates(ud1) - sdnum + 1;
%         ix2 = UD.dates(ud2) - sdnum + 1;                
        [Lia, Locb] = ismember(tbl.stnID, ncTbl.stnID);
        if (any(~Lia)), error("oops.  some sites not found!"); end
            % move the new data into the right locations in a temporary array
        mydata(Locb, :) = tbl.data;
            % find the valid data
        keepers = ~isnan(mydata);
            % copy the valid data into the output array.
        ncTbl.data(keepers) = mydata(keepers);
    end
    
        % update the start,end dates.
 
    okflags = ~isnan(ncTbl.data);
    missingData = sum(~okflags(:))>0;
    ncTbl.Properties.UserData.missingData = missingData;
    if (any(~okflags(:,1)))
        for i=1:nsites
            ix1 = find(okflags(i,:),1);
            ncTbl.startDate(i) = dates(ix1);
        end
    else
        ncTbl.startDate = repmat(dates(1),nsites,1);
    end
    if (any(~okflags(:,end)))
        for i=1:nsites
            ix2 = find(okflags(i,:),1,'last');
            ncTbl.endDate(i) = dates(ix2);
        end
    else
        ncTbl.endDate = repmat(dates(end),nsites,1);
    end
end

function [ncTbl, ngrids, varName, npts] = ARRM_V2_get_grid_table_sub(ncName, varargin)

    [ncName, startDvec, endDvec, showDates, loadTblVars, loadData, ...
     varNames, outcalendar, do_random, maxmem, tryMatfile, latrange, lonrange] = initParams(ncName, varargin{:});
        % make sure filename has path if in local folder, so we can put full path info into UserData
        
    if (isNCtbl(ncName))
        ncTbl = ncName;
%       ncName = ncTbl.Properties.UserData.ncName;
    else
        [ncTbl, nc]  = read_netcdf_table(ncName, tryMatfile);        % will also read .mat file if ncName is matfile. 
    end
    
    if (loadTblVars)
         ncTbl = load_table_vars(ncTbl, nc, loadTblVars);
%       siteTbl = load_table_vars(siteTbl, nc, loadTblVars, siteTbl.Properties.UserData.nctblNames, siteTbl.Properties.UserData.varName);
    end   
        
    varName      = ncTbl.Properties.UserData.varName;
    calendar     = ncTbl.Properties.UserData.calendar;
    day1         = ncTbl.Properties.UserData.day1;
    dates        = ncTbl.Properties.UserData.dates;
    fileDates    = ncTbl.Properties.UserData.fileDates;
    missingDates = ncTbl.Properties.UserData.missingDates;
    lats         = ncTbl.Properties.UserData.fileLats;
    lons         = ncTbl.Properties.UserData.fileLons;
    
    if (calendar_length(calendar)==360 && endDvec(2)==12 && endDvec(3)==31), endDvec(3)=30; end % change Dec 31 to Dec 30 for 360-day calendar
    
        % we need these as datenums, not datevecs, but needed to know the calendar before we could convert them.
    modifiedDates = false;
    if (~isempty(startDvec))
        start_dnum = datenum_cal(startDvec,calendar);
        if (start_dnum ~= fileDates(1))
            modifiedDates = true;
            day1 = min(day1, start_dnum);
        end
    else
        start_dnum = dates(1);
    end
    if (~isempty(endDvec))
        end_dnum = datenum_cal(  endDvec,calendar); 
        if (end_dnum ~= fileDates(end))
            modifiedDates = true;
        end
    else
        end_dnum = dates(end);
    end
    
    if ((start_dnum > fileDates(end) || end_dnum < fileDates(1)))
        error('no valid data in date range %s to %s', datestr_cal(start_dnum,calendar,'yyyy-mm-dd'),datestr_cal(end_dnum,calendar,'yyyy-mm-dd')); 
    end
    
        % if we've modified the dates, check for missing dates only in the range specified.
        % we check the actual dates, in case the data is all there but is out of order, which we treat as
        % "missingDates"
    if (modifiedDates)
        dates = (start_dnum:end_dnum)';
        ncTbl.Properties.UserData.day1 = day1;
        ncTbl.Properties.UserData.dates = dates;
        ix1 = find(fileDates >= start_dnum,1);
        ix2 = find(fileDates <= end_dnum,1,'last');
        if (~isequal(fileDates(ix1:ix2), dates)), missingDates = true; end
    end
    
        % in case start, end dates are in text strings instead of datenums.
    if (ischars(ncTbl.startDate)), ncTbl.startDate = datenum_cal(datevec(ncTbl.startDate), calendar); end
    if (ischars(ncTbl.endDate)),   ncTbl.endDate   = datenum_cal(datevec(ncTbl.endDate), calendar); end
    
    if (start_dnum ~= fileDates(1))
        ncTbl.startDate = max(ncTbl.startDate, start_dnum);
    end
    if (end_dnum ~= fileDates(end))
        ncTbl.endDate = min(ncTbl.endDate, end_dnum);
    end
        
                % limit sites to lat/lon region if lats or lons specified.

    allSites = true;
    if (~isempty(latrange) || ~isempty(lonrange))
        if (~isempty(latrange))
            latkeepers = ncTbl.lat >=latrange(1) & ncTbl.lat <= latrange(2);
            if (sum(latkeepers)==0), error("error:  no gridcells in lat range %.4f - %.4f)", latrange); end 
            allSites = all(latkeepers);
            lats = lats(lats >= latrange(1) & lats <= latrange(2));
        end        

        if (~isempty(lonrange))
            lonkeepers = ncTbl.lon >=lonrange(1) & ncTbl.lon <= lonrange(2);
            if (sum(lonkeepers)==0), error("error:  no sites in lon range %.4f - %.4f)", lonrange); end 
            allSites = allSites & all(lonkeepers);
            lons = lons(lons >= lonrange(1) & lons <= lonrange(2));
        end
        keepers = latkeepers & lonkeepers;
        ncTbl = ncTbl(keepers,:);
    end

    ngrids = size(ncTbl,1);
        
%     npts = endDate - startDate + 1;    
%     siteTbl.Properties.UserData.npts = npts;
    ncTbl.Properties.UserData.allSites = allSites;
    ncTbl.Properties.UserData.nsites = ngrids;
    ncTbl.Properties.UserData.lats = lats;
    ncTbl.Properties.UserData.lons = lons;
    ncTbl.Properties.UserData.modifiedDates = modifiedDates;
    ncTbl.Properties.UserData.missingDates = missingDates;
%   ncTbl.Properties.UserData.fileDateRange = sprintf("%s %s", string(datestr(datevec_cal([fileDates(1); fileDates(end)],calendar),'yyyy-mm-dd')));
    ncTbl.Properties.UserData.dateRange = sprintf("%s %s", string(datestr(datevec_cal([start_dnum; end_dnum],calendar),'yyyy-mm-dd')));
    
    [ncTbl.stnID, ncTbl.stnName] = create_ids_and_names(ncTbl);
    
    if (showDates)
        try
            ncTbl.startDate = string(datestr_cal(ncTbl.startDate,calendar, "yyyy-mm-dd"));
            ncTbl.endDate   = string(datestr_cal(ncTbl.endDate,  calendar, "yyyy-mm-dd"));
        catch
        end
    end
    
    if (loadData)
        ncTbl = gridded_get_data(ncTbl, datevec_cal(start_dnum,calendar), datevec_cal(end_dnum,calendar), "varNames", varNames, "calendar", outcalendar, "ncdf", nc, "do_random", do_random, "maxmem", maxmem);
        npts = size(ncTbl.data,2);
    else
        npts = [];
    end
    
    
    ncTbl.Properties.UserData = fix_field_order(ncTbl.Properties.UserData);
end

function userdata = fix_field_order(userdata)
%   re-orders the UserData fields so they're easy to scan

    neworder = { 'ncName', 'varName', 'longName', 'units', 'NAFlag', 'nsites', 'allSites', 'npts', ...
                    'calendar', 'timeunits', 'day1', 'isUTC',     'dateRange',    'dates', 'modifiedDates', ...
                'fileCalendar',                               'fileDateRange', 'fileDates','missingDates', ...
                'lats','lons','fileLats','fileLons', ...
                'isNCtbl', 'isQCtbl', 'nctblNames', 'table_source'};
            
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

function [ncName, startDvec, endDvec, showDates, loadTblVars, loadData, ...
          varNames, outcalendar, do_random, maxmem, tryMatfile, latrange, lonrange] = initParams(ncName, varargin)
    
    p = inputParser;
    p.StructExpand = true;
    p.KeepUnmatched = false;     % throw error if any unexpected input parameters.

    
    addRequired(p,"ncName",                 @(s) ischar_s(s) || isQCstntbl(s));
    addOptional(p,"startDate",    [],       @(s) isempty(s) || isnumeric(s) && isrow(s) && any(length(s)==[1,3,6]));        % take either version...
    addOptional(p,"endDate",      [],       @(s) isempty(s) || isnumeric(s) && isrow(s) && any(length(s)==[1,3,6]));
    addParameter(p,"showDates",  false,     @(s) islogical(s) || (isnumeric(s) && any(s==[0,1])));
    addParameter(p,"loadTblVars", false,    @(s) islogical(s) || (isnumeric(s) && any(s==[0,1])));
    addParameter(p,"loadData",    false,    @(s) islogical(s) || (isnumeric(s) && any(s==[0,1])));
    addParameter(p,"varNames",    "",       @(s) ischars(s));
    addParameter(p,"calendar",    "",       @(s) ischar_s(s));      % output calendar
    addParameter(p,"do_random",   true,     @(s) islogical(s) || (isnumeric(s) && s>=1));
    addParameter(p,"maxmem",      1e9,      @(s) isnumeric(s) && s>=10e6);
    addParameter(p,"tryMatfile",  1,        @(s) islogical(s) || (isnumeric(s) && s>=0 && s<=2));
    addParameter(p,"latrange",    [],       @(s) isempty(s) || (isnumeric(s) && numel(s)==2 && s(1)>= -90 && s(1)<  90 && s(2)>s(1) && s(2)<=90));
    addParameter(p,"lonrange",    [],     @(s) isempty(s) || (isnumeric(s) && numel(s)==2 && s(1)>=-180 && s(1)< 360 && s(2)>s(1) && s(2)<=360));

    parse(p, ncName, varargin{:});
    Parms = p.Results;
    
    ncName       = Parms.ncName;
    startDvec   = Parms.startDate;
    endDvec     = Parms.endDate;
    showDates    = Parms.showDates;
    loadTblVars  = Parms.loadTblVars;
    loadData     = Parms.loadData;
    varNames     = Parms.varNames;
    outcalendar  = Parms.calendar;
    do_random    = Parms.do_random;
    maxmem       = Parms.maxmem;
    tryMatfile   = Parms.tryMatfile;
    latrange     = Parms.latrange;
    lonrange     = Parms.lonrange;

    if (~isempty(lonrange))    
        if (all(lonrange==[0,360]))
            lonrange = [-180,180];
        elseif (all(lonrange~=[-180,180]))
            lonrange = mod(lonrange+180,360)-180;   % make sure lonrange is in [-180,180], not [0,360].   problem here specifying lon range crossing 0 or 180...
            if (lonrange(2) == -180), lonrange(2) = 180; end  
        end
    end
        
    if (ischar_s(ncName))
                % make sure we have full pathname to file
        [dir,~,~] = fileparts(char(ncName));
        if (isempty(dir)), ncName=fullfile(pwd(),char(ncName)); end

        if (~isfile(ncName)), error("error:  file does not exist:  %s", ncName); end
    else
        error("ncName is not a filename or a valid ncTbl");
    end
    
    if (strcmp(startDvec,"all"))
        startDvec = [];
        endDvec = [];
    end
    if (length(startDvec)==1), startDvec = [startDvec, 1, 1]; end
    if (length(  endDvec)==1),   endDvec = [  endDvec,12,31]; end
    
end

function [siteTbl, matname, ncName] = read_netcdf_matfile(ncName,  tryMatfile)
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

function [ncTbl,nc] = read_netcdf_table(ncName, tryMatfile)
%   Reads full list of stations from netcdf file and returns an NCtbl.
%
%   This could probably be sped up considerably, Ian, with single call to ncdf once list of variables needed is
%   determined.


    ncName = canonical_path(ncName);
    [~,~,ext] = fileparts(ncName);
    ncTbl = [];
    
    if (tryMatfile || strcmpi(ext,".mat"))  % reading matfile is much faster, if it's available.
        
        [ncTbl, ~, ncName] = read_netcdf_matfile(ncName, tryMatfile);
        if (~isempty(ncTbl))
            fn = ncTbl.Properties.UserName.ncName;
            if (~isfile(fn)), error("error:  netcdf file %s doesn't exist", fn); end
            nc=ncdf(fn);
        end
    end
    
    if(isempty(ncTbl))
    
        [yn, nc] = isnetcdf(ncName);
        if (~yn), error("error:  %s is not a netcdf file %s can work with", ncName, mfilename); end

        vars = ["Tmax","Tmin","Prec",'tasmax','tasmin','pr', "temp_F","rh_F","temperature","rhsmax","rhsmin","relhum"];
        varName="";
        for i=1:length(vars)
            vix = find(strcmpi(vars(i), {nc.Variables.Name}),1);
            if (~isempty(vix))
                varName = string(nc.Variables(vix).Name);
                break;
            end
        end   

        if (strlength(varName)~=0)


    %       ix = find(strcmp({nc.Dimensions.Name},'time'),1);
    %       npts = nc.Dimensions(ix).Length;
    %       timeunits = ncreadatt(ncName, 'time','units');   % should be:  "days since 1850-01-01"
            timeunits = nc.getattvalue('time/units');
            [day1vec, timescale, isUTC] = nc_parse_date_str(timeunits);
    %         calendar = ncreadatt(ncName,'time','calendar');
            filecalendar = nc.getattvalue('time/calendar');
            day1 = datenum_cal(day1vec, filecalendar);
    %         dates = day1 + ncread(ncName,'time')/timescale;
            tstamps = nc.getvardata('time');
            if (isempty(tstamps))
                tstamps = nc.loadvar('time');
            end
            fileDates = day1 + tstamps/timescale;
            sdatevec = datevec_cal(min(fileDates), filecalendar);       % using min & max in case data is not in date order.
            edatevec = datevec_cal(max(fileDates), filecalendar);
            
                % truncate to 0000 hrs if not hourly data, so we don't lose dec. 31st by mistake.
            if (timescale <= 1)
                fileDates = floor(fileDates);
                sdatevec(4:end)=[];
                edatevec(4:end)=[];
            end

            sdnum = datenum_cal(sdatevec, filecalendar);
            ednum = datenum_cal(edatevec, filecalendar);
            ndates = round((ednum-sdnum)*timescale) + 1;
            dates = linspace(sdnum, ednum, ndates)';                     % we'll make a full time series, even if the input file has missing dates.
            missingDates = ~isequal(fileDates, dates);

            var     = nc.get(varName);
            try
                longName= var.getattvalue('long_name');
            catch
                longName = '';
            end
            try
                NAFlag  = var.getattvalue('_FillValue');
            catch
                NAFlag = [];
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

          tblNames = ["lat","lon","elev",     ];    % possible names to put in the table
        ncvarNames = ["lat","lon","elevation" ];    % equivalent names in the netcdf file.
        ntbl = length(tblNames);
        stringNames =  "merged_ncnames" ;           % list of possible char-based variables in netcdf files.
                                                    % only likely one so far is an array of filenames if the netcdf is
                                                    % a merging of either ARRM_V2 output grid files, or a merging of
                                                    % downloaded model files.

        gotlat = false;
        gotlon = false;
                % read QC table variables from nc file
        vals = cell(ntbl,1);
        varkeepers = false(ntbl,1);         % flags which variables were in the netcdf file.
        fileVars = string({nc.Variables(:).Name});
        for i=1:ntbl
            vname = ncvarNames(i);
            oname = tblNames(i);
            if (ismember(vname, fileVars))
    %             vals{i} = ncread(ncName, vname);                
                if (strcmp(oname, "lat"))
                    gotlat = true;
                    fileLats = nc.readvar(vname);
                    nlats = length(fileLats);
                elseif (strcmp(oname, "lon"))
                    gotlon = true; 
                    fileLons = nc.readvar(vname);
                    fileLons = mod(fileLons+180,360)-180;  % make sure lons are -180 to 180, not 0-360.
                    nlons = length(fileLons);
                elseif (ismember(vname, stringNames)) 
                    myvals = nc.readvar(vname);
                    varkeepers(i) = true;
                    vals{i} = strtrim(string(myvals')); 
                    varkeepers(i) = true;
                else
                    vals{i} = nc.readvar(vname);
                    varkeepers(i) = true;
                end
            end
        end

        if (~gotlat), error("error:  cannot find variable for lat in %s", ncName); end
        if (~gotlon), error("error:  cannot find variable for lon in %s", ncName); end


        nsites = nlats * nlons;

        stnID = strings(nsites,1);
        stnName = strings(nsites,1);
        index = nan(nsites,2);
        lat  = nan(nsites,1);
        lon  = nan(nsites,1);
        ix=0;
        for i=1:nlons
            for j=1:nlats
                ix=ix+1;
                index(ix,:) = [j,i];
                lat(ix)    = fileLats(j);
                lon(ix)    = fileLons(i);
            end
        end
        
        startDates = repmat(dates(1), nsites,1);
        endDates   = repmat(dates(end), nsites,1);

                % put variables into siteTbl
        ncTbl = table(stnID, stnName, index, lat, lon, startDates, endDates, 'VariableNames',{'stnID','stnName','index','lat','lon','startDate','endDate'});
        for i=1:ntbl
            if (varkeepers(i))
                vname = char(tblNames(i));
                try
                ncTbl.(vname) = vals{i};
                catch
                    fprintf(2, "warning:  problem inserting variable %s into table", vname);
                end
            end
        end

        ncTbl.Properties.UserData.timeunits = timeunits;
        ncTbl.Properties.UserData.calendar = filecalendar;
        ncTbl.Properties.UserData.day1 = day1;
        ncTbl.Properties.UserData.isUTC = isUTC;
        ncTbl.Properties.UserData.isNCtbl=true;       % flags this as a QC stn table.
        ncTbl.Properties.UserData.dates = dates;  
        
        ncTbl.Properties.UserData.filecalendar = filecalendar;
        ncTbl.Properties.UserData.fileDay1 = day1;
        ncTbl.Properties.UserData.fileDates = fileDates;    
        ncTbl.Properties.UserData.missingDates = missingDates;    
        ncTbl.Properties.UserData.fileDateRange=sprintf("%s %s", string(datestr_cal([fileDates(1); fileDates(end)],filecalendar,"yyyy-mm-dd")));
        ncTbl.Properties.UserData.dateRange    =sprintf("%s %s", string(datestr_cal([    dates(1);     dates(end)],filecalendar,"yyyy-mm-dd")));

        ncTbl.Properties.UserData.ncName = ncName;                    % capturing output.  add Userdata info to table
        ncTbl.Properties.UserData.varName = varName;
        ncTbl.Properties.UserData.npts = length(dates);
        ncTbl.Properties.UserData.NAFlag = NAFlag;
        ncTbl.Properties.UserData.longName = longName;
        ncTbl.Properties.UserData.units = varUnits;
        ncTbl.Properties.UserData.nsites = nsites;
        ncTbl.Properties.UserData.fileLats = fileLats;
        ncTbl.Properties.UserData.fileLons = fileLons;

        table_source = dir_canonical(ncName);
        ncTbl.Properties.UserData.table_source=table_source;

        ncvarNames = ncvarNames(varkeepers);
        ncTbl.Properties.UserData.nctblNames = ncvarNames;

        ncTbl.Properties.UserData.globalAttributes = ncdf_get_attributes(nc);
        try
            ncTbl.Properties.UserData.RunParams = ncdf_get_attributes(nc,"RunParams");
        catch
            ncTbl.Properties.UserData.RunParams = [];
        end
        try
            ncTbl.Properties.UserData.DownscalingParams = ncdf_get_attributes(nc,"DownscalingParams");
        catch
            ncTbl.Properties.UserData.DownscalingParams = [];
        end

    end
end  

function siteTbl = load_table_vars(siteTbl, nc, loadTblVars)
% loads any non-standard 1-D variable from the table which uses stnnum as the dimension. 
%   or any character variable with stnnum as the 2nd dimension and 1st dimension < 256. 

    if (islogical(loadTblVars) && ~loadTblVars), return; end
    
    if (isempty(nc)), nc=netcdf(siteTbl.ncName); end
    
    nctblNames = siteTbl.Properties.UserData.nctblNames;
    varName = siteTbl.Properties.UserData.varName;
    
    nvars = length(nc.Variables);
    for i=1:nvars
        v=nc.Variables(i);
        if (strcmpi(v.Name, varName) || any(strcmpi(v.Name,nctblNames))), continue; end    % we already loaded these.
        [isStnVar, isStnText] = is_site_variable(nc.Variables(i));
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

function  [isStnVar, isStnText] = is_site_variable(ncVariable)
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

function [stnID, stnName] = create_ids_and_names(ncTbl)
%   creates strings for stnID and stnName.
%   stnID   is "(lat,lon)"
%   stnName is "[lat_index,lon_index]"

    nsites  = size(ncTbl,1);
    stnID   = strings(nsites,1);
    stnName = strings(nsites,1);

    for i=1:nsites
        stnID(i)   = sprintf("(%8.4f,%9.4f)", ncTbl.lat(i), ncTbl.lon(i));
        stnName(i) = sprintf("[%4d,%4d]", ncTbl.index(i,:));
    end

end

