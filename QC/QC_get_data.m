function [tbl, nc] = QC_get_data(tbl, nc, start_date, end_date, outcal, maxmem, varnames, use_varNames, do_random, load_singly, latlon_keepers)
% QC_get_data:  returns updated site table with data from underlying netcdf file
%
%       This should be redone with an initParams(....) function, Ian!
%
%   Inputs:  
%       tbl             station table with sites to be read  
%                           from a call to QC_get_site_table(...), which reads
%                           the important metadata from the netcdf file given in tbl's UserData
%       nc              ncdf object, or name of nc file to use (if empty, will look in tbl's UserData to get ncname)
%       start_date,     datevec's of start and end dates  [yyyy, mm, dd] of data to retrieve.  
%                           'true' means start at beginning of file
%                           set to year only to start at beginning of year
%       end_date            set to true to retrieve to end-of-file.
%                           set to year only to stop and end of year.
%       outcal          string.  convert dates to calendar outcal ("day-365","julian", etc.) see function calendar_length() for valid calendar strings.
%                           for backward compatibility:  can be logical true/false for older version with variable "removeleaps".
%                               if true, sets outcal to "365-day".  If false, outcal = tbl's calendar.                         
%       maxmem          max amount of memory to use.  Will abort if too much data is attempted to be read. [ 1 GB ]
%       varnames        list of variables to load.  Data will be loaded into columns with varnames as the column names.
%                           If empty or missing, tbl.Properties.UserData.varName's  data is loaded into  column 'data'.
%                           unless use_varName is set to true
%                           Currently, varnames should be a single variable, the name of the main data variable in the nc
%                           file (Tmin, Tmax, Prec, etc.).
%       use_varNames    boolean.  If true, use the climate variables names for the columns.
%                                 If false, missing or empty, first climate variable column is named "data", and first zvals column is named "zvals"  
%       load_singly     boolean.  If true, loads each site singly.  If false, reads each site separately.  If true, loads all data (see keepers)
%       latlon_keepers  boolean or numerical array of which sites' data to keep.  If numeric, is an array of indexes where to load from full array.
%                                 if present:  tbl.(varname) = data(keepers,:);
%                                 Use keepers if keeping only some sites.  If sites to keep are in same order as data
%                                 file, keepers can be a boolean array.  If sites are not in the same order, then
%                                 keepers must be a numerical array telling where each output line should come from in
%                                 the full input data.
%                           
%   Outputs
%       tbl             updated Site Table, with data column added.
%           table fields:
%                   stnID
%                   lat
%                   lon
%                   elev
%                   stnName
%                   startDate   (adjusted to be true datenums for calendar-length -- see notes below) 
%                   endDate     (adjusted to be true datenums for calendar-length -- see notes below)
%                               (adjusted to date range if original data for site covered wider period than specified
%                                by start_date and end_date)
%                   pctValid
%                   (data)      And data for all climate variables or for climate variables passed in in varnames. 
%                                   NOTE:  *_zvals are considered climate variables.
%                                   If use_varNames is true, then data columns have original variable names 
%                                       if false, than data columns will be names "data" and "zvals".
%
%
%       and additional info in the table's Properties.UserData struct:
%           ncName  name of netcdf file where data originated
%           npts    # of time values
%           nsites  # of sites in file
%           timeunits   timeunits string from time variable
%           calendar    standard, day-365 or day-360, based on the data.
%                           calendar is adjusted to day-365 if remove_leaps && orig calendar was standard.
%           day1        datenum of timeunits' "days since" date.
%           dates   matlab datenums of dates read from file. or datenum365's or datenum360's. 
%                       note:  datenums are adjusted to 365-day calendar if leap-days removed.
%                       Note:  timestamps in netcdf file are days-since-somedate, but are returned here 
%                              as datenums (days since 0-0-0) (datenum_cal(..., calendar), days based on calendar).
%                       Note:  calendar NOT adjusted to day-365 from day-360, regardless of setting of remove_leaps.
%                              QC files should all be standard;
%                       NOTE:  If start_date or end_date are outside the range of the file, the table's data starts or
%                              ends at the file's date range, not the specified range.  If you need the full range of
%                              dates, you'll need to expand it yourself (should make this an option, Ian...)
%
% QC data files have no gaps
% fills in the tbl's data fields and UserData
%
%   ToDo:  allow nc file to be gridded;  in that case, extract the nearest lat/lon points to each station location from the nc file
%
%
%   2022-09-19  icsf    added tolerance to lats & lons when specifying sites by lat/lon so can get closest location for given lat/lon.
%_____________

    ncName = tbl.Properties.UserData.ncName;
    varName = tbl.Properties.UserData.varName;
    calendar = tbl.Properties.UserData.calendar;
    day1 = tbl.Properties.UserData.day1;
    dates = tbl.Properties.UserData.dates;
    nsites = size(tbl,1);
    datevec_range = datevec_cal(strsplit(tbl.Properties.UserData.dateRange,' '), calendar);
    if (~exist('nc','var')            || isempty(nc)),            nc = []; end
    if (~exist('start_date','var')    || isempty(start_date)),    start_date = datevec_range(1,:); end % use table's specified start & end date if not given here.
    if (~exist('end_date','var')      || isempty(end_date)),      end_date   = datevec_range(2,:); end
    if (~exist('outcal','var')        || isempty(outcal)),        outcal = strings(0);   end
    if (~exist('maxmem','var')        || isempty(maxmem)),        maxmem = 1e9;   end
    if (~exist('varnames','var')),  varnames = strings(0); else, varnames = string(varnames); end       % We'll assume all climate variables unless user specifies a list to load.
    if (~exist('use_varNames','var')  || isempty(use_varNames)),  use_varNames = false; end
    if (~exist('do_random','var')     || isempty(do_random)),     do_random = false; end
    if (~exist('load_singly','var')   || isempty(load_singly)),   load_singly = false; end
    if (~exist('latlon_keepers','var')),                          latlon_keepers = []; end
    
        % handle old version where outcal was "removeleaps".
    if (islogical(outcal) || (isnumeric(outcal) && any(outcal==[0,1])))
        if (outcal && calendar_length(calendar) == 365.25)
            outcal = "day-365";
        else
            outcal = calendar;
        end
    elseif (isempty(outcal) || strlength(outcal)==0)
        outcal = calendar;
    end
    if (calendar_length(outcal) ~= calendar_length(calendar))
        adjust_calendar = true;
    else
        adjust_calendar = false;
    end
    
    if (length(start_date)==1), start_date = [start_date,1,1]; end
    if (length(end_date)==1)
        if (calendar_length(calendar)==360)
            end_date = [end_date,12,30];
        else
            end_date = [end_date,12,31];
        end
    end
            % make specified start & end dates datenums, and limit them to the dates in the file.
    start_dnum = max(dates(1),   datenum_cal(start_date,calendar));
    end_dnum   = min(dates(end), datenum_cal(end_date,  calendar));
    
    if (isempty(nc))
        nc=ncdf(ncName);
    elseif(isstring(nc) || ischar(nc))
        nc = ncdf(nc);
    end
    
    if (isempty(varnames))
        varnames = nc.varlist;
    end
    varnames = varnames(is_climate_variable(varnames));
            
        % translate startDate, endDate columns back to dnums if they were text.
    if (~isnumeric(tbl.startDate))
        dates_were_text=true;
        dlen = max(strlength(tbl.endDate));
        if (dlen <=12)
            calfmt="yyyy-mm-dd";
        else
            calfmt="yyyy-mm-dd HH:MM";
        end
        tbl.startDate = datenum_cal(tbl.startDate,calendar);
        tbl.endDate   = datenum_cal(tbl.endDate,  calendar);
    else
        dates_were_text = false;
    end

            % if all data requested, get all data for the sites
    if (start_dnum > dates(1) || end_dnum < dates(end) || adjust_calendar)
        recalc_pcts = true;
    else
        recalc_pcts = false;
    end
    start_ix = start_dnum - dates(1)+1;
    end_ix   = end_dnum   - dates(1)+1;
    dates = dates(start_ix:end_ix);
    mpts = end_ix - start_ix + 1;
    
        % before loading the file's data variable, figure out the memory needed and make sure we can use that much memory.
        %
    varnames = to_column(string(varnames));
%   varName_included=true;
    if (~any(strcmp(varName, varnames)))
        varnames = [varnames; string(varName)]; 
%       varName_included = false;
    end
    nvars = length(varnames);
    ddtype = strings(nvars,1);
    totmem = 0;
    for v=1:nvars
        vv=nc.get(varnames(v));
        dtype = vv.Datatype;
        if (any(strcmp(dtype,["single","double"])))
            ddtype(v) = dtype;
        elseif (any(strcmp(dtype,["int8","int16","int32","uint8","uint16","uint32"])))
            ddtype(v) = "single";
        elseif (any(strcmp(dtype,["int64","uint64"])))
            ddtype(v) = "double";
        elseif (strcmp(dtype,'char'))
            error("error:  can't load char data yet");
        else
            error("error:  can't load %s data yet", dtype);
        end
        
        if (strcmp(ddtype(v),"single"))
            dbytes=4;
        else
            dbytes=8;
        end        
        totmem = totmem + nsites * mpts * dbytes;
    end
    
    if (totmem > maxmem)
        error("too many stations to load data (%.1f Gbytes total); specify larger maxmem (%g) or reduce number of stations", totmem/1e9,  maxmem); 
    end
    
%         %   Load the file's varName's data into column called "data".
%         %
%     if (~varName_included)
%         tbl.data = nan(nsites,mpts, ddtype(end));   % we stored varNames type at the end of ddtype(...)
%             
%         for i=1:nsites
%             if (~isnan(tbl.index(i)))  % index is nan if no data for site in file
%                 tbl.data(i,:) = ncread(ncName,varName,[start_ix,tbl.index(i)],[mpts,1],[1,1]);
%             end
%         end
%     end

        % load the list of variables wanted.
        
    [all_lats,all_lons] = ncdf_get_latlons(nc);
    latlons = [tbl.lat,tbl.lon];
    nlatlons = size(latlons,1);
    if (isempty(latlon_keepers))
        tol = .5;   % accept closest point within .5 degrees, about 55 km
        [latlon_keepers, llkeepers_in_tolerance] = closest_latlon(latlons, all_lats, all_lons, tol);
        if (sum(llkeepers_in_tolerance) ~= nlatlons), error("error:  cannot match latlons to lats & lons in file with tolerance of %.3 degrees", tol); end
            % if we're keeping everything, clear out latlon_keepers.
        if (all(to_column(latlon_keepers) == to_column(1:nlatlons)) && nlatlons == numel(all_lats))
            latlon_keepers = [];
        end
    end 
    
        % check to see if we're keeping all the points, and in their original order.  If so, no need to use
        % latlon_keepers.
    if (~isempty(latlon_keepers))
        if (islogical(latlon_keepers))
            if (numel(latlon_keepers) == numel(all_lats) && sum(latlon_keepers) == length(latlon_keepers))
                latlon_keepers = [];
            end
        else
            if (numel(latlon_keepers) == numel(all_lats) && all(to_column(latlon_keepers) == to_column(1:nlatlons)))
                latlon_keepers = [];
            end
        end
    end

    for i=1:length(varnames)
        vname = varnames(i); %char(varnames(v));
        if (load_singly)
            tbl.(vname) = nan(nsites, mpts, ddtype(i));
            for v=1:nsites
                if (~isnan(tbl.index(v)))  % index is nan if no data for site in file
                    tbl.(vname)(v,:) = ncread(ncName,vname,[start_ix,tbl.index(v)],[mpts,1],[1,1]);
                end
            end
        else
            data = ncread(ncName, vname, [start_ix,1],[mpts,inf],[1,1]);
            if (~isempty(latlon_keepers))
                tbl.(vname) = data(:, latlon_keepers)';
            else
                tbl.(vname) = data';
            end
            clear data;
        end
            
            % if user didn't want the data called by the varName, then rename the data column to data.
        if (~use_varNames)
            if (strcmp(vname,varName))
                vix = find(strcmp(tbl.Properties.VariableNames, vname),1);
                if (~isempty(vix))
                    tbl.Properties.VariableNames{vix} = 'data';
                end
%               tbl.data = tbl.(vname);
            elseif (strcmp(vname, sprintf("%s_zvals", varName)))
                vix = find(strcmp(tbl.Properties.VariableNames, vname),1);
                if (~isempty(vix))
                    tbl.Properties.VariableNames{vix} = 'zvals';
                end
            end
        end
    end

        % Ian:  make sure dates are datenums!  
        % and make sure varnames are all data variables, including zvals.
    if (adjust_calendar)
        [out_dnums, inserted, deleted, to_ix, keepers] = calendar_convert(dates, calendar, outcal, "do_random",do_random);
        for v=1:length(varnames)
            if (strcmp(varnames(v),varName) && ~use_varNames)
                vname = 'data';
            else
                vname = varnames(v);
            end
            if (~isempty(deleted))
                tbl.(vname) = tbl.(vname)(:,keepers);
            elseif (~isempty(inserted))
                tbl.(vname)(:,to_ix) = tbl.(vname);
                tbl.(vname)(:,inserted) = nan;
            end
        end
%         for i=1:nsites
%             if (~isempty(tbl.data{i}))
%                 tbl.data{i} = tbl.data{i}(keepers); 
%             end
%         end
            % update all the calendar variables to reflect the change in calendar
        calendar = outcal;
        dates = out_dnums;
        start_dnum = dates(1);
        end_dnum = dates(end);
%         day1_vec = datevec(day1);
%         day1 = datenum_cal(day1_vec, calendar); 
%         dvecs = datevec_cal(dates);
%         dates = datenum365(dvecs);
%         start_dnum = datenum365(datevec(start_dnum));
%         end_dnum   = datenum365(datevec(end_dnum));
        recalc_pcts = true;
    end
    
        % this should exclude zvals, Ian!
        % recalculate pctValid if # of data points has changed
    if (recalc_pcts)
            % get new date Ranges for each station
        for j=1:nsites
            vnames = varnames(is_climate_variable(varnames,true));          % get list of climate variables (true -> skip zvals)
            if (isempty(vnames))
                vnames = varnames(is_climate_variable(varnames,false));    % if empty, then include zvals
            end
            if (~isempty(vnames)) 
                vname1 = vnames(1);                                         % and recalculate pct_valid for 1st climate variable we find.
                if (strcmp(vname1,varName) && ~use_varNames)                % (Ian:  this should be smarter...look for varname first?
                    vname = 'data';
                else
                    vname = vname1;
                end
                good = ~isnan(tbl.(vname)(j,:));
                ix1=find(good,1);
                if (isempty(ix1))
                    tbl.startDate(j)=0;
                    tbl.endDate(j)=0;
                    tbl.pctValid(j)=0;
                else
                    ix2=find(good,1,'last');
                        tbl.startDate(j)= dates(ix1);
                        tbl.endDate(j)  = dates(ix2);
                    ndays               = tbl.endDate(j) - tbl.startDate(j) + 1;
                    ngood = sum(good);
                    if (ndays == 0)
                        tbl.pctValid(j)=0;
                    else
                        tbl.pctValid(j) = 100.0*ngood/ndays;
                    end
                end
            end
        end
        keepers = tbl.pctValid > 0;
        if (sum(keepers)==0), error("error: no stations with valid data in criteria specified"); end
        tbl = tbl(keepers,:);
    end    
    
            % update UserData in table Properties.
            
    nsites = size(tbl,1);
    tbl.Properties.UserData.nsites = nsites;
    tbl.Properties.UserData.modifiedDates = true;
    tbl.Properties.UserData.npts = mpts;
    tbl.Properties.UserData.day1 = day1;
    tbl.Properties.UserData.dates = dates;
    tbl.Properties.UserData.calendar = outcal;
    tbl.Properties.UserData.dateRangeLimit = sprintf("%s %s", string(datestr(datevec_cal([start_dnum; end_dnum],calendar),'yyyy-mm-dd')));
    tbl.Properties.UserData.originalCalendar = calendar;    
        
    if (dates_were_text)
        tbl.startDate = datestr_cal(tbl.startDate,calendar, calfmt);
        tbl.endDate   = datestr_cal(tbl.endDate,calendar,   calfmt);
    end
end

