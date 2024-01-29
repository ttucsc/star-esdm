function tbl = gridded_get_data(tbl, varargin)
% function tbl = gridded_get_data(tbl, start_date, end_date, varargin)
% gridded_get_data:  returns updated site table with data from underlying netcdf file
%
%  tbl = gridded_get_data(tbl, start_date, end_date, "latrange", latrange, "lonrange", lonrange, ...
%                         "calendar", outcalendar, "varname", varname, "ncdf", nc, "maxmem", maxmem)
%
%   Inputs:  
%       tbl                 station table with sites to be read  
%                               from a call to QC_get_site_table(...), which reads
%                               the important metadata from the netcdf file given in tbl's UserData
%       start_date,         datevec's of start and end dates  [yyyy, mm, dd] of data to retrieve.  
%                               set to empty to start at beginning of file
%                               set to year only to start at beginning of year
%       end_date                set to empty to retrieve to end-of-file.
%                               set to year only to stop and end of year.
%     optional name/value arguments
%       "latrange", [lat1,lat2] keep only locations in half-open latitude  range [lat1,lat2);
%       "lonrange", [lon1,lon2] keep only locations in half-open longitude range [lon1,lon2);
%       'calendar', calendar    string.  If not empty, change calendar length to new calendar.
%       'do_random', do_random  true/false.  If true  and out-calendar ~= calendar, then inserts or deletes days randomly
%                                            If false and out-calendar ~= calendar, then inserts or deletes days on the
%                                            same day for every year.
%       'varNames', varNames    list of variables to load.  Data will be loaded into columns with varNames as the column names.
%                                   If empty or missing, table-varname's data is loaded into  column 'data'.
%                                   otherwise, if table's varname appears in varNames, then column name is the varname. 
%       'ncdf', nc              ncdf object from which tbl was read;  if present, avoids time needed to read netcdf metadata.
%       'maxmem', maxmem        max amount of memory to use.  Will abort if too much data is attempted to be read. [ 1 GB ]
%                           
%   Outputs
%       tbl             updated Site Table, with data column added.
%           table fields:
%                   stnID       "(lat, lon)"     string with lat/lon indexes for each site (from input)
%                   stnName     "[latix, lonix]" string with lat/lon values for each site (from input)
%                   lat         latitude & longitude for each site as numeric
%                   lon
%                   elev        (may not be present)
%                   startDate   (adjusted to correct datenums for output calendar-length -- see notes below) 
%                   endDate     (adjusted to correct datenums for output calendar-length -- see notes below)
%                               (adjusted to date range if original data for site covered wider period than specified
%                                by start_date and end_date)
%                   pctValid    % of valid data points for each site.
%
%       and additional info in the table's Properties.UserData struct:   THIS NEEDS UPDATING, IAN!
%           ncName  name of netcdf file where data originated
%           npts    # of time values
%           nsites  # of sites in file
%           timeunits   timeunits string from time variable
%           calendar    standard, day-365 or day-360, based on the data (or calendar input parameter).
%                           calendar is adjusted to new calendar if calendar input parameter is provided.
%           day1        datenum of 1st date in file
%           dates   matlab datenums of dates read from file. or datenum365's or datenum360's. 
%                       Note:  timestamps in netcdf file are  days-since-day1, but are returned here 
%                              as datenums (days since 0-0-0).
%           fileDay1
%           fileDates
%           fileCalendar
%           missingDates
%           modifiedDates
%           and some others...fill this in when it's all working right, Ian.
%
% ncTbl data files have no gaps - this function sets any missing data to NAs.
%           

    [tbl, start_date, end_date, latrange, lonrange, outcalendar, varNames, nc, do_random, maxmem] = initParams(tbl, varargin{:});

    ncName        = tbl.Properties.UserData.ncName;
    fileVarName   = tbl.Properties.UserData.varName;
    calendar      = tbl.Properties.UserData.calendar;
%   day1          = tbl.Properties.UserData.day1;
%   dates         = tbl.Properties.UserData.dates;
%   fileDay1      = tbl.Properties.UserData.fileDay1;
    fileDates     = tbl.Properties.UserData.fileDates;
%   missingDates  = tbl.Properties.UserData.missingDates;
%   modifiedDates = tbl.Properties.UserData.modifiedDates;
    allSites      = tbl.Properties.UserData.allSites;
    lats          = tbl.Properties.UserData.lats;
    lons          = tbl.Properties.UserData.lons;
    fileLats      = tbl.Properties.UserData.fileLats;
    fileLons      = tbl.Properties.UserData.fileLons;
    
    nsites = size(tbl,1);
    
            % make specified start & end dates datenums, and limit them to the dates in the file.
    start_dnum = datenum_cal(start_date,calendar);
    end_dnum   = datenum_cal(end_date,  calendar);
    
    if (start_dnum ~= fileDates(1) || end_dnum ~= fileDates(end))
        modifiedDates = true;
    else
        modifiedDates = false;
    end

    if (~isempty(latrange) || ~isempty(lonrange))
        if (~isempty(latrange))
            latkeepers = tbl.lat >=latrange(1) & tbl.lat <= latrange(2);
            if (sum(latkeepers)==0), error("error:  no gridcells in lat range %.4f - %.4f)", latrange); end
            allSites = allSites & all(latkeepers);
            lats = lats(lats >= latrange(1) & lats <= latrange(2));
        end        

        if (~isempty(lonrange))
            lonkeepers = tbl.lon >=lonrange(1) & tbl.lon <= lonrange(2);
            if (sum(lonkeepers)==0), error("error:  no sites in lon range %.4f - %.4f)", lonrange); end 
            allSites = allSites & all(lonkeepers);
            lons = lons(lons >= lonrange(1) & lons <= lonrange(2));
        end
        keepers = latrange & lonrange;
        tbl = tbl(keepers,:);
        nsites = size(tbl,1);
    end    
    
    tbl = sortrows(tbl,{'lon','lat'});
    
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

    date_start_ix = find(fileDates >= start_dnum,1);
    date_end_ix   = find(fileDates <= end_dnum, 1,'last');
    dates = (start_dnum:end_dnum)';
    mpts = length(dates);
    file_extract_dates = fileDates(date_start_ix:date_end_ix);
    mpts2 = length(file_extract_dates);
    missingDates = ~isequal(dates, fileDates(date_start_ix:date_end_ix));
    
    if (missingDates)
        out_ix = file_extract_dates - dates(1)+1;
    end
    
    lat_ix1  = find(fileLats == min(tbl.lat), 1);
    lat_ix2  = find(fileLats == max(tbl.lat), 1);
    lon_ix1  = find(fileLons == min(tbl.lon), 1);
    lon_ix2  = find(fileLons == max(tbl.lon), 1);
    
    nlats = lat_ix2 - lat_ix1 + 1;
    nlons = lon_ix2 - lon_ix1 + 1;
    
    read_start = [date_start_ix, lat_ix1, lon_ix1];
    read_count = [mpts2,         nlats,   nlons];
    
    totmem = calc_mem_usage(varNames, nc, nsites, mpts);
    if (totmem > maxmem)
        error("too many stations to load data; specify larger maxmem (%g) or reduce number of stations", maxmem); 
    end
   
        % read data and put into column in table.
    tbl.pctValid = nan(nsites,1); 
    nvars = length(varNames);
    for i=1:nvars
        mydata = ncread(ncName, varNames(i), read_start, read_count);
        mydata = (reshape(mydata,mpts2,nlats*nlons))';
        if (missingDates)    % if file has missing data, we need to re-position all the data and leave nan's where the data is missing.
            md=mydata;
            mydata=nan(nsites,mpts);
            mydata(:,out_ix) = md;
        end
        if (strcmp(varNames(i),fileVarName))
            tbl.data = mydata;
        end
        if (nvars > 1 || ~strcmp(varNames(i),fileVarName))
            tbl.(varNames(i)) = mydata;
        end
    end
    
    if (length(varNames)==1 && strcmp(varNames, fileVarName))
        allvars = "data";
    else
        allvars = [to_column(varNames); "data"];
    end
        
            % update UserData in table Properties.
            
    tbl.Properties.UserData.nsites = nsites;
    tbl.Properties.UserData.modifiedDates = modifiedDates;
    tbl.Properties.UserData.npts = mpts;
    tbl.Properties.UserData.dates = dates;
    tbl.Properties.UserData.dateRange = sprintf("%s %s", string(datestr(datevec_cal([start_dnum; end_dnum],outcalendar),'yyyy-mm-dd')));
    tbl.Properties.UserData.allSites = allSites;
    tbl.Properties.UserData.lats = lats;
    tbl.Properties.UserData.lons = lons;
        
    if (calendar_length(outcalendar) ~= calendar_length(calendar))        
        tbl = adjust_calendar(tbl, outcalendar, do_random, allvars);
        dates = tbl.Properties.UserData.dates;
    end
    
        % calculate pctValid for the main data variable

    if (any(ismember('data',tbl.Properties.VariableNames)))
                % calculate % valid for main table variable
        for i=1:nsites
            nnans = sum(~isnan(mydata(i,:)));
            tbl.pctValid(i) = 100*(nnans/size(mydata,2));
        end
    end
    
    tbl.startDate = repmat(dates(  1), nsites,1);
    tbl.endDate   = repmat(dates(end), nsites,1);
    if (~any(tbl.pctValid < 100.0))
        for i=1:nsites
            if (tbl.pctValid(i) < 100)
                ix1 = find(~isnan(tbl.data(i,:)),1);
                ix2 = find(~isnan(tbl.data(i,:)),1,'last');
                tbl.startDate(i) = dates(ix1);
                tbl.endDate(i)   = dates(ix2);
            end
        end
    end
                    

    if (dates_were_text)
        tbl.startDate = datestr_cal(tbl.startDate,calendar, calfmt);
        tbl.endDate   = datestr_cal(tbl.endDate,calendar,   calfmt);
    end
end

function [tbl, start_date, end_date, latrange, lonrange, outcalendar, varNames, nc, do_random, maxmem] = initParams(tbl, varargin)

    p = inputParser;
    p.StructExpand = true;
    p.KeepUnmatched = false;     % throw error if any unexpected input parameters.
    
    addOptional(p,"start_date",   [],       @(s) isempty(s) || isnumeric(s) && isrow(s) && any(length(s)==[1,3,6]));
    addOptional(p,"end_date",     [],       @(s) isempty(s) || isnumeric(s) && isrow(s) && any(length(s)==[1,3,6]));
    addParameter(p,"latrange",    [],       @(s) isempty(s) || (isnumeric(s) && numel(s)==2 && s(1)>= -90 && s(1)<  90 && s(2)>s(1) && s(2)<=90));
    addParameter(p,"lonrange",    [],       @(s) isempty(s) || (isnumeric(s) && numel(s)==2 && s(1)>=-180 && s(1)< 180 && s(2)>s(1) && s(2)<=180));
    addParameter(p,"varNames",    "",       @(s) ischar_s(s));
    addParameter(p,"calendar",    "",       @(s) ischar_s(s));
    addParameter(p,"ncdf",        [],       @(s) isempty(s) || strcmp(class(s), "ncdf"));
    addParameter(p,"do_random",   true,     @(s) islogical(s) || (isnumeric(s) && s>=1));
    addParameter(p,"maxmem",      1e9,      @(s) isnumeric(s) && s>=10e6);
    parse(p, varargin{:});
    Parms = p.Results;
    
    start_date   = Parms.start_date;
    end_date     = Parms.end_date;
    latrange     = Parms.latrange;
    lonrange     = Parms.lonrange;
    varNames     = Parms.varNames;
    outcalendar  = Parms.calendar;
    do_random    = Parms.do_random;
    nc           = Parms.ncdf;
    maxmem       = Parms.maxmem;
    
    calendar = tbl.Properties.UserData.calendar;
    dates    = tbl.Properties.UserData.dates;
    ncName   = tbl.Properties.UserData.ncName;
%   fileDates= tbl.Properties.UserData.fileDates;
    
    if (isempty_s(outcalendar)), outcalendar = calendar;                          end
    if (isempty(start_date)),    start_date  = datevec_cal(dates(1),   calendar); end
    if (isempty(end_date)),      end_date    = datevec_cal(dates(end), calendar); end
    
%         % limit start, end dates to date range of tbl or file, whichever is less..        SKIP FOR NOW...
%     start_date = datevec_calendar(max([datenum_calendar(start_date), dates(1),   fileDates(1)]),   calendar);
%       end_date = datevec_calendar(min([datenum_calendar(  end_date), dates(end), fileDates(end)]), calendar);
    
    varNames = string(varNames);  
    if (isempty_s(varNames))
        varNames = tbl.Properties.UserData.varName;
    end
        
    if (length(start_date)==1), start_date = [start_date,1,1]; end
    if (length(end_date)==1)
        if (calendar_length(calendar)==360)
            end_date = [end_date,12,30];
        else
            end_date = [end_date,12,31];
        end
    end
    

    if (isempty(nc)), nc = ncdf(ncName, "create",false); end

end

function ncTbl = adjust_calendar(ncTbl, outcalendar, do_random, varNames)

        % make sure data is in the list of varNames to adjust.
    varNames = string(varNames);
    if (~ismember("data",varNames)), varNames = [to_column(varNames); "data"]; end
        
        % switch random number generators so we use the same one for serial or parallel processing.
    if (isnumeric(do_random))
        orig_rng = rng(do_random, 'multFibonacci');
    else
        orig_rng = rng(2,'multFibonacci');
    end
            
    outlen = calendar_length(outcalendar);
    calendar = ncTbl.Properties.UserData.calendar;
    inlen  = calendar_length(calendar);
    if (outlen == inlen), return; end
    dates = ncTbl.Properties.UserData.dates;
    dvecs = datevec_cal(dates, calendar);
    ndates = length(dates);
    yrs = unique(dvecs(:,1));
    nyrs = length(yrs);
        % figure out where we want to insert or delete days
    changers = false(ndates,1);       % list of dates we want to change.
        % places to insert for 360-day starting calendar or delete for 360-day ending calendar
    if (outlen == 360 || inlen == 360)
        if (~do_random)
            dns = repmat(datenum_cal([yrs,ones(nyrs,1),zeros(nyrs,1)],calendar),5,1) + to_column(repmat([36,109,182,255,328],nyrs,1));
            dns = sort(dns);
            dns = dns(dns >= dates(1) & dns <= dates(end));
            ix  = dns - dates(1)+1;
            changers(ix) = true;
        else
            change_days = (0:4)*73;     % start of each 1/5 of the year.
            yrstarts = datenum_cal([yrs,ones(nyrs,1),ones(nyrs,1)], calendar) - dates(1);    % datenums for start of each year (minus 1)
            for i=1:nyrs
                for j=1:5
                    kx = yrstarts(i) + change_days(j) +randi(73);   % choose a random day in the jth part of the year.
                    changers(kx) = true;
                end
            end
        end
    end
    if (inlen == 365.25 || outlen == 365.25)     % and flag leap days if needed
        if (inlen == 365.25)
            leaps = find(dvecs(:,2)==2 & dvecs(:,3) == 29);
        else
            leapyrs = find(leapyear(yrs));
            leaps = [];
            for i=1:length(leapyrs)
                ix = find(dvecs(:,1)==yrs(leapyrs(i)) & dvecs(:,2)==2 & dvecs(:,3)==28, 1);
                if (~isempty(ix)), leaps(end+1) = ix; end %#ok<AGROW>
            end
        end
        if (do_random)
            for i=1:length(leaps)
                if (changers(leaps(i))==true)      % if we randomly chose a leap day to change, then chose a different random day.
                    ix = leaps(i)-60 + randi(73);
                    while (changers(ix)==true), ix = leaps(i)-60 + randi(73); end   % in case we chose feb. 29th again!
                    leaps(i) = ix;  % flag this one, since the leap day had already been randomly flagged.
                end
            end
        end
        changers(leaps) = true;        
    end
    
            % if first (last) point is being changed, switch it with the next (previous) so we  keep the full date range.
    if (changers(1)==true)
        ix=find(~changers,1);
        changers(1)=false;
        changers(ix)=true;
    end
    if (changers(end)==true)
        ix=find(~changers,1,'last');
        changers(end)=false;
        changers(ix)=true;
    end
    
            % now set up the arrays of what to keep or where to move things.
    if (outlen < inlen)
        out_ix = [];
        keepers = true(1,ndates);       % flag the ones we want to keep
        keepers(changers) = false;
        if (outlen == 360 && dvecs(end,2)==12 && dvecs(end,3)==31)    % change last day to Dec 30th if 360-day calendar and ending at years' end.
           dvecs(end,3)=30;
        end
    else
        keepers = [];
       if (inlen == 360 && dvecs(end,2)==12 && dvecs(end,3)==30)
           dvecs(end,3)=31; % include last december 31'st if changing from 360-day calendar            
       end
        out_ix = 1:length(dates);
        changers = find(changers);  % changers now is list of where to insert extra days.
        for i=1:length(changers)
            out_ix(changers(i):end) = out_ix(changers(i):end)+1;    % indexes of where to move them to.
        end
    end
    
    dn1 = datenum_cal(dvecs(  1,:), outcalendar);
    dn2 = datenum_cal(dvecs(end,:), outcalendar);    
    outDates = dn1:dn2;   % new dates
    ndates_out = length(outDates);
    %   make sure we have the right number of data points left:
    
    if (~isempty(keepers))
        if (sum(keepers) ~= ndates_out)
            error("problem:  nkeepers (%d) doesn't match new output size (%d)", sum(keepers), ndates_out);
        end
    else
        if (out_ix(end) ~= ndates_out)
            error("problem:  # of output days (%d) doesn't match new output size (%d)", out_ix(end), ndates_out);
        end
    end
    
        % now shrink the ncTbl's variables to just the keepers, or expand them by moving them to out_ix.
    nsites = size(ncTbl,1);
    for i=1:length(varNames)
        vname = varNames(i);
        if (~ismember(vname, ncTbl.Properties.VariableNames)), continue; end    % will skip over main varName if it is only in the file as "data";
        if (~isempty(keepers))
            ncTbl.(vname) = ncTbl.(vname)(:,keepers);
        else
            dtype = class(ncTbl.(vname));
            out_data = nan(nsites, ndates_out, dtype);
            out_data(:,out_ix) = ncTbl.(vname);
            ncTbl.(vname) = out_data;
        end
    end
       
            % update the UserData 
    ncTbl.Properties.UserData.modifiedDates = true;
    ncTbl.Properties.UserData.dates = outDates;   % new dates
    ncTbl.Properties.UserData.npts = ndates_out;
    ncTbl.Properties.UserData.calendar = outcalendar;
    ncTbl.Properties.UserData.dateRange = sprintf("%s %s", string(datestr(datevec_cal([outDates(1); outDates(end)],outcalendar),'yyyy-mm-dd')));        

        % and put the random number generator back where it was.
    rng(orig_rng);
end

function totmem = calc_mem_usage(varNames, nc, nsites, mpts)

    nvars = length(varNames);
    totmem = 0;
    for i=1:nvars
        try
            vv=nc.get(varNames(i));
        catch
            error("error:  %s:  cannot locate variable %s in %s", mfilename, varNames(i), nc.Filename); 
        end
        
        dtype = vv.Datatype;
        if (any(strcmp(dtype,["single","double"])))
            ddtype = dtype;
        elseif (any(strcmp(dtype,["int8","int16","int32","uint8","uint16","uint32"])))
            ddtype = "single";
        elseif (any(strcmp(dtype,["int64","uint64"])))
            ddtype = "double";
        elseif (strcmp(dtype,'char'))
            error("error:  can't load char data yet");
        else
            error("error:  can't load %s data yet", dtype);
        end
        
        if (strcmp(ddtype,"single"))
            dbytes=4;
        else
            dbytes=8;
        end
        
        totmem = totmem + nsites * mpts * dbytes;
    end

end