function [vals, dates, min_date, max_date, nvals, lat, lon, elev, calendar, timescale, isUTC] = read_QC_data_nc(stnID, start_date, end_date, ncName, varName, remove_leaps)
% reads QC data for a single station stnID from QC netcdf file
% returns dates, values, min & max dates, & # of valid datapoints, & lat, lon & elevation of station.
%   all missing values are filled with NAs, so you are guaranteed to get entire
%   date range requested with no gaps.
%
%   Inputs:
%       stnID           station ID to extract from nc file fname
%                           can be stnID, or index of station within file.
%       start_date      1st date to extract.  To start at beginning of year, use just year, 1950, or full date [1950,1,1].   
%                           If empty or logical-true, starts at beginning of data.
%       end_date        last date or year to extract.  To end at end of year, use just year, 1950, or full date [1950,12,31].  (ok to use 12/31 for 360-day year.)
%                           set start_date and end_date to empty or true to read all data.
%                           If empty or logical-true, starts at beginning of data.
%                       NOTE:  array of data returned spans entire start_date to end_date, even if outside the bounds of
%                              the file.  data will be padded with NAs. 
%                               IAN:  perhaps start_date & end_date could also possibly be keepers and indexes, so could
%                                avoid rereading time units, calculating keepers and indexes, etc. every time
%                                when reading large array of sites?
%                       
%       ncName          name of netcdf file to read
%       varName         name of variable to read
%       remove_leaps    true/false.  if true, removes leap days from data and from dates.
%
%   Outputs as labeled, with following notes:
%       dates       are matlab datenums (days since Jan 0, year 0), in calendar-length specified in the time_units attribute int nc file.
%                   NOTE:  dates are in units of days, even if units given as "hours since"... or "minutes since"...
%       min_date    date of 1st valid (non-nan) data point for station
%       max_date    date of last valid data point for station
%                       mindate, maxdate are returned as strings, yyyy-mm-dd so are not dependent on calendar.

            % read the data from the file.
    if (~exist('remove_leaps','var') || isempty_s(remove_leaps)), remove_leaps = false; end
    
    if (isempty(start_date)), start_date = true; end
    if (isempty(end_date)),   end_date   = true; end
    
    ncName = char(ncName);
    varName = char(varName);
    if (isnumeric(stnID))
        stnix = stnID;
    else
        stnID = string(stnID);
        z=ncread(ncName,'stnID')';
        stnIDs=strtrim(string(z));
        stnix = find(stnIDs == stnID,1);
        if (isempty(stnix)), error('readQC_data_nc:  error:  stnID %s doesn''t exist in file %s', stnID, ncName); end
    end
    try
        calendar  = ncreadatt(ncName, 'time','calendar');
    catch
        calendar  = 'standard';  % if no calendar in file, we'll assume 365.25-day.  
    end
    timeunits = ncreadatt(ncName, 'time','units');      % should be something like "days since 1850-01-01"

    [from_vec, timescale, isUTC] = nc_parse_date_str(timeunits);
    if (timescale < 1), error('error:  invalid timeunits: %s', timeunits); end
    
    day1 = datenum_cal(from_vec, calendar);
    
%     file_dnums = day1 + ncread(ncName,'time')/timescale;
    tstamps = nc.getvardata('time');
    if (isempty(tstamps))
        tstamps = nc.loadvar('time');
    end
    file_dnums = day1 + tstamps/timescale;
    
    if (islogical(start_date)), start_date = datevec_cal(file_dnums(1), calendar); end
    if (islogical(end_date)),   end_date   = datevec_cal(file_dnums(end), calendar); end
    
    filepts = length(file_dnums);
    v=ncread(ncName,varName,[1,stnix],[filepts,1],[1,1]);       % this reads all varName data for the site

    if (start_date <= 1)        % start, end dates not given (true, false or 0).  Get all dates.
        dates = file_dnums;
        vals = v;
    else
        if (length(start_date) == 1), start_date = [start_date, 1, 1]; end
        if (length(  end_date) == 1),   end_date = [  end_date,12,31]; end
        if (length(start_date) == 3), start_date = [start_date,  0,  0,  0]; end
        if (length(end_date)   == 3), end_date   = [end_date,   23, 59, 59.1]; end

        if (calendar_length(calendar) == 360), end_date(3) = 30; end
    
        start_dnum = datenum_cal(start_date, calendar);
        end_dnum   = datenum_cal(end_date,   calendar); 
                
            % get the time step size from the file datenums.  
        dnum_difs = file_dnums(2:end) - file_dnums(1:(end-1));     % find difference between each time step.
        stepsize = min(dnum_difs);                              % find smallest one.  Assumes time is monotonically increasing.
        stepsPerDay = round(1/stepsize);    % problem here if step size doesn't divide evenly into 24 hours.  We'll settle for within 10 msec, which should be OK even if time is float rather than double..
        if (abs(stepsPerDay * stepsize - 1.0) > 1e-7), error('error:  %s:  file datenums not based on even # of steps per day', ncName); end
        days_since = make_tstamps([start_date; end_date], calendar, timeunits, stepsPerDay);
        dates = days_since + day1;
        npts = length(dates);
        vals = nan(npts,1);
%         ix = round((dates - dates(1))*stepsPerDay) + 1;     
                                            % indexes of each date relative to day 0 (1 day before 1st day we want to keep)
                                            % shouldn't be any gaps in the data, but this will put data into
                                            % correct locations if there are.  Will be a slower than just copying
                                            % range of values.
        keepers = (file_dnums >= start_dnum & file_dnums <= end_dnum);  % boolean array, true if date is within desired range.
        keep_dnums = file_dnums(keepers);
        
        keep_ix = round((keep_dnums - start_dnum) * stepsPerDay) + 1;   % location of where each data point belongs in output array                  
        vals(keep_ix) = v(keepers);              % put the keeper data points where they belong.

    end
    
    if (remove_leaps && calendar_length(calendar) == 365.25)
        dvecs = datevec(dates);
        leaps = isleapday(dates);
        dvecs(leaps,:) = [];
        vals(leaps) = [];
        dates = datenum365(dvecs);
        calendar = '365_day';
    end
        
    
    if (nargout <= 2), return; end
    
    good = ~isnan(vals);
    nvals = sum(good);

    if (nvals == 0)
        min_date=[];
        max_date=[];
    else
        minix = find(good,1);
        maxix = find(good,1,'last');
        if (calendar_length(calendar) == 360)
            mindvec = datevec360(dates(minix));
            maxdvec = datevec360(dates(maxix));
        elseif (calendar_length(calendar) == 365)
            mindvec = datevec365(dates(minix));
            maxdvec = datevec365(dates(maxix));
        elseif (calendar_length(calendar) == 365.25)
            mindvec = datevec(dates(minix));
            maxdvec = datevec(dates(maxix));
        else
            throw(MException('ARRMV2:BAD_CALENDAR',sprintf('error:  calendar not 360-day, standard/julian/gregorian or 365_day/noleap: %s', calendar)));
        end
        mindnum = dates(minix);
        maxdnum = dates(maxix);
        if (mod(mindnum,1)~= 0 || mod(maxdnum,1)~= 0)
            min_date = datestr(mindvec,'yyyy-mm-dd HH:MM:SS'); 
            max_date = datestr(maxdvec,'yyyy-mm-dd HH:MM:SS'); 
        else
            min_date = datestr(mindvec,'yyyy-mm-dd'); 
            max_date = datestr(maxdvec,'yyyy-mm-dd'); 
        end
    end
    
    % get lat, lon & elevation if requested.
    if (nargout > 5)
        lats=ncread(ncName,'lat');
        lat=lats(stnix);
    else
        lat=[];
    end
    if (nargout > 6)
        lons=ncread(ncName,'lon');
        lon=lons(stnix);
    else
        lon=[];
    end
    if (nargout > 7)
        elevs=ncread(ncName,'elevation');
        elev=elevs(stnix);
    else
        elev=[];
    end
        
end

