function [days_since, time_units, start_vec] = make_tstamps(dateRange, caltype, tunits, stepsPerDay)
% generates time_units and days_since from yr_range, t_units and calendar type.
% tunits is optional.  If present, parses tunits to get starting date, and sets day_since 
% used for files which do not have valid time info.

    if (nargin < 4), stepsPerDay=1; end
    
    calendar_len = calendar_length(caltype);
    if (length(dateRange)==2)
        dateRange=[dateRange(1),1,1; dateRange(2),12,31];
        if (calendar_len == 360), dateRange(2,3) = 30; end
    end
    
    if (~exist('tunits','var') || isempty_s(tunits))
        time_units = strtrim(sprintf('days since %4d-%02d-%02d %02d:%02d:%02d', dateRange(1,:)));
        start_vec = dateRange(1,:);
    else
        time_units = tunits;
        [start_vec, timescale] = nc_parse_date_str(tunits);
        if (timescale > 1 && size(dateRange,2) == 3), dateRange(end,4:6)=[23,59,59]; end            
    end

    timescale = timescale * stepsPerDay;
    
    if (calendar_length(caltype) == 360)
        dnum1     = datenum360(start_vec);
        dnumfirst = datenum360(dateRange(1,:));
        dnumlast  = datenum360(dateRange(end,:));
    elseif (calendar_length(caltype) == 365)
        dnum1     = datenum365(start_vec);
        dnumfirst = datenum365(dateRange(1,:));
        dnumlast  = datenum365(dateRange(end,:));
    elseif (calendar_length(caltype) == 365.25)
        dnum1     = datenum(start_vec);
        dnumfirst = datenum(dateRange(1,:));
        dnumlast  = datenum(dateRange(end,:));
    else
        throw(MException('ARRMV2:BAD_CALENDAR',sprintf('error:  calendar not standard/julian/gregorian or 365_day/noleap: %s', caltype)));
    end
    nyrs = dateRange(end,1)-dateRange(1,1) + 1;
    dnums = dnumfirst + (0:(366*nyrs*timescale))/timescale;  % more than we need...
    dnums(dnums>dnumlast)=[];                                % get rid of the ones we don't need
    days_since = dnums-dnum1;                                % and change to offset from start_vec;
    days_since = days_since';                                % and make it a column vector.
end

