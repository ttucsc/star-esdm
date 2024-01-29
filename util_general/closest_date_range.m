function [keepers_ix, prenans, postnans] = closest_date_range(date_range, calendar, start_vec, days_since)
% returns time indexes of data for date_range, given calendar, start_vec and timevals
%
%   Inputs:
%       date_range  either [start_year, end_year] or [start_yyyy_mm_dd; end_yyyy_mm_dd] of desired range of data
%       calendar    calendar string of time units
%       start_vec   datevec of day 0 for timevals
%       timevals    vector of days-since-start_vec
%       
%   Outputs:
%       keepers_ix  indexes of times within the desired date range.  (indexes of keepers).
%       prenans     # of days in date_range which precede earliest date in inputs
%       postnans    # of days in date_range which are after last date in inputs

    if (length(date_range) == 2), date_range = [date_range(1),1,1; date_range(2),12,31]; end
    if (calendar_length(calendar) == 360)
        date_range(:,3) = min(30, date_range(:,3));     % in case day-of-month is 31...
        start_num = datenum360(start_vec);
        nc_dnums = floor(start_num + days_since); % dnums in nc file
        dnums = datenum360(date_range(1,:)):datenum360(date_range(2,:));  % dnums of desired output range.
        
    elseif (calendar_length(calendar) == 365)
        start_num = datenum365(start_vec);
        nc_dnums = floor(start_num + days_since);
        dnums = datenum365(date_range(1,:)):datenum365(date_range(2,:));        
    elseif (calendar_length(calendar) == 365.25)
        start_num = datenum(start_vec);
        nc_dnums = floor(start_num + days_since);
        dnums = datenum(date_range(1,:)):datenum(date_range(2,:));        
    else
        throw(MException('ICSF:BAD_CALENDAR',sprintf('error:  unknown calendar type: %s\n', calendar)));
    end
            
    prenans = max(0, nc_dnums(1) - dnums(1));
    postnans = max(0, dnums(end) - nc_dnums(end));

        % remove dnums past end of file
    nc_dnums(nc_dnums > dnums(end)) = [];
    nc_dnums(nc_dnums < dnums(1)) = [];
    
    keepers_ix = floor(nc_dnums - start_num - floor(days_since(1))) + 1;
end

