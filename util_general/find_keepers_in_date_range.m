function [keepers, keepers_ix, dest_ix, outnpts] = find_keepers_in_date_range(date_range, calendar, start_vec, days_since)
% returns info on which values to read from input file, and where to put them in the output data.
%   This will allow you to read a limited set of data from the input file, and will handle input data
%   that isn't in sequential date order, has gaps in it or has overlapping (duplicate) dates.
%
%   Inputs:
%       date_range  either [start_year, end_year] or [start_yyyy_mm_dd; end_yyyy_mm_dd] of desired range of data
%       calendar    calendar string of time units
%       start_vec   datevec of day 0 for timevals
%       days_since  vector of days-since-start_vec
%       
%   Outputs:
%       keepers     boolean array of input dates within specified date range (relative to the days_since vector)
%                       user should read the data and retain only the keepers.
%       keepers_ix  indexes of keepers (in days_since).
%       dest_ix     indexes of where to put keepers in output array so that output array covers the entire date_range                       
%                       assumption is that user wants data specified in date range, and will create an array covering 
%                       the full date_range.  If data covers the full date range, then dest_ix will be 1:npts. 
%                       If available data doesn't cover the full date range, then dest_ix tells user where to put the 
%                       data in their output array.
%       outnpts     length of desired output array (i.e., long enough to cover the entire date_range, even if the 
%                       available data doesn't cover the full date range.)
%
%           If you read the entire dataset (all dates), then:
%               data(keepers) should go into dest(dest_ix), where dest is initialized to nan(outnpts,1);
%
%           If you want to read only the data between min(keepers_ix) and max(keepers_ix)
%               ix1 = find(keepers,1);
%               ix2 = find(keepers,1,'last');
%                   or
%               ix1 = keepers_ix(1);
%               ix2 = keepers_ix(end);
%
%               indata = read(data between ix1 and ix2)
%               dest(dest_ix) = indata(keepers(ix1:ix2))

    if (length(date_range) == 2), date_range = [date_range(1),1,1; date_range(2),12,31]; end
    if (calendar_length(calendar) == 360)
        date_range(:,3) = min(30, date_range(:,3));     % in case day-of-month is 31...
        day0 = datenum360(start_vec);              % datenum of start_vec -- datenum of day 0 for days_since.
        nc_dnums = floor(day0 + days_since); % dnums in nc file
        dnum1 = datenum360(date_range(1,:));
        dnum2 = datenum360(date_range(2,:));  % dnums of desired output range.
        
    elseif (calendar_length(calendar) == 365)
        day0 = datenum365(start_vec);
        nc_dnums = floor(day0 + days_since);
        dnum1 = datenum365(date_range(1,:));
        dnum2 = datenum365(date_range(2,:));        
    elseif (calendar_length(calendar) == 365.25)
        day0 = datenum(start_vec);
        nc_dnums = floor(day0 + days_since);
        dnum1 = datenum(date_range(1,:));
        dnum2 = datenum(date_range(2,:));        
    else
        throw(MException('ICSF:BAD_CALENDAR',sprintf('error:  unknown calendar type: %s\n', calendar)));
    end
     
    outnpts = dnum2 - dnum1 + 1;
    keepers = (nc_dnums >= dnum1) & (nc_dnums <= dnum2);
    keepers_ix = find(keepers);
    
    dest_ix = nc_dnums - dnum1 + 1;     % destination for every input date.  dest's < 1 or >outnpts will be out of range
    dest_ix = dest_ix(keepers);         % destination for only the keepers.
                                        % data(keepers) should go into dest(dest_ix).
    
end

