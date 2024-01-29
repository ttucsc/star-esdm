function timevals = get_timevals_cal(yrs, pivot_yr, calendar)     
% timevar = get_timevals_cal(yrs, pivot_yr)
%
% returns time variable for given years years, w/ base of 'days since pivot_yr-01-01, 00:00:00'
%
%   Inputs:
%       yrs         [start_year, end_year] , years given as yyyy
%       pivot_yr    base year yyyy (for netcdf time attribute 'days since yyyy-01-01 00:00:00'
%       calendar    'standard','360-day','365-day', 'julian',etc..
%
    
    dnum0 = datenum_cal([pivot_yr,1,1],calendar);
    dnum1 = datenum_cal([yrs(1)  ,1,1],calendar);
    if (calendar_length(calendar) == 360)
        dnum2 = datenum_cal([yrs(end),12,30], calendar);
    else
        dnum2 = datenum_cal([yrs(end),12,31], calendar);
    end
    timevals = (dnum1:dnum2) - dnum0;
end
