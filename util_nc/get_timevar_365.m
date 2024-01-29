function timevar = get_timevar_365(yrs, pivot_yr)     
% timevar = get_timevar_365(yrs, pivot_yr)
%
% returns time variable for 365-day years, w/ base of 'days since pivot_yr-01-01, 00:00:00'
%
%   Inputs:
%       yrs         [start_year, end_year] , years given as yyyy
%       pivot_yr    base year yyyy (for netcdf time attribute 'days since yyyy-01-01 00:00:00'
%

    istart = (yrs(1)-pivot_yr)*365;
    nyrs = yrs(2) - yrs(1)+1;
    days = (1:(365*nyrs))'-1;
    timevar = istart + days;
end
