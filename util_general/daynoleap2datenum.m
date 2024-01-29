function dn = daynoleap2datenum(day, pivotyr)
%DAYNOLEAP2DATENUM Convert from days since, no leap to serial date number
%
% dn = daynoleap2datenum(day, pivotyr)
%
% A lot of model output is saved with a time scale of "days since
% YYYY-01-01 00:00:00", where every year has 365 day.  This function
% converts those dates into a serial date number.
%
% Input variables:
%
%   day:        number of days since pivot year
%
%   pivotyr:    pivot year, i.e. year that day count begins (on Jan 1)

% Determine which years in range are leap years

nyr = max(day./365);
yrs = pivotyr:(pivotyr + nyr);

isleap = @(x) (mod(x,4)==0 & mod(x,100)~=0) | mod(x,400) == 0;
leapyrs = yrs(isleap(yrs));

% Calculate date numbers

dayofyear = rem(day, 365);
yr = floor(day/365) + pivotyr;

dn = datenum(pivotyr, 1, 1) + day;
for ileap = 1:length(leapyrs)
  needsbump = yr > leapyrs(ileap) | (yr == leapyrs(ileap) & dayofyear >= 59);
  dn(needsbump) = dn(needsbump) + 1;
end
