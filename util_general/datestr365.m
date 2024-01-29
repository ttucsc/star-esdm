function dstr = datestr365(dnums365, pivotyear, fmt)
%returns datestr for datenums based on a 365-day calendar.
%   alternately, dnums365 can be days since Jan 1 of pivotyear
%   so can work with netcdf dates, which are given in days-since some start date.


    if (~exist('fmt','var') || isempty(fmt)), fmt = 'yyyy-mm-dd'; end
    if (~exist('pivotyear','var'))
        dstr = datestr(daynoleap2datenum(dnums365-1,0),fmt);           % dnums are days since 0-0-0;
    else
        if (isnumeric(pivotyear))
            dstr = datestr(daynoleap2datenum(dnums365,pivotyear),fmt);     % dnums are days since 1-1-pivotyear.
        else
            fmt = pivotyear;
            dstr = datestr(daynoleap2datenum(dnums365-1,0),fmt);     % dnums are days since 1-1-pivotyear.
        end
    end
end

