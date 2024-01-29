function dstr = datestr_cal(dnums, calendar, pivotyear_or_fmt, fmt)
% function dstr = datestr_cal(dnums, calendar, pivotyear_or_fmt, fmt)
%
% if pivotyear is char or string, it contains the format.
% if numeric, then dnums is treated as "days since" start of pivotyear.
% if pivotyear is a single number, if less than 10000 is treated as year, otherwise it's treated as a datenum in the
% calendar space given by calendar.
% this allows it to work with netcdf dates, given as "days since yyyy-mm-dd."
    if (nargin < 3)
        pivotyear_or_fmt='yyyy-mm-dd';
    elseif (nargin < 4)
        fmt = 'yyyy-mm-dd';
    end
    if (~isa(dnums,"double"))
        dnums=double(dnums);    
    end
    
    if (ischar(pivotyear_or_fmt) || isstring(pivotyear_or_fmt))
        fmt = pivotyear_or_fmt;
    else
        if (length(pivotyear_or_fmt) == 1)
            if (pivotyear_or_fmt < 10000)
                pivotyear = pivotyear_or_fmt;
            else
                day1 = pivotyear_or_fmt;
                dnums = dnums + day1;
            end
        else
            dnums = dnums + datenum_cal(pivotyear_or_fmt, calendar);
        end
    end
    
    if (calendar_length(calendar) == 360)
        if (exist('pivotyear','var') && exist('fmt','var'))
            dstr = datestr360(dnums, pivotyear, fmt);
        else
            dstr = datestr360(dnums, fmt);
        end
    elseif (calendar_length(calendar) == 365)
        if (exist('pivotyear','var') && exist('fmt','var'))
            dstr = datestr365(dnums, pivotyear, fmt);
        else
            dstr = datestr365(dnums, fmt);
        end
    else
        fmt = pivotyear_or_fmt;
        dstr = datestr(dnums, fmt);
    end    
    dstr = string(dstr);
end

