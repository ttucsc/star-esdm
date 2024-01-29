function dstr = datestr360(dnums360, pivotyear, fmt)
%returns datestr for datenums based on a 360-day calendar.
%   alternately, dnums365 can be days since Jan 1 of pivotyear
%   so can work with netcdf dates, which are given in days-since some start date.

    if (~exist('fmt','var') || isempty(fmt)), fmt = 'yyyy-mm-dd'; end
    if (exist('pivotyear','var') && ~isempty(pivotyear))
        if (isnumeric(pivotyear))
            dnums360 = dnums360 + datenum360([pivotyear,1,1]);
        else
            fmt = pivotyear;
        end
    end
    
    
    dvecs = datevec360(dnums360);
    dstr = strings(size(dnums360));
    
    if (strcmp(fmt,"yy-mm-dd"))
        fmt = "%02d-%02d-%02d";
        dvecs(:,1) = mod(dvecs(:,1),100);
        cols=1:3;
    elseif (strcmp(fmt,"yy-mm-dd HH"))
        fmt = "%02d-%02d-%02d %02d";
        dvecs(:,1) = mod(dvecs(:,1),100);
        cols=1:4;
    elseif (strcmp(fmt,"yy-mm-dd HH:MM"))
        fmt = "%02d-%02d-%02d %02d:%02d";
        dvecs(:,1) = mod(dvecs(:,1),100);
        cols=1:5;
    elseif (strcmp(fmt,"yy-mm-dd HH:MM:SS"))
        fmt = "%02d-%02d-%02d %02d:%02d:%02d";
        dvecs(:,1) = mod(dvecs(:,1),100);
        cols=1:6;
    elseif (strcmp(fmt,"yyyy-mm-dd"))
        fmt = "%04d-%02d-%02d";
        cols=1:3;
    elseif (strcmp(fmt,"yyyy-mm-dd HH"))
        fmt = "%04d-%02d-%02d %02d";
        cols=1:4;
    elseif (strcmp(fmt,"yyyy-mm-dd HH:MM"))
        fmt = "%04d-%02d-%02d %02d:%02d";
        cols=1:5;
    elseif (strcmp(fmt,"yyyy-mm-dd HH:MM:SS"))
        fmt = "%04d-%02d-%02d %02d:%02d:%02d";
        cols=1:6;
    else
        error("error:  datestr360 only accepts formats of ""yy-mm-dd"" or ""yyyy-mm-dd"" with optional HH:MM:SS");
    end
    for i=1:length(dnums360)
        dstr(i) = sprintf(fmt, dvecs(i,cols));
    end
end

