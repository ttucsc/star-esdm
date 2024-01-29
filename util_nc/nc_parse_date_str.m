function [from_vec, timescale, isUTC] = nc_parse_date_str(timeunits, calendar)
% returns starting datevec and the timescale from a date_str of the form "days since 1950-01-01 00:00:00"

%   from_vec    [yyyy, mm, dd, HH, MM, SS]
%   timescale   time steps per day.  1 for days, 24 for hours, 24*60 for minutes, etc.  0 if no info available
%               from date_str
%               NOTE:  timescale is approximate, and will depend on calendar, if "years since" or "months since"....
%                       so use timescale as a flag, rather than a conversion factor if timescale is < 1.

    dt=strsplit(timeunits);
    if (~exist('calendar','var') || isempty(calendar)), calendar='365_day'; end
    
    if (strcmp(dt{1},'days') )       
        timescale=1.0;
    elseif (strcmp(dt{1},'hours'))
        timescale=24.0;
    elseif (strcmp(dt{1},'months'))
        timescale=1/12.0;
    elseif (strcmp(dt{1},'years'))
        timescale=1/calendar_length(calendar);
    elseif (strcmp(dt{1},'minutes'))
        timescale=24*60;
    elseif (strcmp(dt{1},'seconds'))
        timescale=24*60*60;
    else
        timescale=0;
    end

    dat=sscanf(dt{3},'%d-%d-%d')';
    if (length(dt)>3 && ~isempty(regexp(dt{4},'\d{2}:\d{2}:\d{2}','ONCE')))
        tim=sscanf(dt{4},'%d:%d:%d')';
    else
        tim=[0,0,0];
    end

    from_vec = [dat,tim];

    if (strcmpi(dt{end},'UTC'))
        isUTC=true;
    else
        isUTC=false;
    end
end

