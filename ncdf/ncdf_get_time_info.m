function [time_units, calendar, start_vec, days_since, timename, tstamps] = ncdf_get_time_info(ncobj, varName, do_load)
%
%    
%   Returns start and count vectors to read
%
%   Inputs:
%       ncobj           either netcdf filename or ncdf object or gpid of group inside a netcdf file
%       varName         variable name or varid of variable in fname
%                           if not 'time', then will look up time dimension of varName 
%       do_load         optional.   if true, then updates ncobj with time's vdata;  
%                                   if false, reads timevals from file but does not update time variable
%                                   if not present, returns existing timevals from ncobj.
%
%   Outputs:
%           time_units  units attribute of time
%           calendar    calendar attribute of time
%           start_vec   [year,month,day, hh,mm,ss] of timebase in time_units.
%           days_since  days_since start_vec
%           timename    actual time variable name.
%           tstamps     days-since as Matlab time values (days, and fraction of days, since 0/0/0000
%
%               NOTE:  days_since and tstamps are returned as doubles, regardless of the datatype in the file.


    if (~isa(ncobj,'Group'))
        ncobj = ncdf(ncobj);
    end
    
    if (~exist('varName','var') || isempty_s(varName))
        varlist = ncobj.varlist();
        if (any(strcmp(varlist,'time')))
            timename='time';
        elseif (any(strcmp(varlist,'Time')))
            timename='Time';
        else
            timename = 'TIME';
        end
    else
        timename = ncdf_get_timename(ncobj, varName);
    end
    
    time_units=strings(0,0);
    calendar=strings(0,0);
    start_vec = [];
    days_since = [];
    tstamps = [];
    try
        t=ncobj.get(timename);
        time_units = string(t.getattvalue('units'));
        try
            calendar = string(t.getattvalue('calendar'));
        catch
        end
        start_vec = nc_parse_date_str(time_units);
        if (nargout > 3 || do_load)
            if (~exist('do_load','var') || ~do_load)
                if (isempty(t.vdata))
                    days_since = double(ncobj.readvar(timename));
                else
                    days_since = double(t.vdata);
                end
            elseif (do_load)
                days_since = double(ncobj.loadvar(timename));
            else
                days_since = double(ncobj.readvar(timename));
            end
            if (nargout > 5)
                tstamps = datenum_cal(start_vec, calendar) + days_since; 
            end
        end
    catch
    end
end
