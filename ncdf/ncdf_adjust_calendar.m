function [out_start_vec, out_days_since, out_data, varNames] = ncdf_adjust_calendar(ncobj, out_calendar, varNames, do_update_time, out_start_vec, do_random, calendar, time_units, timevals)
%
%   Adjusts data and calendar in an ncdf object to a new calender, inserting or deleting days as needed, to switch
%   between standard, 365-day and 360-day calendars.
%
%   Required inputs:
%       ncobj           ncdf object, usually an ncdf of a file
%       out_calendar    desired output calendar
%       
%   optional inputs (read from file)
%       varnames                            list of variables in file to adjust. if missing or blank, adjusts every 
%                                               variable with timename as a dimension.
%                                               NOTE:  may not work right if multiple time variables and/or groups with
%                                               time-dependent data...
%       do_update_time  [false]             updates time variable as well as deleting or inserting days.
%                       EXCEPT:                 set to true if varNames includes time.
%       out_start_vec   [orig_start_vec]    "days_since" date for time's units.  
%       do_random       [false]             used to insert random days going from 360 to 365- or 365.25-day years
%       calendar
%       time_units
%       time_vals
%
%   NOTE:
%       do_update_time should probably be set to false if not adjusting all time-dependent variables, or if adjusting
%       them with multiple calls.
%       If adjusting them with sequential calls, set it to true for the last one.
%       To update ALL time dependent variables in 1 call, leave varnames blank. In
%       that case, the code will find all time dependent variables and adjust them, including the time variable.
%       Alternately, manually find all time-dependent variables and pass in a list of variables,  
%       
%       Alternately, pass in the original time info for specified variables on each call.
%
%   Outputs:
%       ncobj           is modified in place (ncobj is a handle, so doesn't need to be returned to modify contents)
%       out_start_vec 
%       out_days_since
%       out_data        adjusted variables' data.  
%                           If 1 variable only, is variable's data.
%                           If multiple variables, is cell array, 1 for each variable
%       varNames        string array adjusted variable names, in same order as cell array out_data.
    
%   [tunits, cal, start_vec, tvals] = ncdf_get_time_info(ncobj, 'time');
    [tunits, cal, ~,         tvals, timename] = ncdf_get_time_info(ncobj, 'time');
    if (~exist('calendar','var') || isempty_s(calendar))
        calendar = cal;
    end
    if (isempty_s(calendar))
        throw(MException('ARRMV2:NO_CALENDAR','ncdf_adjust_calendar: error:  no calendar info available')); 
    end
    
    if (~exist('time_units','var') || isempty_s(time_units))
        time_units = tunits;
        start_vec = nc_parse_date_str(time_units);      % already got it above.
    else
        start_vec = nc_parse_date_str(time_units);
    end
    
    if (~exist('timevals','var') || isempty_s(timevals))
        timevals = tvals;
    end
    if (~exist('do_random','var') || isempty_s(do_random))
        do_random = false;       % remove/insert at random locations in file so same day isn't always removed or inserted
    end
    if (~exist('do_update_time','var') || isempty_s(do_update_time))
        do_update_time = false;
    end        
    if (~exist('out_start_vec','var') || isempty_s(out_start_vec))
        out_start_vec = start_vec;
    end
    if (~exist('varNames','var') || isempty(varNames))
        varNames = find_time_dependent_variables(ncobj, timename);
        do_update_time = true;
    else
        varNames = string(varNames);
        if (any(strcmp(varNames,timename))), do_update_time = true; end
    end
        
        % get time dimension's ID
        
    nvars = length(varNames);
    if (nvars > 1)
        out_data = cell(nvars,1);
    end
    for i=1:nvars
        varName = varNames(i);
        v = ncobj.get(varName);
        dimid = v.findix('time','Dimensions');

        data = ncobj.getvardata(varName);
        if (isempty_s(data)), data = ncobj.readvar(varName); end
        if (is_time_dependent(ncobj, varName, timename))
            [ out_start_vec, out_days_since, adjusted_data ] = ncadjust_calendar( start_vec, timevals, data, calendar, out_calendar, do_random, out_start_vec, dimid  );
        end
        v.put('vdata',adjusted_data);
        if (nvars == 1)
            out_data = adjusted_data;
        else
            out_data{i} = adjusted_data;
        end
    end
    
    if (do_update_time)
        out_units = sprintf('days since %s', datestr(datenum(out_start_vec), 'yyyy-mm-dd HH:MM:SS'));
        ncobj.put('/Variables/time/Attributes/units',out_units);
        ncobj.put('/Variables/time/Attributes/calendar',out_calendar);
        ncobj.put('Variables/time/vdata',out_days_since);
        ncobj.put('Dimensions/time/Length',length(out_days_since));
        
            % check for "calendar_type" attribute;  some netcdf files (GFDL perfect model, e.g.) have both calendar and
            % calendar_type attributes.
        try
            ncobj.get("/Variables/time/Attributes/calendar_type");
            ncobj.put('/Variables/time/Attributes/calendar_type',out_calendar);
        catch
        end
    end    
end

function varNames = find_time_dependent_variables(ncobj, timename)
% returns string array containing all time-dependent variables EXCEPT timename.
%
    myvarlist = string(ncobj.varlist());
    nv = length(myvarlist);
    keepers = false(nv,1);
    
    for i=1:nv
        keepers(i) = is_time_dependent(ncobj, myvarlist(i), timename); 
    end
    varNames = myvarlist(keepers);
end

function tf = is_time_dependent(ncobj, varname, timename)

    v = ncobj.get(varname);
    if (~isempty(v.Dimensions))
        dimnames = string({v.Dimensions.Name});
        tf = any(strcmp(timename, dimnames));
    else
        tf = false;
    end

end

