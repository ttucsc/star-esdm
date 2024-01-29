function [timevals, calendar, units] = ncget_timeinfo_ic( fname, group)
%
%   Returns time info from fname (& group is specified).  
%   Assumes variable is named 'time', has an attribute 'calendar' and an attribute 'units'
%
%   Inputs:
%       fname           either netcdf filename or ncid of open file or gpid of group inside a netcdf file
%       group           optional.  group name within fname
%
%   returns the raw timevals.
%
%   calendar, units or timevals are set to empty if it cannot retrieve the information.
%
%   
%
%   To get datenums from the  timvals, calendar and units:
%
%           [from_vec, timescale, isUTC] = nc_parse_date_str(units);
%           dnums = datenum_cal(from_vec, calendar) + timevals;



    if (~exist('group','var'))
        group=[];
    end
    if (ischar_s(fname))
        ncid = ncopen_ic(fname,'NOWRITE');
        opened = true;
    else
        ncid = fname;
        opened = false;
    end    
    id = ncid;
    if (~isempty_s(group))
        gpid = netcdf.inqNcid(ncid,group);
        id = gpid;
    end
    
    timevals = [];
    calendar = '';
    units = '';
    
    try
    
        timid = ncget_varinfo_ic( id, [], 'time' );

        timevals = netcdf.getVar(id, timid); 
        
            % get the units, calendar, etc, and timevals from the netcdf file.
        units = netcdf.getAtt(id, timid,'units');
        calendar = netcdf.getAtt(id, timid, 'calendar');
    catch
    end
 
    if (opened)
        ncclose_ic(ncid);
    end    
end

