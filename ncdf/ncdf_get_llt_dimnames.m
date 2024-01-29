function [latname, lonname, timename, doyname, latix, lonix, timix, doyix] = ncdf_get_llt_dimnames(nc)
%   [latname, lonname, timename] = ncdf_get_llt_dimnames(nc)
%
%   returns actual dimension name for latitude, longitude and time dimensions in netcdf file nc
%
%   inputs
%       nc      either an ncinfo struct, an ncdf object, or the filename of a netcdf file.
%   Outputs
%       latname     first match to the obvious possibilities below.  
%       lonname         if not found, name returned is empty.
%       timename 

    if (isa(nc,'Variable'))
        v=nc;
        [latix,~,~,latname]  = v.diminfo('lat');
        [lonix,~,~,lonname]  = v.diminfo('lon');
        [timix,~,~,timename] = v.diminfo('time');
        [doyix,~,~,doyname]  = v.diminfo('doy');
    else
        latix = [];
        lonix = [];
        doyix = [];
        
        if (ischar_s(nc)), nc=ncdf(nc,'create_ok',false); end
        latname  = ncdf_getvarname(nc,'lat');
        lonname  = ncdf_getvarname(nc,'lon');
        timename = ncdf_getvarname(nc,'time');
        doyname  = ncdf_getvarname(nc,'doy');
        dims=nc.dimlist();
        if (~isempty(latname))
            latix = find(strcmp(dims, latname),1);
        end
        if (~isempty(lonname))
            lonix = find(strcmp(dims, lonname),1);
        end
        if (~isempty(timename))
            timix = find(strcmp(dims, timename),1);
        end
        if (~isempty(doyname))
            doyix = find(strcmp(dims, doyname),1);
        end
    end
end



