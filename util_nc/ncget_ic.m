function [data] = ncget_ic(fname, varName, group, start, count)
%
%   
%   Returns attributes and data associated with group and varName in netcdf file fname
%
%   Inputs:
%       fname           either netcdf filename or ncid of open file or gpid of group inside a netcdf file
%       varName         variable name or varid of variable in fname
%                           can be cell array of possible variables e.g., {'lat','latitude'}
%                           note:  for latitude & longitude, will look for any of the versions in latnames or lonnames
%                           to look for a specific instance, pass it in as a single cell array {'lat'}
%       group           group within fname to get results for.  
%       start           vector of 3 indexes, ([x_index,y_index,time_index]), defining where to start.
%       count           vector of 3 counts, ([x_count,y_count,time_count]) defining how many values for that dimension to read.
%                           for netcdf files, x is usually longitude, y is usually latitude.
%
%                       % note:  if longitude wraps around 360 degrees, do you need to do 2 separate reads.

            % to look for various forms of lat, lon
    latnames={'lat','latitude','lats','Latitude','LATITUDE','LAT','LONGS'};
    lonnames={'lon','longitude','lons','long','longs','LONGITUDE','LON','LONG','LONGS'};

    if (ischar_s(varName) && (sum(strcmpi(latnames,varName))>0) )
        varName=latnames;
    elseif (ischar_s(varName) && (sum(strcmpi(lonnames,varName))>0) )
        varName=lonnames;
    end

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
    
    if (isnumeric(varName))        
        varid = varName;
    else
        if (ischar(varName)), varName={varName}; end
        
        for i=1:length(varName)
            ok=false;
            vname=varName{i};
            try
                varid = netcdf.inqVarID(id,vname);
                ok=1;
                break;
            catch
            end
        end
        if (~ok)
            if (opened), ncclose_ic(ncid); end
            allvars = sprintf('%s ', varName{:});
            if (ischar_s(fname))
                throw(MException('ARRMV2:BADVARNAME', sprintf('error:  can''t find var %s in netcdf file %s', allvars, fname)));
            else
                throw(MException('ARRMV2:BADVARNAME', sprintf('error:  can''t find var %s in netcdf file', allvars)));
            end
        end
    end
    if (exist('start','var'))
        if (exist('count','var'))
            data = netcdf.getVar(id,varid, start, count);
        else
            data = netcdf.getVar(id,varid, start);
        end
    else
        data = netcdf.getVar(id,varid); 
    end
    
    try
        ty=netcdf.inqFormat(ncid);
        if (strcmp(ty,'FORMAT_NETCDF4'))
            [~, fillVal] = netcdf.inqVarFill(id, varid);
        else
            fillVal      = netcdf.getAtt(id, varid,'_FillValue');
        end
        data(data==fillVal) = nan;
    catch
    end
    if (opened)
        ncclose_ic(ncid);
    end
end
