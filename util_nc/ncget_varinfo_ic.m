function [ varid, varName, xtype, dimids, natts ] = ncget_varinfo_ic( fname, group, varName )
%
%   
%   Returns attributes and data associated with group and varName in netcdf file fname
%
%   Inputs:
%       fname           either netcdf filename or ncid of open file or gpid of group inside a netcdf file
%       varName         variable name or varid of variable in fname
%                           if empty or missing, returns attributes for file or group
%       group           group within fname to get results for.  

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
    if (ischar_s(varName))
        varid = netcdf.inqVarID(id, varName);
    else
        varid = varName;
    end
    [varName, xtype, dimids, natts] = netcdf.inqVar(id, varid);
    
    if (opened)
        ncclose_ic(ncid);
    end
end

