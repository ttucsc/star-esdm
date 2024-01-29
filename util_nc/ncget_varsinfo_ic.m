function [ varNames, varIDs, xtypes, dimids, natts ] = ncget_varsinfo_ic( fname, group)
%
%   
%   Returns varName and varIDs of all variables in netcdf file fname
%
%   Inputs:
%       fname           either netcdf filename or ncid of open file or gpid of group inside a netcdf file


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
    
    varIDs = netcdf.inqVarIDs(id);
    nvars = length(varIDs);
    varNames = cell(nvars,1);
    xtypes   = nan(nvars,1);
    dimids   = cell(nvars,1);
    natts    = cell(nvars,1);
    for i=1:nvars
        [varNames{i}, xtypes(i), dimids{i},natts{i}] = netcdf.inqVar(ncid,i-1);
    end
    if (opened)
        ncclose_ic(ncid);
    end
end

