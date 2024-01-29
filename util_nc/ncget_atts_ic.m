function [ attrs ] = ncget_atts_ic( fname, varName, group, atts )
%
%   
%   Returns attributes and data associated with group and varName in netcdf file fname
%
%   Inputs:
%       fname           either netcdf filename or ncid of open file or gpid of group inside a netcdf file
%       varName         variable name or varid of variable in fname
%                           if empty or missing, returns attributes for file or group
%       group           group within fname to get results for.  
%       atts            optional list of attributes to return.  If missing or empty, returns all attributes

    if (~exist('varName','var'))
        varName=[];
    end
    if (~exist('group','var'))
        group=[];
    end
    if (~exist('atts','var'))
        atts={};
    elseif (ischar(atts))
        atts={atts};
    end
    
    NC_GLOBAL = netcdf.getConstant('NC_GLOBAL');
    
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
    
    attrs = struct();
    if (isempty_s(varName))       % just return global attributes
        if (isempty_s(atts))
            [~,~, ngatts] = netcdf.inq(id);
            for i=0:(ngatts-1)
                [atname,attval] = getAttByNum(id, NC_GLOBAL, i);
                attrs.(atname) = attval;
            end
        else
            for i=1:length(atts)
                try
                    attrs.(atts{i})=netcdf.getAtt(ncid, NC_GLOBAL, atts{i});
                catch
                    throw(MException('ICSF:BAD_ATTR',sprintf('netcdf error:  no such attribute: %s',atts{i})));
                end
            end
        end
            
    else
        if (~ischar_s(varName))        
            varid = varName;
        else
            varid = netcdf.inqVarID(id,varName);
        end
        attrs = getAttrs(id,varid,atts);
            
    end
    
        % now copy requested atts to output attrs
        
    if (opened)
        ncclose_ic(ncid);
    end
end
function [attname, attval] = getAttByNum(id, varid, attid)
    attname = netcdf.inqAttName(id, varid, attid);
    attval = netcdf.getAtt(id,varid,attname);
end

function attrs = getAttrs(id,varid, atts)
    if (~isempty_s(atts))
        for i=1:length(atts)
            try
                attname = matlab.lang.makeValidName(atts{i},'Prefix','ZZZ');
                attrs.(attname)=netcdf.getAtt(id, varid, atts{i});
            catch
                throw(MException('ICSF:BAD_ATTR',sprintf('netcdf error:  no such attribute: %s',atts{i})));
            end
        end
    else            
        [~,~,~,natts] = netcdf.inqVar(id,varid);
        for i=0:(natts-1)
            [attname,attval] = getAttByNum(id,varid,i);
            attname=matlab.lang.makeValidName(attname,'Prefix','ZZZ');
            attrs.(attname) = attval;
        end
    end
end

