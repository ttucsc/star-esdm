function [timename, timix] = ncdf_get_timename(ncobj, varName)
% returns name and index of time dimension for varName.  Will look for case-insentive match for 'time'
%   If varName is not present, or varName is 'time', returns the time dimension's name and index for ncobj
%   If varName is not 'time', then timix is the index for time in the variable's list of dimensions.

    if (~exist('varName','var'))
        vname = "time"; 
    else
        path = split(string(varName),"/");
        vname = path(end);    
    end
    
    if (strcmpi(vname,"time"))
        dlist = ncobj.dimlist();
        timix = find(strcmpi("time", dlist),1);
        if (isempty(timix))
            timename = [];
        else
            timename = dlist(timix);
        end
    else    
        v = ncobj.get(varName);
        dlist = v.dimlist();
        timix = find(strcmpi("time", dlist),1);
        if (isempty(timix))
            timename = [];
        else
            timename = dlist(timix);
            if (length(path) > 1)
                path(end+1) = timename;
                timename = join(path,"/");
            end
        end
    end
end

