function [varname, var] = ncdf_getvarname(ncobj, varpath)
% ncdf_gervarnamne:  returns variable's name, & optionally the variable object as well, for varpath.
    % ncobj         is ncdf object or filename.
    % varpath       variable name (or qualified ncdf path).  For time, lat, lon or day_of_year, finds actual variable.

    latnames={'lat','latitude','lats','latitudes'};
    lonnames={'lon','longitude','lons','longitudes'};
    timenames={'time','times'};
    doynames={'doy','day_of_year','dayofyear'};
    
    path = split(varpath,'/');
    varname=path{end};
%   vp_orig = varpath;
    
    if (ischar_s(ncobj))
        ncobj = ncdf(ncobj,"create_ok",false);
    end        
    
    varlist = ncobj.varlist();
    vix=find(strcmpi(varpath,varlist),1);        % see if variable exists.  
    if (~isempty(vix))
        varname = varlist{vix};
        if (nargout > 1)
            var = ncobj.getvar(varname);
        end
        return;
    end
        % no exact match.  look for alternate name for standard dimensions.
    if     (any(strcmpi(varname,latnames))), names = latnames;
    elseif (any(strcmpi(varname,lonnames))), names = lonnames;
    elseif (any(strcmpi(varname,timenames))), names = timenames;
    elseif (any(strcmpi(varname,doynames))), names = doynames;
    else,  names = {};
    end
    varname = [];
    for i=1:length(names)
        path{end} = names{i};
        varpath = join(path,'/');
        varpath=varpath{1};
        vix=find(strcmpi(varpath,varlist),1);        % see if variable exists.  
        if (~isempty(vix))
            varname = varlist{vix};
            if (nargout > 1)
                var = ncobj.getvar(varname);
            end
            return;
        end
    end
%   error('error:  no match for variable %s',vp_orig);

end


