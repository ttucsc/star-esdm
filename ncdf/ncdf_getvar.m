function [vals, varpath, varname, ncobj] = ncdf_getvar(ncobj, varpath, read_ok)
% ncdf_gervar:  returns data for varpath in ncobj.  If ncobj hasn't read the data yet, reads it from the file.
% ncobj         is ncdf object, or netcdf filename.
% varpath       dimension name (or qualified ncdf path)
% read_ok       bool.  if true, loads from file if data is not present.  Otherwise throws an error.
    
%       NOTE:  this should use ncdf_gervarname(...), ian!
%
% 6/7/22 icsf modified to return lats & lons as doubles.

    latnames={'lat','latitude','lats'};
    lonnames={'lon','longitude','lons'};
    timenames={'time','times'};
    
    
    if (ischar_s(ncobj))
        ncobj = ncdf(ncobj, "do_create",false);
        read_ok = true;
    elseif (~exist('read_ok','var') || isempty(read_ok))
        read_ok = false;
    end
    
    path = split(varpath,'/');
    varname=path{end};
    vp_orig = varpath;
    
    varlist = ncobj.varlist();
    if (any(strcmpi(varpath,varlist)))     % see if variable exists.
        try
            vix=find(strcmpi(varpath,varlist),1);
            vals = ncobj.getvardata(varlist{vix});
            if (isempty(vals) && read_ok)
                vals = ncobj.readvar(varlist{vix});
            end            
            return;
        catch me
            if (~read_ok), rethrow(me);end
        end
        try
            vals = ncobj.readvar(varlist{vix});
            return;
        catch
        end
        throw(MException('NCDF:NODATAYET',sprintf('error:  variable %s exists but not read in yet',varpath)));
    end    
        % no exact match.  look for alternate name for standard dimensions.
        
    is_latlon = false;
    if     (any(strcmpi(varname,latnames)))
        names = latnames;
        is_latlon = true;
    elseif (any(strcmpi(varname,lonnames)))
        names = lonnames;
        is_latlon = true;
    elseif (any(strcmpi(varname,timenames)))
        names = timenames;
    else
        names = {};
    end
    for i=1:length(names)
        path{end} = names{i};
        varpath = join(path,'/');
        varpath=varpath{1};
        if (any(strcmpi(varpath,varlist)))
            try
                varname = names{i};
                vals = ncobj.getvardata(varpath);
                if (isempty(vals) && read_ok)
                    vals = ncobj.readvar(varpath);
                end
                if (is_latlon)
                    vals = double(vals);
                end
                return;
            catch me
                if (~read_ok), rethrow(me);end
            end
            if (read_ok)
                try
                    vals = ncobj.loadvar(varpath);
                    if (is_latlon)
                        vals = double(vals);
                    end
                    return;
                catch
                end
            end
        end
    end 
    throw(MException('NCDF:NOSUCHVAR',sprintf('error:  no match for variable %s',vp_orig)));
end


