function [yn, nc] = isnetcdf(fname, do_validate)
%   returns boolean for whether fname is a netcdf file.
%   if output nc is requested, or do_validate is true, then fname is opened and schema is read in to validate that the
%   file really is a netcdf file.  Otherwise, it just checks that the extension is ".nc" or ".NC".
%
%   if extension is not .nc,  or if nargout > 1, then tries to open it as an ncdf object.
%   This may take a little time if it is a complex netcdf file, as ncdf(...) reads the file's schema in. (does not read
%   the actual data, though.)
%
    nc=false;
    yn=false;
    try
        if (~isfile(fname)), yn=false; return; end
        [~,~,ext] = fileparts(fname);
        if ((~exist('do_validate','var') || ~do_validate) && strcmpi(ext,".nc") && nargout < 2), yn=true; return; end
        try
            nc=ncdf(fname);
            yn=true;
        catch
        end
    catch
    end
end