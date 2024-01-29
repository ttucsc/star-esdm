function [units, ncvar, varname, atts] = ncdf_getvar_info(ncobj, varname)
% ncdf_gervar_info:  returns data for varpath in ncobj.  If ncobj hasn't read the data yet, reads it from the file.
%   Inputs:
%       ncobj          is ncdf object, or netcdf filename.
%       varname        variable name
%   Outputs
%       units           variable's units attribute value
%       ncvar           ncdf Variable object
%       varname         varname found
%       atts            variable's attributes array.
    
    if (ischar_s(ncobj))
        ncobj = ncdf(ncobj, "do_create",false);
    end
    
            % this needs to parse the path and go down the path to work for embedded groups, Ian
%     path = split(varname,'/');
%     varname=path{end};
    
    varlist = ncobj.varlist();
    if (any(strcmpi(varname,varlist)))     % see if variable exists.
        vix=find(strcmpi(varname,varlist),1);
        if (~isempty(vix))
            ncvar=ncobj.Variables(vix);
            attlist = ncvar.attlist();
            ix = find(strcmp(attlist,"units"),1);
            if (~isempty(ix))
                units=ncvar.getattvalue("units");
            else
                units=strings(0);
            end
        end
        atts = ncvar.Attributes;
        return;
    else
        error("error:  variable %s not found\n", varname)
    end    
end


