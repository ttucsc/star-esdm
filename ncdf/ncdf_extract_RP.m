function RP = ncdf_extract_RP(ncobj, varargin)
%
%       If ncobj is a Group, then varargin{1} must be the varname.
%       Otherwise, code assumes ncobj is either a ncdf object or netcdf filename.
%
    RP = [];
   if (strcmp(ncobj.Name,"RunParams"))
       if (~is_climate_variable(varargin{1})), error("error:  1st parameter must be a climate variable name"); end
       RP = extract_RP(ncobj, varargin{1}, varargin{2:end});
   else
       if (isstring(ncobj) || ischar(ncobj))
           ncobj = ncdf(ncobj);
       end
       if (~isempty(varargin) && ((isstring(varargin{1}) || ischar(varargin{1})) && is_climate_variable(varargin{1})))
            varname = varargin{1};
            if (length(varargin)>1)
                varargin = varargin{2:end};
            else
                varargin = cell(0);
            end
       else
            varname=find_climate_varname(ncobj.varlist,1);
       end
       
       for i=1:length(ncobj.Groups)
           RP = ncdf_extract_RP(ncobj.Groups(i), varname, varargin{:});
           if (~isempty(RP))
               return;
           end
       end
   end               
end

function RP = extract_RP(ncobj, varname, varargin)

    
    RP = ARRM_V2_RunParams(varname);      % create a DataParams object with default values.

        % copy all DataParam fields from the group's attributes
    for i=1:length(ncobj.Attributes)
        nam = ncobj.Attributes(i).Name;
        val = ncobj.Attributes(i).Value;
%       fprintf("%s: %s\n", nam, string(val));
        if (isprop(RP,nam))
            try
                RP.(nam) = val;     % in try/catch block in case nam is a dependent property.  
                                    % there is a way to test if something is a dependent property, but it only works for
                                    % non-hidden properties, so we just use a try/catch block here instead...
                                    % Ref: https://www.mathworks.com/help/matlab/matlab_oop/getting-information-about-properties.html
            catch
            end
        end
    end
    
        % then apply any needed updates specified in varargin.
    RP = RP.update(varargin{:});
    
end
