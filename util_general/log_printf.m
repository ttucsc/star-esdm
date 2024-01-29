function log_printf(logname, console_id, tmplate, varargin)
% log_printf(logname, do_console, tmplate, varargin)
%
%   writes msg to log file identified by logname, and optionally to console as well.
%
%   Inputs:
%       logname     filename or filehandle to log file
%       console_id  true/false or 1 or 2.  1->stdout;  2->stderr
%       template    fprintf template
%       varargin    any arguments needed for template.

    if (console_id)
        conid = 1*console_id;   % in case it's logical true.
        fprintf(conid, tmplate, varargin{:});
    end
    if (ishandle(logname))
        fid = logname;
    else
        fid = fopen(logname, "a");
        if (~fid)
            error("error:  cannot open logname %s for writing", logname);
        end
        opened = true;
    end
    fprintf(fid, tmplate, varargin{:});
    if (opened)
        fclose(fid);
    end
end

