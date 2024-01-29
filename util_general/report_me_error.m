function report_me_error(me, srcname, fid, do_rethrow)
%       for displaying caught error messages.  Use with:
%       try
%           ...
%       catch me
%           report_me_error(me,mfilename [, fid, do_rethrow])
%       end
%
%   Inputs:
%       me          error caught
%       scrname     name of file that caught the error.  mfilename() returns name of current file.
%       fid         fileID to write error message to.  [default: 2  stderr].
%                       can also be the name of a file to write to.  
%       do_rethrow  if true, rethrows error after reporting info. [false].
    if (~exist('fid','var') || isempty(fid)), fid=2; end
    st=dbstack();
    if (numel(st)>1)
        myprog = st(2).name;
        myline = st(2).line;
        fprintf(fid, "%s:%d ", myprog, myline);
    end
    if (exist("srcname","var") && ~isempty(srcname))
        srcname = basename(srcname);
        fprintf(fid, '%s:  caught exception:  %s\t\t%s\n', srcname, me.identifier, me.message);
    end
    
    if (isstring(fid) || ischar(fid))
        fid = fopen(fid, "a");
        do_close = true;
    else
        do_close = false;
    end
    msgtext=getReport(me);
    fprintf(fid, 'error traceback:\n%s\n', msgtext);
    fprintf(fid, '------\n');
    if (exist('do_rethrow','var') && do_rethrow)
        rethrow(me);
    end
    if (do_close), fclose(fid); end
    
end

