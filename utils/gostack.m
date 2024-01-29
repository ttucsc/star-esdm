function gostack(n)
%   jumps IDE to current line in program execution.
%   optional input n jumps to current line in nth stack entry.
%   gostack(1) is the same as gostack  (i.e., jump to current line)
%       NOTE:  1 is added to n

    if (~exist('n','var') || isempty(n)), n=1; end
    if (ischar(n) || isstring(n))
        nn=str2double(n);
        if (~isnan(nn)), n=nn; end
    end
    dbstack(1);
    ST=dbstack(n);
    if (isempty(ST))
        if (n==1)
            fprintf("no program running\n");
        else
            fprintf("no program at stack depth %d\n", n);
        end
        return
    end
    fname=ST(1).file;
    linenum=ST(1).line;
        
    matlab.desktop.editor.openAndGoToLine(which(fname), linenum);
end