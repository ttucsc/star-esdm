function [dirinfo, fullnames] = dirmatch(fns)
% returns output of dir for fns.  fns can be:
%   directory name
%   one or more filenames  (as char array, string array or cell array of chars or strings)
%   one or more linux-style file templates.
%       or any combination of the above
%
% Also returns fully qualified filenames of each file in fullnames;
%
%   NOTE 1:  This will match constructs such as "abc[1-3]*" and "filename with spaces*"  
%               (do NOT escape spaces in the template;  the code will escape any isspace(...) characters.
%   NOTE 2:  If a file matches more than one of the filenames/templates passed in in fns, it will appear multiple times in the output.
%
%   NOTE 3:  this routine works only on Unix-like systems (linux, MacOS, etc.).  Will probably not work on Windows  machines.

    fns = string(fns);
    for i=1:length(fns)
        tmpl = fns(i);
        if (~contains(tmpl,'?') && contains(tmpl,'[') && contains(tmpl,'{'))
            if (~exist('dirinfo','var'))
                dirinfo = dir(tmpl);
            else
                dirinfo = [dirinfo; dir(tmpl)]; %#ok<AGROW>
            end
        else
            tmpl = esc_spaces(tmpl);    % escape all spacees in template.
            [~,myfns] = system("ls -1 "+tmpl);
            myfns = string(strsplit(myfns,"\n"));
            myfns(strlength(myfns)==0)=[];
            for j=1:length(myfns)
                if (~exist('dirinfo','var'))
                    dirinfo = dir(myfns(j));
                else
                    dirinfo = [dirinfo; dir(myfns(j))]; %#ok<AGROW>
                end
            end
        end
    end
    if (nargout > 1)
        nfiles=length(dirinfo);
        if (nfiles==0)
            fullnames=strings(0,0);
        else
            fullnames = strings(nfiles,1);
            for i=1:length(dirinfo)
                fullnames(i) = fullfile(dirinfo(i).folder, dirinfo(i).name);
            end
        end
    end            
end

function outstr = esc_spaces(instr)

    instr = char(instr);
    outstr = "";
    for i=1:length(instr)
        c = instr(i);
        if (isspace(c))
            outstr=outstr+"\"+c;
        else
            outstr = outstr+c;
        end
    end
end
            
    


