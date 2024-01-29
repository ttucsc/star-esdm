function mystruct = keyvalue_read(fname)
%   reads a key/value file into a struct.
%       numeric values are returned as doubles
%       text values are returned as strings.
%       text can span lines by ending a line with \  (newline will be included in text)
%       or by using a quoted string (newline will be replaced by space).  
%       leading and trailing blanks are removed, unless the string is quoted.
%       quoted numeric strings are returned as strings, not numbers.
%
%       Needs a way to add comments to a file.  e.g., ignore blank lines, ignore everything after unquoted # or % or //
%       probably add input parameter comment_seq to define comment start.
%

    fid=fopen(fname);
    mystruct=[];
    if (fid<0), error('error opening file %s',fname); end
    while (~feof(fid))
        kwd = fscanf(fid, '%s',1);
%        fprintf('kwd: %s :', kwd);
        val = rest_of_line(fid);
%        fprintf('\tval: %s\n', val);
        if (strlength(kwd)>0)
            mystruct.(kwd) = val;
        end
    end
end

function val = rest_of_line(fid)

    val='';
    done = false;
    quoting = false;
    quoted = false;
    qchar='';
    while (~done)
        txt = fgets(fid);
        if (feof(fid)), break; end
        txt =strtrim(txt);
        if (strlength(txt)==0), return; end            
        if (txt(end)=='\')
            txt=txt(1:end-1);
            if (strlength(val)==0)
                val=sprintf('%s\n',txt);
            elseif (val(end)==newline)
                val=sprintf('%s%s\n',val,txt);
            else
                val=sprintf('%s %s\n',val,txt);
            end
            done = false;
        elseif (quoting)
            if (txt(end)==qchar)
                txt=txt(1:end-1);
                quoting = false;
                quoted  = true;
                done    = true;
            end
            if (strlength(val)==0)
                val = txt;
            else
                val = sprintf('%s %s',val,txt);
            end
        else
            val=sprintf('%s%s',val,txt);
            if (val(1)=='''' || val(1)=='"')
                quoting = true;
                qchar   = val(1);
                if (val(end)==qchar)
                    val=val(2:end-1);
                    quoting = false;
                    quoted  = true;
                    done    = true;
                else
                    val=val(2:end);
                    done = false;
                end
            else
                done = true;
            end
        end
    end
    
    if (~quoted)
        dval = str2double(val);
        if (~isnan(dval))
            val = dval;
        end
    end
end
