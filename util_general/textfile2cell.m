function mytext = textfile2cell(fname,skip_empty)
%   reads a text file into a cell array, 1 line of the file per element..
%   lines ending with a '\' will be appended to the previous line (with newline), rather than becoming its own element.
%   
%   Inputs:
%       filename    name of file to read
%       skip        boolean.  if true, removes any lines which are empty or whitespace-only.  default:  false

    if (~exist('skip_empty','var') || isempty(skip_empty)), skip_empty = false; end
    
    fid=fopen(fname);
    try
        mytext=[];
        nlines=0;
        if (fid<0), error('error opening file %s',fname); end
        while (~feof(fid))
            txt = get_one_line(fid);
            if (txt == -1) 
                break; 
            end
%             fprintf('\tval: %s\n', txt);
            if (~skip_empty || nonwhitelength(txt)>0)
                nlines=nlines+1;
                mytext{nlines} = txt;
%                fprintf('%2d : %s\n', nlines, mytext{nlines});
            end
        end
        if (isrow(mytext)), mytext=mytext'; end
    catch m
        fclose(fid);
        rethrow(m);
    end
    fclose(fid);
end

function val = get_one_line(fid)

    val='';
    done = false;
    while (~done)
        txt = fgets(fid);
        if (txt(end)==newline), txt=txt(1:end-1); end
        if (txt==-1)
            if (isempty(val)), val=-1; end
            break; 
        end
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
        else
            val=sprintf('%s%s',val,txt);
            done = true;
        end
    end
end
