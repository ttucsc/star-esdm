function d = clean_string(instr, maxlen,replchar)
%   clean up string to remove non-alphanumeric chars
%   replaces them with underscores, then removes any double underscores
%   Useful for using cleaned up Station Name in a filename

    if (~exist('replchar','var') || isempty(replchar))
        replchar='_';
    end
    dupchar=repmat(replchar,1,2);
    if (iscell(instr))
        wascell = true;
        allstrings = false;
        d=instr;
    else
        allstrings=isstring(instr);
        d = cellstr(instr);
        wascell=false;
    end
    npts=length(d);
    for i=1:npts
        dd = d{i};
        if (isstring(dd))   % in case we had a cell array of strings
            wasstring = true;
            dd=char(dd);
        else
            wasstring = false;
        end
        dd(~isstrprop(dd, 'alphanum')) = replchar;
        p=strfind(dd,dupchar);
        while (~isempty(p)) 
            dd(p(1))=[];
            p=strfind(dd,dupchar);
        end
        if (exist('maxlen','var') && ~isempty(maxlen))
            dd(maxlen+1:end)=[];
        end
        if (wasstring)      % make it back into a string if changed to char.
            d{i}=string(dd); 
        else
            d{i}=dd;
        end
    end
    if (npts==1 && ~wascell)    % make single item if it was only 1 element long, and not a cell array.
        d=d{1};
    end
    if (allstrings)     % put back to strings if it was originally.
        d = string(d);
    end
end
