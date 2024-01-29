function NC_attributes = read_global_attributes(fname, NC_attributes)
%   Reads fname and returns struct w/ keyword/values. 
%   All values returned as strings (note: not char arrays), even if they contain numeric-only data.
%
%   Lines beginning with # are treated as comments.
%   first word is treated as the keyword, rest of line is treated as the value.
%   Multi-line values can be created using single or double quotes.  Embedded single-quotes must be quoted outside with
%   double quotes, and vice-versa.
%   Quoted strings have leading and trailing spaces (outside the quotes) removed.  
%   Non-quoted strings retain trailing spaces.

    if (iscell(fname))
        if (mod(length(fname),2) ~= 0)
            error('error:  attributes must be cell array of keyword,values pairs');
        end
        for i=1:2:length(fname)
            NC_attributes.(fname{i})=fname{i+1};
        end
        return; 
    elseif (isstruct(fname))
        kwds=fieldnames(fname);
        for i=1:length(kwds)
            NC_attributes.(kwds{i}) = fname.(kwds{i});
        end
        return;
    end
    fid = fopen(fname);
    if(fid == -1), error("error: %s:  cannot open attributes file %s for reading", mfilename, fname); end

    nl=0;
    while true    
        tline=fgetl(fid);
%         fprintf('line read: ''%s''\n', tline);
        if (tline == -1), break; end
        nl=nl+1;
        [keyword,val] = split_word(tline);
        if (isempty_s(keyword) || extractBefore(keyword,2) == '#'), continue; end
        if (isempty_s(val)), error('read_global_attributes: no value for keyword in file %s at line %d: ''%s''\n', fname, nl, tline); end
        if (val(1) == '"' || val(1) == '''')
            ll=deblank(val);
            if (ll(end)==ll(1))
                val = ll(2:end-1);
            else
                [rest, nl] = rest_of_quote(fid,val(1), nl);
                v = [val, rest];
                val = v(2:end-1);    % gets rid of warning...
            end
        end
                % either append to existing value (so we add to comments, rather than overwrite), or create new if not
                % keyword present.
        fields = fieldnames(NC_attributes);
        if (any(strcmp(keyword,fields)) && ischar_s(val) && ischar_s(NC_attributes.(keyword)))
            NC_attributes.(keyword) = sprintf('%s\n%s', NC_attributes.(keyword), val);
        else
            NC_attributes.(keyword) = string(val);
        end
%         fprintf('keyword: %s  val: ''%s''\n', keyword, NC_attributes.(keyword));
    end
end
function [rest, nl] = rest_of_quote(fid, delim, nl)
    rest = '';
    n=0;
    while true
        tline = fgetl(fid);
%        fprintf('in roq:  line is ''%s''\n', tline);
        if (tline == -1), error('read_global_attributes: unterminated quote starting at line %d\n', nl); end
        n=n+1;
        ll = deblank(tline);
        if (~isempty_s(ll) && ll(end) == delim)
            r = [rest,newline,ll];
            rest = r;   % gets rid of warning...
            nl=nl+n;
            return;
        end
        r = [rest,newline,tline];
        rest = r;
%         fprintf('rest: ''%s''\n', rest);
    end
end

