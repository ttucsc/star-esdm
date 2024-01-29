function s = underscore(s, maxlen)
% replaces spaces in s with underscores
% useful for makeing filenames from strings that might contain spaces.
% NOTE:  replaces slashes with underscores, too.

    if (~exist('maxlen','var') || isempty(maxlen)), maxlen=0; end

    if (iscell(s))
        for i=1:length(s)
          s{i} = underscore(s{i},maxlen);
        end
    elseif (isstring(s))
        if (length(s) == 1)
            s = string(underscore(char(s), maxlen));
        else
            for i=1:length(s)
                s(i) = underscore(s(i),maxlen);
            end
        end
    elseif (ischar(s))
%         ix = isspace(s);
        if (maxlen>0 && strlength(s)>maxlen), s=extractBefore(s,maxlen+1); end
        ix = ~isstrprop(s,'alphanum');
        s(ix) = '_';
    end
    
end