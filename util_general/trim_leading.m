function l = trim_leading(l)
%   removes leading blanks from l.  (does not remove trailing blanks.)
%   To trim both leading and trailing, use strtrim(...);
%   To remove trailing only, use deblank(...)
    if (isstring(l))
        was_string = true;
        l = char(l);
    else
        was_string = false;
    end
    p = find(~isspace(l),1);
    if (isempty(p))
        l='';
    else
        l=l(p:end);
    end
    if (was_string)
        l = string(l);
    end
end

