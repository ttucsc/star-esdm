function s=fix_spaces(s)
%   replaces spaces with underscores.
%   for filenames, see function underscore(...), which replaces all non-alphanum chars with underscores.

    wasstring=false;
    if (isstring(s))
        wasstring=true;
        s=char(s);
    end    
    s(isspace(s))='_';
    if (wasstring)
        s=string(s);
    end
end