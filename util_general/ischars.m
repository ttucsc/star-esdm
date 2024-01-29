function tf = ischars(s)
% returns true if s is a char array or s is a string or string array or cell array of char arrays.
%
%   Like ischar_s(...), except also returns true for array of strings or for cell array or char stuff.
    if (iscell(s))
        for i=1:length(s)
            if (~ischars(s{i}))
                tf=false;
                return;
            end
        end
        tf=true;
    else
        tf = ischar(s) || isstring(s); 
    end
end

