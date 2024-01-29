function tf = ischar_s(s)
% returns true if s is a one-row char array or s is an empty or single string
% NOTE:  To test if s is a multi-row array of chars (all the same length) or an array of strings, use ischars(...).
    tf = (ischar(s) && size(s,1)<=1) || (isstring(s) && length(s)<=1); 
end

