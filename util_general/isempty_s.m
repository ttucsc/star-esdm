function tf = isempty_s(x)
% extension of isempty() which returns true if x is a string with zero length.
% Matlab's isempty(...) unfortunately returns false if s is a string object but has no chars in it.
     tf = (isstring(x) && length(x) == 1 && strlength(x)==0) || isempty(x);
end 