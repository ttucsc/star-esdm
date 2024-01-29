function [word, rest] = split_word(l)
% returns 1st word from l, and the rest of the contents.
% word has leading and trailing spaces removed.
% rest has leading spaces removed.
% To split a string into many tokens with the same delimiter, use split(...).

%     l = trim_leading(l);
    if (isstring(l))
        l = char(l);
        was_string = true;
    else
        was_string = false;
    end
    p = find(~isspace(l),1);        % code from trim_leading(...), so split_word(...) is standalone.
    if (isempty(p))
        l='';
    else
        l=l(p:end);
    end

    p = find(isspace(l),1);
    if (isempty(p))
        word=l;
        rest='';
    else
        word=l(1:p-1);
        rest=trim_leading(l(p+1:end));
    end
    if (was_string)
        word = string(word);
        rest = string(rest);
    end
end
        