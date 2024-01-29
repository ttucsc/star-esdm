function [w, punc1, punc2] = trim_punctuation(w, puncs)
% removes trailing punctuation from w.
% puncs should be a char array of punctuation chars to strip off

    if (ischar(w))
        w = string(w);
        wasChar = true;
    else
        wasChar = false;
    end
    if (~exist('puncs','var') || isempty_s(puncs))
        puncs='(.,\n;:! \t)';
    end
    punc1="";
    if (~isempty_s(w))
        l = strlength(w);
        first = extractBefore(w,2);
        while (contains(puncs, first))
            punc1 = punc1 + first;
            w=extractAfter(w,1);
            l = strlength(w);
            if (l==0), break; end
            first = extractBefore(w,2);
        end
    end
    punc2="";
    if (~isempty_s(w))
        l = strlength(w);
        last = extractAfter(w,l-1);
        while (contains(puncs, last))
            punc2 = last + punc2;
            w=extractBefore(w,l);
            l = strlength(w);
            if (l==0), break; end
            last = extractAfter(w,l-1);
        end
    end
    if (wasChar), w=char(w); end
end
