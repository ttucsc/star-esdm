function [mkeepers, unmatched] = find_matches(names, matchlist, insensitive)
% returns a bool array showing which strings in matchlist are in array "names".
%
%   insensitive, if present and true, does case-insensitive matching.
%   
%   if matchlist contains the string "-" , then list returned is inverted (not-matches)
%   if a matchlist entry starts with "~", then matching is partial:  "~gfdl" will match all names containing "gfdl"
%   if matchlist contains the string "~" , then all matching is partial...

    matchlist = string(matchlist);
    names = string(names);
    nnames = length(names);

    unmatched = strings(0,0);
    partial = false;
    mkeep = true;
    ix = find(matchlist=="-");
    if (~isempty(ix)),mkeep = false; matchlist(ix) = []; end   % exclude models, rather than include.
    ix = find(matchlist=="~");
    if (~isempty(ix)),partial = true; matchlist(ix) = []; end   % do partial matching.
    
    origlist = matchlist;
    if (insensitive)
        names = upper(names);
        matchlist = upper(matchlist);
    end
    if (~isempty(matchlist))
        matches = false(nnames,1);
        for i=1:length(matchlist)
            myword = matchlist(i);
            mypartial = false;
            if (extractBefore(myword,2) == '~')
                mypartial = true;
                myword = extractAfter(myword,1);
            end
            if (partial || mypartial)
                matched = contains(names, myword);
            else
                matched = strcmp(names, myword);
            end
                
            if (~any(matched))
                unmatched = cat(1,unmatched,origlist(i));
            else
                matches = matches | matched;
            end
        end
        if (mkeep)
            mkeepers = matches;
        else
            mkeepers = ~matches;
        end
    else
        mkeepers = true(nnames,1);
    end
end
