function show_progress(nloops, maxloops, barlen)
% show_progress(nloops, maxloops, barlen)
%   outputs the next char of a progress bar if it is time to do so.
%
%   Inputs:
%       nloops      current loop counter
%       maxloops    max loop count expected
%       barlen      length of bar to be output.  default:  100
%
%   This accomplishes the same as matlab's waitbar(...), but works in text mode so doesn't require a matlab GUI.
%
%   Normally outputs something along the lines of:
%
%....+....1....2....+....3....+....4....+....5....+....6....+....7....+....8....+....9...+....*
%
%   but results will be a little squirrelly if maxloops is less than barlen.
%
    if (~exist('barlen','var') || isempty(barlen)), barlen = min(100, maxloops); end
    if (nloops > maxloops)
        nloops = mod(nloops-1, maxloops)+1;
    end
    if (nloops==maxloops)
        fprintf("*\n");
    else
        delta = .5/maxloops;           % tolerance to be right on a boundary.
        outpos = floor(barlen/maxloops * nloops);      % current  location along the bar, in chars
        prvpos = floor(barlen/maxloops * (nloops-1));  % previous location along the bar, in chars

        if (outpos ~= prvpos)   % if true, time to output a char.
            tenpos = outpos/10;
            del = mod(tenpos,1);
            if (del > .5), del = 1-del; end
            if (del < delta)
                fprintf('%d', round(tenpos));
            else
                fivepos = outpos/5;
                del = mod(fivepos,1);
                if (del > .5), del = 1-del; end
                if (del < delta && barlen >= 60)
                    fprintf('+')
                else
                    fprintf('.');
                end
            end
        end
    end
end
