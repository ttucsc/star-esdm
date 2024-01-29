function  movefigs(figrange, newpos, absolute)
    % pulls figures to top, and optionally shifts them by (dx,dy) pixels.

    if (~exist("absolute","var") || isempty(absolute))
        absolute = false;
    elseif (isnumeric(absolute))
        absolute = logical(absolute);
    elseif (strncmpi(string(absolute),"abs",3))
        absolute = true;
    else
        absolute = false;
    end

    if (~exist("newpos","var") || isempty(newpos))
        newpos = [0,0,0,0];
    else
        newpos(end+1:4)=0;
    end
        
    
    fignos = [];
    if (length(figrange)==1)
        fignos = 1:figrange;
    elseif (length(figrange)==2)
        if (figrange(2)>=figrange(1))
            fignos = figrange(1):figrange(2);
        else
            fignos = figrange(1):-1:figrange(2);
        end
    elseif (length(figrange)>2)
        fignos = figrange;
    end    
    
    for i=fignos
        if (ishandle(i))
            h = figure(i);
            if (~all(newpos==0))
                uu = h.Units;
                if (any(mod(newpos,1) ~= 0))
                    h.Units = "normalized";
                else
                    h.Units = "pixels";
                end
                if (~absolute)
                    h.Position = h.Position + newpos;
                else
                    h.Position = newpos;
                end
                h.Units = uu;
            end
        end
    end
end

