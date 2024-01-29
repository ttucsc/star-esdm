function  upfigs(maxfig, direction)
    % pulls figures to top.  direction is 1 or -1.
    if (~exist('direction','var') || isempty(direction)), direction = 1;end
    if (~any(direction==[1,-1])), error("direction must be 1 or -1, not %s\n", string(direction)); end
    
    fignos = [];
    if (length(maxfig)==1)
        minfig = 1;
    elseif (length(maxfig)==2)
        minfig = min(maxfig);
        maxfig = max(maxfig);
    else
        fignos = maxfig;
    end
    
    if (~isempty(fignos))
        for i=fignos
            if (ishandle(i)), figure(i); end
        end
    elseif (direction == 1)
        for i=minfig:maxfig
            if (ishandle(i)), figure(i); end
        end
    else
        for i=maxfig:-1:minfig
            if (ishandle(i)), figure(i); end
        end
    end            
end

