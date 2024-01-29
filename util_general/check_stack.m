function [caller, main, depth, stack] = check_stack()
%   returns name of main routine (starting function), stack depth, and calling stack.
%       Outputs:
%           main        name of main function
%           depth       # of functions deep
%           stack       call stack
%
%   useful to determine environment (where current program got started)
%
%
    stack = dbstack;
    depth = length(stack);
    if (depth==1)
        caller=strings(0,0);
        main=strings(0,0);
    else
        main = stack(end).name;
        if (depth<3)
            caller = strings(0,0);
        else
            caller = stack(3).name;
        end
    end
    depth = depth-1;
    if (nargout > 3)
        stack=stack(2:end);
    end
end

