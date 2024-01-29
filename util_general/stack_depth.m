function depth = stack_depth()
%   returns call-stack depth, excluding call to stack_depth() itself.
%   if stack_depth() is called from the console, returns 0
%   if called from main function, returns 1
%   etc.
%   if stack_depth() > 1, then you're in a child function.
%
%   See also:  [caller, main, depth, stack] = check_stack()
%
    stack=dbstack();
    depth=max(0,length(stack)-1);
end

