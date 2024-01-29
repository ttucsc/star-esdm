function pull_workspace(names)
%   pulls variables in cell array names{} into caller's workspace without 
%   popping the workspace stack
%
%   pulled variable will be a local copy of the stack's variable, 
%   so modifying it will leave the stack's variable untouched.
%
    if (~exist('names','var') || isempty(names))
        pull_all = true;
    else
        pull_all = false;
%           if varnames is not a cell array, then user gave us just 1 varname as a string.  make it a cell array.
        if (~iscell(names))
            names = {names};
        end
    end
    
    % Peek last workspace on stack
    c = getappdata(0, 'WORKSPACE_STACK');
    if isempty(c)
        warning('Nothing on workspace stack');
        return;
    end
    s = c{end};

    % Stick vars back in caller's workspace
    if (pull_all)
        names = fieldnames(s);
    end
    for i = 1:numel(names)
        assignin('caller', names{i}, s.(names{i}));
    end

end

