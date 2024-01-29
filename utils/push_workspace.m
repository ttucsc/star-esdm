function push_workspace()
%saves all current workspace variables on a stack.  
%
%   Careful using this if you're got large arrays in memory!
%

    c = getappdata(0, 'WORKSPACE_STACK');
    if isempty(c)
        c = {};
    end

    % Grab workspace
    w = evalin('caller', 'whos');
    names = {w.name};
    s = struct;
    for i = 1:numel(w)
        s.(names{i}) = evalin('caller', names{i});
    end

    % Push it on the stack
    c{end+1} = s;
    setappdata(0, 'WORKSPACE_STACK', c);


end

