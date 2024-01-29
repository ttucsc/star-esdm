function pop_workspace(keep_cur)
%restores previously pushed workspace
%
%   inputs:
%       keep_cur    bool    false: clear current variables
%                           true:  keep current workspace variables as well as restore stack
%                                       note:   if a current var name is in stack top, 
%                                               stack top will overwrite current var
%
    if (~exist('keep_cur','var') || isempty(keep_cur))
        keep_cur = false;
    end
    % Pop last workspace off stack
    c = getappdata(0, 'WORKSPACE_STACK');
    if isempty(c)
        warning('Nothing on workspace stack');
        return;
    end
    s = c{end};
    c(end) = [];
    setappdata(0, 'WORKSPACE_STACK', c);

    if (~keep_cur)
    % Do this if you want a blank slate for your workspace
        evalin('caller', 'clear');
    end

    % Stick vars back in caller's workspace
    names = fieldnames(s);
    for i = 1:numel(names)
        assignin('caller', names{i}, s.(names{i}));
    end

end

