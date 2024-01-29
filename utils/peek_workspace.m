function peek_workspace(names)
%lists variables in top entry of workspace stack, or reports whether names are in top entry of workspace stack
%
    if (~exist('names','var') || isempty(names))
        peek_all = true;
    else
        peek_all = false;
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
    if (peek_all)
        fprintf('stack top peek:\n');
        names=fieldnames(s);
        for i = 1:numel(names)
            fprintf('\t%20s: %s\n', names{i}, my_info(s.(names{i})));
        end
    else
        fprintf('stack top peek:\n');
        fldnames = fieldnames(s);
        for i = 1:numel(names)
            if (find(strcmp(fldnames,names{i})))
                fprintf('\t%20s: %s\n', names{i}, my_info(s.(names{i})));
            else
                fprintf('\t%20s: (missing)\n', names{i});
            end
        end
    end
end
function info = my_info(myvar)
    sz = size(myvar);
    info = [ '[',sprintf('%d ', sz), sprintf('%s',class(myvar)),']'];
end
