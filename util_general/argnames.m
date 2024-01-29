function args = argnames(varargin)
    n = length(varargin);
    args = {};
    i=1;
    while(i<=n)
        arg = varargin{i};
        if (isstruct(arg))
            args=[args;fields(arg)]; %#ok<AGROW>
            %args = targs;
        else
            args=[args; {arg}]; %#ok<AGROW>
            %args = targs;
            i=i+1;
        end
        i=i+1;
    end
end

