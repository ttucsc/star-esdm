function len = nonwhitelength(str)
%returns length of string str, excluding whitechars.

    if (isstring(str) || iscell(str))
        len=zeros(length(str),1);
        for i=1:length(str)
            len(i)=nonwhitelength(str{i});
        end
    else
        len = sum(~isspace(str));
    end
end

