function dirs = dirname(fnames)
%returns directory portion of a filename
    
    if (ischar(fnames) || iscell(fnames))
        waschar = true; 
    else
        waschar = false;
    end
    
    fnames=string(fnames);
    sz=size(fnames);
    dirs=strings(sz);
    for i=1:numel(fnames)
        dirs(i) = fileparts(fnames(i));
    end

    if (waschar)
        if (numel(dirs)==1)
            dirs=char(dirs);
        else
            dr2=dirs;
            dirs=cell(sz);
            for i=1:numel(dr2)
                dirs{i}=dr2{i};
            end
        end
    end
    
end

