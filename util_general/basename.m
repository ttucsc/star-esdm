function [bnames, exts] = basename(fnames, no_ext)
%returns basename of a file
%If 2 args out or no_ext is true, returns extension separately
    if (~exist("no_ext","var") || isempty(no_ext))
        no_ext = false;
    end
    if (nargout ==2), no_ext = true; end

    if (ischar(fnames) || iscell(fnames))
        waschar = true; 
    else
        waschar = false;
    end
    
    fnames=string(fnames);
    sz=size(fnames);
    bnames=strings(sz);
    exts  =strings(sz);
    for i=1:numel(fnames)
        [~,bnames(i), exts(i)] = fileparts(fnames(i));
        if (~no_ext), bnames(i)=sprintf('%s%s',bnames(i), exts(i)); end
    end

    if (waschar)
        if (numel(bnames)==1)
            bnames=char(bnames);
            exts=char(exts);
        else
            bn2=bnames;
            ex2=exts;
            bnames=cell(sz);
            exts=cell(sz);
            for i=1:numel(bn2)
                bnames{i}=bn2{i};
                exts{i}=ex2{i};
            end
        end
    end
end

