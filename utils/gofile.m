function gofile(fname, linenum, varargin)

    fl=strsplit(fname,':');
    if (length(fl)>1)
        fname=fl{1};
        linenum=floor(str2double(fl{2}));
    elseif (~exist('linenum','var'))
        linenum=1;
    end
    
    if (isempty(which(fname))), error("no such file; %s", fname); end
        
    matlab.desktop.editor.openAndGoToLine(which(fname), linenum);
end