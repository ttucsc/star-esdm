function s = fix_filename(s, no_path)
% replaces spaces in s with underscores
% useful for makeing filenames from strings that might contain spaces.
%   no_path:  set to true if filename contains no path info.
%             this allows replacing filesep (/ or \) with underscores.
%

    if (~exist("no_path","var") || ~no_path)
        s = strsplit(string(s),filesep);
    end

    if (length(s) == 1)
        cs = char(s);
        ix = ~isstrprop(cs,'alphanum');
        cs(ix) = '_';
        s = string(cs);
    else
        for i=1:length(s)
            s(i) = fix_filename(s(i));
        end
        s = join(s,filesep);
    end    
end