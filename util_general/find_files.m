function [fnames, not_found, finfos] = find_files(fns, basedir, findargs)
% [fnames, not_found, finfos] = find_files(fns, basedir, findargs)
% uses linux find to locate files in path basedir
%
%   Inputs:
%       fns         string array of files to look for.  May contain any usual linux wildcards and other search stuff
%       basedir     base folder to start from.  default:  "."
%       findargs    optional string variable with any additional arguments to find, such as "-newer somefile"
%                       if empty, search just does   "find basedir -name fns(i) -print" on each filename in fns 
%
%   Outputs
%       fnames      string array with just the fully qualified paths to each file
%       not_found   string array of files not found
%       finfos      array of finfo struct.  see matlab function dir(...)
%
%   NOTE:  finfos and fnames will not contain entries for files in not_found.  They will be empty if no files are found.

    if (~exist("basedir","var")  || isempty(basedir)),  basedir  = "."; end
    if (~exist("findargs","var") || isempty(findargs)), findargs = "";  end
    
    fnames    = strings(0);
    not_found = strings(0);
    for i=1:length(fns)
        [st, found] = system(sprintf("find %s %s -name ""%s"" -print", basedir, findargs, fns(i)));
        if (st ~= 0 || isempty(found) || strlength(found)==0)
            not_found = [not_found; fns(i)]; %#ok<AGROW>
        else
            found = splitlines(found);
            fnames = [fnames; found]; %#ok<AGROW>
        end
    end
    
        % remove any empty lines
    fnames(strlength(fnames)==0) = [];
    
    nfiles = length(fnames);
    if (nargout > 2 && nfiles>0)
        finfos = dir(fnames(1));
        for i=2:nfiles
            finfs = dir(fnames(i));
            finfos = [finfos; finfs]; %#ok<AGROW>
        end
    else
        finfos = [];
    end
end
