function [fnstruct, fullpath]  = dir_canonical(fnames)
% function [fnstruct, fullpath]  = dir_canonical(fnames)
%
%   Like dir(...), except follows symlinks and returns  canonical paths.
%   One difference from dir:  if the file does not exist, the fnstruct is not empty, but the filename, folder, etc.
%   are empty strings ( -1 for file bytes and datenum).  
%
%       If you just need the full canonical path to a file, use canonical_path(filename)
%
%   returns file info struct and full path to actual file(s) for fnames, following symlinks as needed.  
%       If a filename is a file, returns dir(...) struct for filename.
%       If fnames is a symbolic link, or contains a symbolic path, returns dir(...) struct for actual file for fname
%       If actual file does not exist, filestruct is empty
%
%   NOTE:  this will be slower than a call to dir(fname) because we have to use a loop and
%           run a python script to get canonical names.
%
%   Inputs:
%       fname       string, array of strings, char array or cell array of chars of filenames
%                       Can also be an array of dir(...) output structs
%                       To get the canonical directory info for all files in current folder, leave fname off
%                       To get the canonical directory info for all files in another folder, use
%                           [[fnstructs, fullpaths] = dir_canonical(dir(somedir));
%
%   Outputs
%       filestruct  file info struct or array of structs from dir(...), 1 for each filename in fnames
%       fullpath    string array of full canonical paths to filenames.
%                       Test for 
%                           strlength(fnstruct.folder) == 0 
%                       to determine if path is valid.
%

            % first, figure out what the input looks like, and make an array of strings from it.
            
    if (~exist('fnames','var') || isempty(fnames))  % no input.  get directory of current folder.
        fnames = dir();  % NOTE:  dir() will contain . and .. 
    elseif (isstruct(fnames))
        if (~isfield(fnames(1),'folder')), error('fnames is struct, but not dir() struct'); end
    elseif (ischar(fnames))         % single char array
        fnames=string(fnames);
    elseif (iscell(fnames))         % cell array, which should just be char arrays.
        fnames=string(fnames); 
    end
    nfns = numel(fnames);
    fnsize = size(fnames);
    
            % create empty arrays for output
    fnstruct = repmat(struct('name',"", 'folder', "", 'date',"", 'bytes',-1, 'isdir',false, 'datenum',-1), fnsize);  
    
    fullpath=strings(fnsize);

    for i=1:nfns
        [fnstr,fullp] = canonical_sub(fnames(i));
        if (~isempty(fnstr))
            fnstruct(i) = fnstr;
            fullpath(i) = fullp;
        end
    end
end

function [fnstruct, fullpath] = canonical_sub(fname)

% readlink doesn't return canonical path on MacOS unless file is in cwd, so we use a python script.
%   python script is called in wrapper function canonical_path.m

    if ( isstruct(fname))    % could have passed in a struct from an earlier call of dir();
        fnstruct = fname;
    elseif (isfolder(fname))
        fnstruct = my_java_filestruct(fname, true);
    else
        fnstruct = dir(char(fname));
    end
    
    if(isempty(fnstruct))
        fullpath="";
        return;
    end
    
        % if folder is relative, make it absolute
    if (strlength(fnstruct.folder)==0)
        fnstruct.folder=pwd;
    else
        if (extractBefore(fnstruct.folder,2) ~= filesep)
            fnstruct.folder = fullfile(pwd, fnstruct.folder);
        end
    end
    
    fullpath = fullfile(fnstruct.folder, fnstruct.name);

    cname = canonical_path(fullpath);
    
    if (~strcmp(fname, cname))
        [fp, ~, ~] = fileparts(cname);
        fnstruct.folder = fp;
        fullpath = cname;
    end

end

function fnstruct =  my_java_filestruct(fname, isdir)
%   Uses java.io.File to construct a file struct for a directory.
%   Using dir(...) on a directory returns its contents, not the info about a file.
        
    f=java.io.File(fname);
    modtime=f.lastModified()/1000/86400;        % days since 1/1/1970 :  msecs/1000 / #secs/day.
    [~,tzs]=system('date +"%z"');
    tz=str2double(tzs)/100/24;
    dnum=datenum([1970,1,1]) + tz + modtime;
    dt=datestr(dnum);

    [fp,fb,fe] = fileparts(fname);
    
    fnstruct=struct('name',sprintf("%s%s",fb,fe), 'folder', fp, 'date',dt, 'bytes',f.length(), 'isdir',isdir, 'datenum',dnum);
end
