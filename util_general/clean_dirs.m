function clean_dirs(base_dir, sub_dirs, varargin)  % options:   "extensions", exts_to_delete, "clean_intermediates", clean_intermediates, "verbose", verbos
%
%   clean_dirs(base_dir, sub_dirs, "extensions", exts_to_delete, "clean_intermediates", clean_intermediates, "verbose", verbose)
%   Cleans out any old .nc or .log files (or other extensions, if specified) from specified subfolders of base_dir. 
%
%   base_dir                base folder.  Use '.' or [] for current folder;  use [] if subdirs are fully qualified.
%                               NOTE:   Normally, base_dir is NOT cleaned, unless subdirs variable is missing or empty.
%                                       If subdirs is empty or missing, base_dir will be cleaned.
%                                       to clean basedir AND subdirs, put an empty string in subdirs.
%   sub_dirs                single subfolder name or string array of subfolder names.
%                                       Example:  ["conus","global","misc/test"]  to clean 3 sudfolders
%                                       Example:  ["", "conus, "global, "misc/test"] to clean basedir and 3 subfolders 
%   Optional keyword/value pairs (can be substring, such as "ext" or "clean":
%       "extensions", exts_to_delete          optional list of extensions to delete.  default:  [".nc",".log"]
%                                                   set to empty to delete .nc and .log files.  
%                                                   if not empty, replaces default (i.e., does not add your extensions
%                                                   to defaul list.)
%       "clean_intermediates", t/f     boolean;  if false, deletes only in specified folders.  default: false
%                                                if true, cleans intermediate folders as well.
%                                                   i.e., if "misc/test" is in the list of subdirs, and
%                                                   clean_intermediates is true, then cleans misc & misc/test.
%                                      if false (default), will clean only misc/test
%       "verbose", t/f                boolean.   outputs name of each file deleted.
%
%-----------------

    if (~exist("base_dir","var") || isempty_s(base_dir)), base_dir = ""; end
    if (~exist("sub_dirs","var") || isempty_s(sub_dirs)), sub_dirs = ""; end

                    % parse input for DA_title
    p = inputParser;        
    addParameter(p,"extensions",[".log",".nc"],  @(s) ischar_s(s));
    addParameter(p,"verbose",             false, @(s) islogical(s) || (isnumeric(s) && any(s==[0,1])));
    addParameter(p,"clean_intermediates", false, @(s) islogical(s) || (isnumeric(s) && any(s==[0,1])));
        
    parse(p, varargin{:});

    exts_to_delete      = p.Results.extensions;
    verbose             = p.Results.verbose;
    clean_intermediates = p.Results.clean_intermediates;
    
    if (isempty_s(base_dir) && isempty_s(sub_dirs))
        error("error:  cowardly refusing to clean anything because no folder specified.  To clean current folder, please specify ""."" ");
    end
    
    if (~isempty_s(base_dir))
        if (~isfolder(base_dir)), error("error:  base_dir %s does not exist", base_dir); end
    end
    
    if (~isstring(sub_dirs)), sub_dirs = string(sub_dirs); end  % in case it's a single char array or cell array of char stuff.
    
    for i=1:length(sub_dirs)
        if (~clean_intermediates)
            my_subdirs = sub_dirs(i);
        else
            my_subdirs = make_subdirlist(sub_dirs(i));
        end
        for j=1:length(my_subdirs)
            dirname = fullfile(base_dir, my_subdirs(j));
            if (~isfolder(dirname))
                fprintf(2,"warning:  no such folder:  %s\n", dirname);
            else
                clean_files(dirname, exts_to_delete, verbose);
            end
        end
    end    
end

function clean_files(dirname, exts_to_delete, verbose)

    finfo = dir(dirname);
    for i=1:length(finfo)
        fn = finfo(i).name;
        c = extractBefore(fn,2);        % look for  ., .. and any hidden files.
        if (c=='.'), continue; end      % skip  ., .. and any hidden files.
        fname = fullfile(dirname, fn);
%       fprintf('fn: %s   c: %c\n', fname, c);
        if (isfile(fname))
            [~,~,ext] = fileparts(fname);          
            if (any(strcmp(ext, exts_to_delete)))
                delete(fname);
                if (verbose)
                    fprintf('deleting %s\n', fname);
                end
            end
        end
    end
end

function my_subdirs = make_subdirlist(subdir)

    if (isempty_s(subdir))
        my_subdirs = "";
        return;
    end
    c=extractBefore(subdir,2);
    fully_qualified = c==filesep;
        
    parts=split(subdir,filesep);
    parts(strlength(parts)==0) = [];   % remove any empty strings.  Can occur if subdir has //
    my_subdirs = strings(size(parts));
    for j=1:length(parts)
        subdir = join(parts(1:j),filesep);
        if (fully_qualified), subdir = strcat(filesep, subdir); end
        my_subdirs(j)=subdir;
    end
end

    


