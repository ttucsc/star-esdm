function setup_dirs(base_dir, sub_dirs, verbose)
%   setup_dirs(base_dir, sub_dirs, verbose)
%
%   Makes sure folders exist, creates them if necessary, and optionally cleans out any old .nc or .log files. 
%
%   base_dir                base folder.  Use '.' or [] for current folder;  use [] if subdirs are fully qualified.
%                               base_dir must exist.    Example:  "/data/downscaled/arrm_v2"
%   sub_dirs                single subfolder or array of subfolder names.
%                                                       Example:  ["conus","global","misc/test"]
%
%   If base_dir is "/data" and subdirs is "abc/def/ghi", will create each subdir in order: 
%                   /data/abc
%                   /data/abc/def
%                   /data/abc/def/ghi
%
%-----------------

    if (~exist('verbose','var') || isempty(verbose))
        verbose = false; 
    elseif (ischar_s(verbose) &&strcmp(verbose,"verbose"))
        verbose=true;        
    end
    
    if (~isempty_s(base_dir))
        if (~isfolder(base_dir)), error("error:  base_dir %s does not exist", base_dir); end
    end

    if (~exist('sub_dirs','var') || isempty(sub_dirs)), return; end
    
    if (~isstring(sub_dirs)), sub_dirs = string(sub_dirs); end  % in case it's a single char array or cell array of char stuff.
    
    my_subdirs = mk_subdirs(sub_dirs);
    for i=1:length(my_subdirs)        
        full_dir=fullfile(base_dir,my_subdirs(i));
        if (~isfolder(full_dir))
            if (verbose), fprintf("creating folder:        %s\n", full_dir); end
            mkdir(char(full_dir));
        elseif (verbose), fprintf("folder already exists:  %s\n", full_dir);  
        end
        if (~isfolder(full_dir)),error('cannot create folder %s', full_dir) ; end
    end
    
end

function my_subdirs = mk_subdirs(subdirs)
    my_subdirs = strings(0,0);
    for i=1:length(subdirs)
        parts=split(subdirs(i),filesep);
        for j=1:length(parts)
            subdir = join(parts(1:j),filesep);
            if (isempty_s(subdir))
                continue    % nothing to do if root folder.  Can't delete, and don't need to create.
            else
                my_subdirs(end+1)=subdir; %#ok<AGROW>
            end
        end
    end
end

    


