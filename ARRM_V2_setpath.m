function [hostname, progname, progext] = ARRM_V2_setpath(progname)
%   set path for ARRM_V2 code, adding the required folders if not already in the path.
%
%   if called with no progname, uses its own location to set the path.
%   makes sure that main folder for progname, as well as subfolders "util_general", "util_nc" and "ncdf" are in path
%   returns hostname & name of calling program.
%   if progname is empty or missing, uses its own location as the base folder, so can be called from the commandline.
%
%   *** NOTE:  THIS FILE SHOULD BE IN THE BASE FOLDER FOR THE ARRM_V2 code ***
%
%   Because this needs to access the file system, don't call this unnecessarily, as it takes measurable time to check
%   that the folders exist.
%
    [~,hostname]=system('hostname');
    hostname = strtrim(hostname);

            % if no filename given, then try current folder
    if (~exist('progname','var'))
        try
            progname=mfilename('fullpath');
            [basedir, progname, progext] = fileparts(char(progname));
        catch
            progname=strings(0);
            basedir = pwd();
        end
    else
        [basedir, progname, progext] = fileparts(char(progname));
        if (isempty(basedir))
            basedir=pwd();
        end
    end
    
%   fid=fopen("whereami.txt","a");
%   fprintf(fid,"%s %s hostname: %s basedir: %s progname: %s\n", pwd, datestr(now,"yyyy-mm-dd HH:MM:SS"), hostname, basedir, progname);
    
    if (~isfile(fullfile(basedir,'ARRM_V2.m')))
        fprintf(2, 'warning:  cannot find ARRM_V2.m to set path properly');
    end
    
           % ok.  we know where ARRM_V2.m is. Add the required subfolders. 
    
    dirs=["";"util_general";"util_nc";"ncdf";"ARRM_V2_subs"; "utils"];
%   dirs=["";"util_general";"util_nc";"ncdf";"ARRM_V2_subs";"prcp"];        % use this with precip version.
    
    p=split(path,':');
    
    for i=1:length(dirs)
        dd=fullfile(basedir,dirs(i));
                %  make sure folders all exist as well.
        if (~isfolder(dd))
           fprintf(2, "error:  folder %s does not exist\n", dd);
        elseif (~(any(p==dd)) )    % if it's not already in the path, add it
            addpath(char(dd));
        end
    end
end

