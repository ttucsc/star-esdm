function version_info = ARRM_V2_version_info(do_update)
%version_info = ARRM_V2_version_info(do_update)
%   returns contents of file ARRM_V2_version_info.sh
%
%   do_update       optional.  if true, runs script g-t_ARRM_V2_version_info.sh to extract info from git archive.
%                                       otherwise, simply returns current contents of file.
%
    if (exist("do_update","var") && do_update)
        if (isfile("util_general/git_ARRM_V2_version_info.sh"))
            system("util_general/git_ARRM_V2_version_info.sh");
        end
    end
    
    if (isfile("ARRM_V2_version.txt"))
        [~,version_info] = system("cat ARRM_V2_version.txt");
        finf = dir("ARRM_V2_version.txt");
        version_info = sprintf("%s\n(status from:  %s as of %s)", version_info, finf.name, finf.date);
        version_info = sprintf("%s\n(currently:                              %s)", version_info, datestr(now));
    else
%       status_file="ARRM_V2.m";
        finf1  = dir("ARRM_V2.m");
        finf2 = dir("ARRM_V2_wrapper.m");
        version_info = sprintf("ARRM_V2 (version unknown;  ARRM_V2_version.txt is missing)\n%-25s %s\n%-25s %s", finf1.name, finf1.date, finf2.name, finf2.date);
    end
    
    matlab_version = version();
    version_info = sprintf("%s\n\nmatlab version: %s\n", version_info, matlab_version);
end

