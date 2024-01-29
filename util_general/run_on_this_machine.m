function [yesno, run_on_hpcc, run_on_hostname, this_hostname] = run_on_this_machine(run_on_hostname, this_hostname)
%
%   returns true if a "match" between this_hostname and run_on_hostname
%   where "match" includes matching across RedRaider, Quanah, and quanah and nocona nodes
%
%   calls system to get this_hostname if not provided.
%   
%   Only looks at 1st part of hostname, up to first period.  (i.e., unqualified hostname). 
%   So it matches "neys.ttu.edu" with "neys".
%
%   Also returns whether or not run_on_hostname is one of the HPCC's
%   systems, so calling code knows whether to set up slurm stuff or a
%   simpler bash script.
%
%   Also returns unqualified values of this_hostname and run_on_hostname.

    if (~exist("this_hostname","var"))
        [~,this_hostname] = system("hostname");
    end

    this_h = strsplit(this_hostname,".");
    runon_h = strsplit(run_on_hostname,".");
    
    hpcc_names = ["login","quanah","cpu-","nocona"];
    run_on_hpcc = false;
    if (strncmp(run_on_hostname,"cpu-",4))
        run_on_hpcc = true;
    else
        for i=1:length(hpcc_names)
            l=strlength(hpcc_names(i));
            if (strncmpi(run_on_hostname, hpcc_names(i), l))
                run_on_hpcc = true;
                break;
            end
        end
    end
    
    
    this_hostname = this_h(1);
    run_on_hostname = runon_h(1);
    
    if (strcmpi(this_hostname, run_on_hostname))
        yesno = true;
    elseif (strcmp(this_hostname,"neys"))
        yesno = strcmp(this_hostname, run_on_hostname);
    elseif (strncmp(this_hostname,"cpu-",4))
        yesno = any(strcmp(run_on_hostname,["nocona","quanah","xlquanah"]));
    elseif (any(strcmpi(this_hostname,["icsf-lmac","icsf-kmac"])))
        yesno = any(strcmpi(this_hostname, run_on_hostname));
    else
        yesno = false;
    end
end
