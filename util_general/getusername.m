function [uname, longname, hostname,homedir] = getusername()
% returns username, user's long name, hostname and home directory from system
% This needs to be updated to support Windows systems as well.  At present
% will only work on unix-like (linux, macos, etc.)
% using dscacheutil to get user's longname only works on MacOS...
    [~,      hostname] = system('hostname');
    [~,      uname]    = system('echo $USER');
    [~,      homedir]  = system('echo $HOME');

    hostname = string(strtrim(hostname));
    uname    = string(strtrim(uname));
    homedir  = string(strtrim(homedir));
    longname="";

    [ls,typ] = system("echo $OSTYPE");  % on MacOS, getenv("OSTYPE") fails, so use system(echo...) to get it.
    if (ls==0)
        if (contains(typ,"darwin"))
            [lstatus,longname] = system("dscacheutil -q user -a name $(whoami) | fgrep gecos | sed -e 's/.*gecos: \(.*\)/\1/' ");
            if (lstatus ~= 0), longname=""; end
        elseif (contains(typ,"linux"))
            [lstatus,pwd_info] = system(sprintf("grep %s /etc/passwd", uname));
            if (lstatus==0)
                p=string(strsplit(pwd_info,":"));
                if (length(p)>=5)
                    longname=p(5);
                end
            end
        end   
    end
    
    if (strlength(longname)>0)
        longname = string(strtrim(longname));
    else
        longname = uname;
    end
end

