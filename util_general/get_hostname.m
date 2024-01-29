function [hostname, on] = get_hostname(name_only)
%[hostname, on] = get_hostname(name_only)
%
%   return hostname of system we're running on.
%   if name_only is present, then the hostname is trimmed to just the first part of the name
%   (i.e., "neys.ttu.edu" gets truncated to "neys")

    [~,hostname]=system('hostname');
    hostname = string(strtrim(hostname));
    
    if (exist("name_only","var") && name_only)
        hh = strsplit(hostname,".");
        hostname = hh(1);
    end
    
    if (nargout > 1)
        on.hostname    = hostname;
        on.neys        = contains(hostname,"neys");
        on.icsf_kmac   = contains(hostname, "kmac");
        on.icsf_lmac   = contains(hostname, "lmac");
        on.office_mac  = on.neys | on.icsf_kmac;
        on.laptop      = contains(hostname,"jcsf") || contains(hostname,"jpro") || contains(hostname, "macair");
        on.dev_system  = contains(hostname,"icsf");
        on.hpcc_system = any(contains(hostname,["quanah","login","cpu","hrothgar"]));   % hrothgar is probably gone now...
        on.KMQ         = contains(hostname,"mara") || contains(hostname,"quetico") || contains(hostname,"killarney");
    end
end
    
