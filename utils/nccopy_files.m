function   nccopy_files(fns, basedir, outbase, varargin)
% function nccopy_files(fns, basedir, outbase, ChunkSize, ChunkCache, DeflateLevel, maxmem)

% Inputs:
%       fns             string array of filenames or name of file containing filenames
%                           filenames may contain "*" for wildcards.
%       basedir         basedir of where to find files.  Will be prepended to each input filename
%                           set to strings(0) for system default.
%                           set to "" if filenames contain full path info
%                           set to "." if filenames contain relative path info.
%       outbase         base destination directory default:  [/lustre/research/hayhoe/cmip6_tllc/SCENARIO]
%
%     Optional key/value pairs:
%
%       chunkesize      (optional) ChunkSize to use on output for each file.  3 values:  [nlons, nlats, ndays].
%                               Default if empty or missing:  [16,8,128].
%       ChunkCache       (optional) vector with ChunkCache to use for each file.  3 values:  [nbytes, nelems, premp].  
%                               nbytes:  # of bytes to use for cache
%                               nelems:  # of chunks to use in cache.  S/b prime #, or >= total # of chunks for file
%                               premp:   Use .25   .  See "doc netcdf.setChunkCache" for more info
%                               Default if empty or missing:  Will do its' best to fit into maxmem, based on each file's 
%                                                             current ChunkSize and specified output ChunkSize
%       DeflateLevel    (optional) compression setting.  Default if empty of missing:  0
%       maxmem          (optional) max mem available for chunking.  Probably should be < total mem - 3 GB.
%                               default if empty or missing:  100 GB on Neys, 
%                                                              60 GB on HPCC node.  *** ASSUMES 64 GB requested!!!
%                               if < 256, is assumed to be Gbytes.
%
%   For default values of optional params, see code in initparams(...).

    if (~exist("basedir","var")), basedir = strings(0); end
    if (~exist("outbase","var")), outbase = strings(0); end

    [fns, basedir, outbase, ChunkSize, ChunkCache, DeflateLevel, shuffle, add_scenario, maxmem] = initparams(fns, basedir, outbase, varargin{:});
        
    fprintf("\n----------------ncccopy_files  %s --------------------\n\n", datestr(now,"yyyy-mm-dd HH:MM:SS"));
    
    nfiles = length(fns);
    

    if (nfiles > 1)
        fprintf("copying %d files:\n", nfiles);        
    end
    fprintf("copying from basedir: %s\n", basedir);
    fprintf("writing to outbase:   %s\n", outbase);
    fprintf("using ChunkCache:    [%12d, %12d, %12d]\n", ChunkCache);
    fprintf("using ChunkSize:     %12s\n", vec2string(ChunkSize, "brackets",'[]'));
    fprintf("DeflateLevel:        %12d\n", DeflateLevel);
    fprintf("maxmem:              %12.1f Gbytes\n", maxmem/1024/1024/1024);
    
    if (~isempty(ChunkCache))
        netcdf.setChunkCache(ChunkCache(1),ChunkCache(2),ChunkCache(3));
    end
            % special case:  if outbase has "/SCENARIO at the end, then extract scenario and replace it in outbase.
    for i=1:length(fns)
        [~, fn, fext] = fileparts(fns(i));
        if (add_scenario)
            fparts = split(fn,"_");
            scenario = fparts(4);
            outdir = fullfile(outbase,scenario);
            if (~isfolder(outdir)), mkdir(outdir); end
        else
            outdir = outbase;
        end
        
        outname = sprintf("%s.tllc%s", fn, fext);
        
        fprintf("copying & rechunking %d of %d:  %s to %s\n", i, length(fns), basename(fns(i)), basename(outname));

        compress_netcdf(fns(i),outname,"DeflateLevel",DeflateLevel,"ChunkSize",ChunkSize,"verbose",true,"overwrite", true,"shuffle",shuffle,"outdir",outdir,"maxmem",maxmem/3);        
    end
    

end

function [fns, basedir, outbase, ChunkSize, ChunkCache, DeflateLevel, shuffle, add_scenario, maxmem] = initparams(fns, varargin)
    

                    % parse input for DA_title
    p = inputParser;
    p.StructExpand = true;
    p.KeepUnmatched = true;        
    
%     dv = datevec(now);
%     curyear = dv(1);
         
                    % these are the params we want to handle in ARRM_V2_wrapper
    addRequired(p,"fns",                     @(s)ischars(s));
    addRequired(p,"basedir",                 @(s)ischars(s));
    addRequired(p,"outbase",                 @(s)ischars(s));
    addParameter(p,"ChunkSize", [16,8,128],  @(s) isempty(s) ||  (isnumeric(s) && length(s)==3));       % nlons, nlats, ndays
    addParameter(p,"ChunkCache",      [],    @(s) isempty(s) ||  (isnumeric(s) && length(s)==3));       % nbytes, nelems, premp
    addParameter(p,"DeflateLevel",     0,    @(s) isnumeric(s) && length(s)==1 && (s>=0 && s <= 9));
    addParameter(p,"maxmem",          [],    @(s) isnumeric(s));                                        % max mem, in GBytes

    

    parse(p, fns, varargin{:});

    fns         = p.Results.fns;
    basedir     = p.Results.basedir;
    outbase     = p.Results.outbase;
    ChunkSize   = p.Results.ChunkSize;
    ChunkCache  = p.Results.ChunkCache;
    DeflateLevel= p.Results.DeflateLevel;
    maxmem      = p.Results.maxmem;
    
    if (isempty(ChunkSize))
        ChunkSize = [16,8,128];
    end
    
        % turning on shuffle improves the compression. 
    if (DeflateLevel > 0)
        shuffle = true;
    else
        shuffle = false;
    end
    
    [hname, on] = get_hostname();
    
    
    if (on.hpcc_system)
        maxgigs=500;
        if (isempty(basedir)), basedir = "/lustre/research/hayhoe/hamed"; end
        if (isempty(outbase)), outbase = "/lustre/research/hayhoe/cmip6_tllc/SCENARIO"; end
        if (isempty(maxmem)), maxmem = 120*1024*1024*1024; end
            
    elseif (on.neys || on.kmac)
        maxgigs=120 - 30*on.kmac;
        if (isempty(basedir)), basedir = "/Volumes/lacie_1/data/gcm_cmip6"; end
        if (isempty(outbase)), outbase = "/Volumes/lacie_1/data/gcm_cmip6/cmip6_tllc/SCENARIO"; end        
        if (isempty(maxmem))
            maxmem = maxgigs*1024*1024*1024; 
        end
    elseif (isempty(outbase) || isempty(maxmem))
        error("not set up yet to do this on %s. Specify basedir,outbase and maxmem, or call Ian", hname);
    end
    
    if (maxmem <= maxgigs)
        maxmem = maxmem * 1024*1024*1024;
    elseif (maxmem > maxgigs*1024*1024*1024)
        error("maxmem (%.1f GB) too large.  Max for this system is %d GB", maxmem/1024/1024/1024, maxgigs);
    end
    
    [outb, scen] = fileparts(outbase);
    if (strcmp(scen,"SCENARIO"))
        outbase = outb;
        add_scenario = true;
    else
        add_scenario = false;
    end
        
    if (isempty(ChunkCache))
        csize  = maxmem;
        nelems = 204803;        % largest gridsize for CMIP6 is 204800.  204803 is 1st prime > 204800.
        premp  =    .25;        % should never need to pre-empt anything from memory as long as maxmem > total memory needed to store entire file in memory.
        netcdf.setChunkCache(csize, nelems, premp);
    elseif (length(ChunkCache)==3)
        netcdf.setChunkCache(ChunkCache(1), ChunkCache(2), ChunkCache(3));
    else
        error("error:  ChunkCache needs 3 values:  csize (#bytes to use), nelems (# of elements in cache's hash table), and premp (# between 0 & 1 that sets the preemption value)");
    end            
    
            % see if fns is a file with list of files to process, instead of list of netcdf files
    if (length(fns)==1 && isfile(fns) && ~isnetcdf(fns))
        tbl = readtable(fns,"FileType", "text", "ReadVariableNames",true,"Delimiter",",");
        try
            fns = string(tbl.name);
        catch
            tbl = readtable(fns,"FileType", "text", "ReadVariableNames",false,"Delimiter",",");
            fns = string(tbl{:,1});
        end
    end 
    
    myfns = strtrim(string(fns));      % in case fns is a char array
    fns = strings(0);
    for i=1:length(myfns)
        finf = dir(fullfile(basedir, myfns(i)));
        for j=1:length(finf)
            fns(end+1) = fullfile(finf(j).folder, finf(j).name); %#ok<AGROW>
        end
    end
    
    
end
    
