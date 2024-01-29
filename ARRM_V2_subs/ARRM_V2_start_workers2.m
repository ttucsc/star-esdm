function [ppool, origCluster] = ARRM_V2_start_workers2(workers_requested, nnodes, cluster, jobname, mb_per_core)
% start maximum number of workers possible for parfor.
% Matlab's "local" profile seems to default to 12 cores, regardless of the
% number available, so query parcluster( ) for max possible and use that.
%
%   This version supports matlab MPS (Matlab Parallel Server)'s SLURM interface for starting workers on "redraider R2020b" cluster.
%
        % start threads if requested.  Otherwise, start worker processes on local or redraider cluster.
    if((isstring(workers_requested) && strcmp(workers_requested,"threads")) || strcmp(cluster,"threads"))
        ppool = ARRM_V2_start_threads();
        return;
    end

    if (~exist('nnodes','var')   || isempty(nnodes)),                           nnodes = 1;                                   end
    if (~exist('cluster','var') || isempty(cluster) || strlength(cluster)==0), cluster='local';                               end
    if (~exist('jobname','var') || isempty(jobname) || strlength(jobname)==0), jobname=sprintf("%s_%d", cluster,randi(1024)); end
    if (~exist('mb_per_core','var')),                                          mb_per_core = [];                              end
    
    myCluster = parcluster(cluster);
    origCluster = myCluster;
        
    fprintf("initial cluster returned:\n");
    disp(myCluster);
    orig_NumWorkers = myCluster.NumWorkers;          % on HPCC redraider , this is something like 10,000 or 100,000.  On local cluster, this *should* be # of actual cores available.
                                                     % On Neys (imac pro), this IS # of actual cores available.  On HPCC, this is always 1 the first time called, regardless of number of
                                                     % cores allocated by SLURM, which is available via getenv("SLURM_NTASKS").
    orig_workers_requested = workers_requested;
    try
        slurm_ntasks = str2double(getenv("SLURM_NTASKS"));
        if (isnan(slurm_ntasks) || slurm_ntasks == 0), slurm_ntasks = workers_requested; end
    catch
        fprintf(2, "NOTE:  could not get SLURM_NTASKS value via getenv(...) assuming %d (workers_requested)\n", workers_requested);
        slurm_ntasks = workers_requested;
    end
%   hostname = get_hostname(true);
    [~,hostname]=system('hostname');
    hostname = string(strtrim(hostname));
    
    fprintf("host:  %s\n", hostname);
    fprintf("SLURM_NTASKS %d\n", slurm_ntasks);
                % Local cluster
    if (strcmp(cluster,'local'))
        try
            fprintf("using local cluster\n");
            if (myCluster.NumWorkers < slurm_ntasks)   % earlier version, this was, incorrectly, <= instead of < .
                fprintf(2, "NOTE:  myCluster.NumWorkers was %d. SLURM_NTASKS: %d  hostname: %s\n", myCluster.NumWorkers, slurm_ntasks, hostname);
                myCluster.NumWorkers = slurm_ntasks;
                fprintf(   "setting myCluster.NumWorkers to %d\n", myCluster.NumWorkers);
            end
            if (workers_requested > slurm_ntasks)
                workers_requested = slurm_ntasks;
                fprintf("workers_requested adjusted to %d from %d\n", workers_requested, orig_workers_requested);
            end
%           cluster_nmax = slurm_ntasks;
%           myCluster.NumWorkers = max(ncores, cluster_nmax);
%           myCluster.saveProfile();
%           yourCluster = parcluster(cluster);
%           fprintf("\n\nafter saveProfile(), local cluster contains:\n");
%           disp(yourCluster);
        catch
            error("error:  cannot get number of cores to start workers, hostname:  %s", hostname);
        end
    else
                % redraider cluster

                    % commented out.  this seems to have changed, now
                    % returns 48, for some reason.  icsf 6/22
%         slurm_maxworkers = orig_NumWorkers;         % on HPCC redraider cluster, this is something like 10,000 or 100,000.
%         if (workers_requested > slurm_maxworkers)
%             workers_requested = slurm_maxworkers;
%             fprintf("NumWorkers adjusted to %d from %d\n", workers_requested, orig_workers_requested);
%         end
            % change the job storage location so we don't overrun the
            % user's disk space allocation.
        try
            username = string(getenv("USER"));
            storageLocation = sprintf("/lustre/scratch/%s/%s", username,fix_filename(cluster));
            if (~isfolder(storageLocation)), mkdir(storageLocation); end
            if (exist('jobname','var') && ~isempty(jobname) && strlength(jobname)>0)
                storageLocation = sprintf("%s/%s", storageLocation, fix_filename(jobname));
                    if (~isfolder(storageLocation)), mkdir(storageLocation); end
            end
            fprintf("job storage location:  %s\n", storageLocation);
        catch
            error("error:  problems setting up storageLocation for SLURM");
        end
        if (~isfolder(storageLocation)), error("error:  cannot create job storage location %s for SLURM's job output", storageLocation); end
        myCluster.JobStorageLocation = storageLocation;
        if (strncmp(cluster,'redraider',9))
            myCluster.AdditionalProperties.Partition = 'nocona';
        elseif (strncmp(cluster,'quanah',6))
            myCluster.AdditionalProperties.Partition = cluster;
        end
        myCluster.AdditionalProperties.NumberOfNodes=nnodes;
        if (~isempty(mb_per_core))
            myCluster.AdditionalProperties.MemUsage = sprintf('%d',floor(mb_per_core));         % memory spec is in MB/core, but mus be passed as a string.
        end
        fprintf("myCluster.NumWorkers:  %d\n", myCluster.NumWorkers);
%         if (cluster_nmax > 512)
%             fprintf("cluster_nmax too large.  Setting NumWorkers to %d and saving\n", ncores);
%             myCluster.NumWorkers = ncores;
%             myCluster.saveProfile();
%             yourCluster = parcluster();
%             fprintf("\n\nafter saveProfile(), %s cluster contains:\n", yourCluster.Profile);
%             disp(yourCluster);
%         end
    end
    
    % debugging:  check how many files are open:
    fids = fopen("all");
    fprintf("total number of files currently open:  %d\n\n", length(fids));
    
    fprintf("\nparCluster:  %d NumWorkers available\n", myCluster.NumWorkers);
    disp(myCluster);
    
        % finally, start the parallel pool.
        % check to see if it's already running.
    ppool = gcp('nocreate');

    if (isempty(ppool))
        fprintf("starting ppool with %d workers;  original cluster workers available=%d\n", workers_requested, orig_NumWorkers);
        
        mystart = tic();
        try
            ppool = parpool(myCluster, workers_requested);
            elapsed = toc(mystart);
            fprintf("time to create ppool with %d workers (%d requested): %.1f secs\n", ppool.NumWorkers, orig_workers_requested, elapsed)
        catch me
            elapsed = toc(mystart);
            fprintf(2,"error:  failed to start ppool with %d workers after attempting for %.1f secs\n", workers_requested, elapsed);
            rethrow(me);
        end
    else
        if (~isa(ppool, "parallel.ProcessPool"))
            error("error:  parpool already running, but is a '%s', not a 'parallel.ProcessPool'", class(ppool));
        end
        fprintf(2, 'NOTE:  Parallel pool already running %d workers.', ppool.NumWorkers); 
    end
    fprintf(2, 'NOTE:  Number of processes: %d\n', ppool.NumWorkers);
end 

function s = fix_filename(s)
% replaces spaces in s with underscores
% useful for makeing filenames from strings that might contain spaces or other chars not allowed in filenames.

    s = strsplit(string(s),filesep);

    if (length(s) == 1)
        cs = char(s);
        ix = ~isstrprop(cs,'alphanum');
        cs(ix) = '_';
        s = string(cs);
    else
        for i=1:length(s)
            s(i) = fix_filename(s(i));
        end
        s = join(s,filesep);
    end    
    s=string(s);
end
