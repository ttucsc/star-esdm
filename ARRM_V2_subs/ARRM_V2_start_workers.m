function ppool = ARRM_V2_start_workers(nmax, partn, jobname)
% start maximum number of workers possible for parfor.
% Matlab's "local" profile seems to default to 12 cores, regardless of the
% number available, so query parcluster( ) for max possible and use that.
%

    if (~exist('partn','var') || isempty(partn) || strlength(partn)==0)
        [~,on] = get_hostname(true);
        if (on.hpcc_system)
            partn = "redraider R2020b";
        else
            partn = "local";
        end
    end
    myCluster = parcluster(partn);
    
    cluster_nmax = myCluster.NumWorkers;
    fprintf("initial cluster returned:\n");
    disp(myCluster);

    if (~strcmp(partn,"local"))
        slurm_ntasks = str2double(getenv("SLURM_NTASKS"));
        fprintf("SLURM_NTASKS %d\n", slurm_ntasks);
        if (~exist("nmax","var") || isempty(nmax) || nmax == 0)
            if (cluster_nmax == 1)
                try
                    if (cluster_nmax == 1)
                        fprintf("NumWorkers was 1.  Using SLURM_NTASKS\n");
                    else
                        fprintf("using local cluster. Using SLURM_NTASKS\n");
                        if (nmax > slurm_ntasks - 4)
                            old_nmax = nmax;
                            nmax = slurm_ntasks-4;
                            fprintf("nmax adjusted to %d from %d\n", nmax, old_nmax);
                        end
                    end
                    cluster_nmax = slurm_ntasks;
                    myCluster.NumWorkers = cluster_nmax;
                    myCluster.saveProfile();
                    yourCluster = parcluster(partn);
                    fprintf("\n\nafter saveProfile(), local cluster contains:\n");
                    disp(yourCluster);
                catch
                    error("error:  cannot get number of cores to start workers");
                end
            else
                nmax = cluster_nmax;
            end
        else
                % change the job storage location so we don't overrun the
                % user's disk space allocation.
                % User must have a folder in their scratch area named nocona,
                % quanah or xlquanah 
                %
            try
                username = string(getenv("USER"));
                storageLocation = sprintf("/lustre/scratch/%s/%s", username,fix_filename(partn));
                if (~isfolder(storageLocation)), mkdir(storageLocation); end
                if (exist('jobname','var') && ~isempty(jobname) && strlength(jobname)>0)
                    storageLocation = sprintf("%s/%s", storageLocation, fix_filename(jobname));
                        if (~isfolder(storageLocation)), mkdir(storageLocation); end
                end
            catch
                error("error:  problems setting up storageLocation for SLURM");
            end
            if (~isfolder(storageLocation)), error("error:  cannot create job storage location %s for SLURM's job output", storageLocation); end
            myCluster.JobStorageLocation = storageLocation;    

            fprintf("valid NumWorkers:  %d\n", cluster_nmax);
            if (cluster_nmax > 512)
                fprintf("cluster_nmax too large.  Setting NumWorkers to %d and saving\n", nmax);
                myCluster.NumWorkers = nmax;
                myCluster.saveProfile();
                yourCluster = parcluster();
                fprintf("\n\nafter saveProfile(), %s cluster contains:\n", yourCluster.Profile);
                disp(yourCluster);
            end
        end
    end
    
    fprintf("\nparCluster:  %d NumWorkers available\n", cluster_nmax);
    disp(myCluster);
    
    if ((~exist('nmax','var') || isempty(nmax)) && cluster_nmax > 35)
        nmax = cluster_nmax - 4; 
    end 
        
    ppool = gcp('nocreate');

    if (isempty(ppool))
        fprintf("starting ppool with %d workers;  cluster workers available=%d\n", nmax, myCluster.NumWorkers);
        
        tic();
        ppool = parpool(myCluster, nmax);
        elapsed = toc();
        fprintf("time to create ppool: %.1f\n", elapsed)
%       ppool = parpool();
    else
        fprintf('Parallel pool already running.  '); 
    end
    fprintf('Number of processes: %d\n', ppool.NumWorkers);
end

function s = fix_filename(s)
% replaces spaces in s with underscores
% useful for makeing filenames from strings that might contain spaces.

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
end