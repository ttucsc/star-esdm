function rerejoin(basedir, newer_str, vnames, do_parfor, dry_run)

    if (~exist("do_parfor","var") || isempty(do_parfor)), do_parfor = false; end
    if (~exist("dry_run","var")   || isempty(dry_run)),   dry_run = false;   end
    
    if (dry_run), do_parfor = false; end
    
    logname = sprintf("%s/rerejoin.log", basedir);
    diary(logname);
    
    if (exist("newer_str","var") && ~isempty(newer_str))
        newer_num = datenum(newer_str);
        newer_str = datestr(newer_num,"yyyy/mm/dd HH:MM:");
        fprintf("\n----------------rerejoin:  %s, %s, skipping where .nc file is newer than %s\n", datestr(now,"yyyy/mm/dd HH:MM"),basedir, newer_str); 
    else
        fprintf("\n----------------rerejoin:  %s, %s, s\n", datestr(now,"yyyy/mm/dd HH:MM"),basedir); 
        newer_str = "";
    end
    
    dirinfo = dir(basedir);
    dirinfo = dirinfo([dirinfo.isdir]);
    dirnames = string({dirinfo.name});

            % find all the folders of possible interest.
    if (exist("vnames","var") && ~isempty(vnames))
        for i=1:length(vnames)
            keeps =  contains(dirnames, vnames(i));
            if (i==1)
                keepers = keeps;
            else
                keepers = keepers | keeps;
            end
        end
    else
        vnames = "";       
    end
    % for i=1:length(dirinfo); fprintf("%2d %d %s\n", i, keepers(i), dirinfo(i).name); end    
    dirinfo = dirinfo(keepers);
    ndirs = length(dirinfo);
    
    keepers = false(size(dirinfo));
    
            % keep only the ones that look like downscaled folders.
    for i=1:ndirs
        tt = split(dirinfo(i).name,"_");
        if (length(tt) ~= 5), continue; end
        keepers(i)  = true;
    end
%   dirnames = dirnames(keepers);
    dirinfo = dirinfo(keepers);
    
    ndirs = length(dirinfo);
    
        
    if (ndirs==0)
        diary off;
        error("error:  no matching folders for %s [%s] %s\n", basedir, join(vnames,", "), newer_str);
    else
        fprintf("found:  %d folders to process (including completed ones):\n", ndirs);
%         for i=1:ndirs
%             fprintf("%2d %s %s\n",i,  datestr(dirinfo(i).datenum, "yyyy/mm/dd HH:MM"), dirinfo(i).name);
%         end        
    end
    
            % find the joined netcdf files newer than the cutoff date.
            % We don't want to re-join these ones.
    finfo = dir(fullfile(basedir, "*.nc"));
    tstamps = [finfo.datenum];
    finfo = finfo(tstamps >= newer_num);
    finished_names = string({finfo.name});
    
    if (~isempty(finished_names))
        fprintf("Completed:\n");
        for i=1:length(finished_names)
            fprintf("%s %s\n", finfo(i).date, finfo(i).name);
        end
        fprintf("\n");
    end
    
%   ppool = ARRM_V2_start_workers2(5, 1, "local");
    if (do_parfor)
        parfor i=1:ndirs        
            run_rejoin(i, ndirs, dirinfo(i), finished_names, dry_run);
        end
    else
        for i=1:ndirs        
            run_rejoin(i, ndirs, dirinfo(i), finished_names, dry_run);
        end        
    end
    
    diary off
end

function run_rejoin(i, ndirs, dirinfo, finished_names, dry_run)

    tt = string(split(dirinfo.name,"_"));
    model    = tt(1);
    ensemble = tt(2);
    varname  = tt(3);
    scenario = tt(4);
    grgn     = tt(5);
    
    outname = sprintf("downscaled.%s.%s.%s.%s.%s.nclimgrid.nclimgrid.1950.2100.tllc_time.nc", model, ensemble, varname, scenario, grgn);
    if (any(strcmp(finished_names, outname)))
        outinfo = dir(fullfile(dirinfo.folder,outname));
        fprintf("%4d of %4d: already rerun.  skipping: %s %9d      %s\n", i, ndirs, datestr(outinfo.datenum, "yyyy-mm-dd HH:MM"), outinfo.bytes, outname);
        return;
    end
    fprintf("getting directory info for %s...", dirinfo.name);
    subfinfo = dir(fullfile(dirinfo.folder, dirinfo.name,"downscaled.*.nc"));
    nncs = length(subfinfo);
    fprintf("done\n");
    try
        fn=subfinfo(1).name;
        ss=string(split(fn,["_","."]));
        gsize=str2double(split(ss(end-1),"x"));
        if (length(gsize)==2)
            gridsize = gsize;
        end
    catch
        if (length(subfinfo) < 500)
            gridsize = [2,2];
        else
            gridsize = [1,1];
        end
    end
    if (all(to_row(gridsize) == [1,1]))
        totfiles = 1508;
    else
        totfiles = 390;
    end
    
    if ( nncs < totfiles)
        fprintf("%4d of %4d: incomplete (%d of %d files files present).  Skipping %s\n", i, ndirs, nncs, totfiles, dirinfo.name);
        return;
    end
    
    

    if (dry_run)
        fprintf("rejoin_downscaled_minifiles(""nclimgrid\"",  ""%s"". ""%s"". ""%s"", ""%s"", ""%s"", ""cmip6"", true, ""nclimgrid"", [%d,%d]. [1950,2100], true)\n", model, ensemble, varname, scenario, grgn, gridsize); 
    else
        try                             % region, model, ensemble, varnames,scenario, grgn, model_set, do_overwrite, obs_src, gridsize, mdl_yrs, include_zvals, pdf_map_method)
            fprintf("\n----------rerejoining: %4d of %4d: %-15s %-10s %-10s  %-10s %4s  nfiles: %d  grid [%d,%d], %s\n\n", i, ndirs, model, ensemble, varname, scenario, grgn, length(subfinfo), gridsize, dirinfo.name ); 
            fprintf("first file: %s\n", subfinfo(1).name);
            rejoin_downscaled_minifiles("nclimgrid",  model, ensemble, varname, scenario, grgn, "cmip6", true, "nclimgrid", gridsize, [1950,2100], true); 
        catch me
            fprintf("********** error:  problem re-joining %s\n", dirnames(i));
            report_me_error(me)
        end
    end
end
