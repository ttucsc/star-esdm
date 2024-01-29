function make_matfiles(ncnames, skip_if_present, do_parfor)

    if (~exist("skip_if_present","var") || isempty(skip_if_present)), skip_if_present = false; end
    if (~exist("do_parfor","var")       || isempty(do_parfor)),       do_parfor       = false; end

    ncnames = string(ncnames);
    
    if (do_parfor && length(ncnames) > 1)
        parfor i=1:length(ncnames)
            make_matfiles(ncnames(i), skip_if_present, false);
        end
    else
        for i=1:length(ncnames)
            fname = canonical_path(ncnames(i));  % canonical_name, in case ncname or basedir is actually a symbolic link.
            if (isfolder(fname))
                basedir = fname;
                ncinfos=dir(fullfile(basedir,"*.nc"));
                mynames = string({ncinfos.name})';
                if (do_parfor)
                    parfor j=1:length(mynames)
                        make_matfiles(fullfile(basedir,mynames(j)), skip_if_present, false);
                    end
                else
                    for j=1:length(mynames)
                        make_matfiles(fullfile(basedir,mynames(j)), skip_if_present, false);
                    end
                end
           else
                if (isQCnetcdf(fname))
                    try
                        if (skip_if_present)
                            [matok, siteTbl, matname] = check_for_QC_matfile(fname);
                            if (matok)
                                fprintf("matfile already exists and timestamp matches:   %s\n", matname);
                            end
                        else
                            matok = false;
                        end
                        if (~matok)
                            fprintf("creating matfile for %s\n",fname);
                            siteTbl = QC_get_site_table(fname,"loadTblvars",true, "tryMatfile", false);
                            [cdir, cbase, ~] = fileparts(fname);
                            matname = fullfile(cdir, sprintf("%s%s", cbase, ".mat"));
                            save(matname,"siteTbl");                   
                            [d,fn,~] = fileparts(canonical_path(fname));
                            matname = fullfile(d,sprintf("%s%s", fn, ".mat"));
                        end
                        system(sprintf("ls -l %s %s", fname, matname));
                        fprintf("%s:  first 3 rows:\n", fname);
                        QC_tbl(siteTbl(1:3,:));
                        fprintf("\n");
                    catch
                        fprintf(2, "\nProblem creating matfile for %s\n\n", fname);
                    end
                else
                    fprintf(2, "error:  %s is not a QC netcdf file\n\n", fname);
                end
            end
        end
    end
end

function [matok, siteTbl, matname, ncName] = check_for_QC_matfile(ncName)
% returns flag (& siteTbl) if ncName has a QC_site_table matfile and timestamps match.  
%   returns false if no.

    matname = [];
    siteTbl = [];
    matok = false;
    if (~isQCnetcdf(ncName)), return; end

    [p,d,~] = fileparts(ncName);
    siteTbl = [];
    
    matname = fullfile(p,sprintf('%s%s',d,".mat"));
    if (~isfile(matname))       % if matfile doesn't exist, check canonical path as well.
        ncNameC = canonical_path(ncName);
        if (~strcmp(ncNameC, ncName))
            [p2,d2,~] = fileparts(ncNameC);
            matname  = fullfile(p2, sprintf('%s%s',d2,'.mat'));
        end
    end

    if (isfile(matname))
        try
                % look for siteTbl in matfile.
            minfo=whos(matfile(matname),'siteTbl');
            if (~isempty(minfo) && strcmp(minfo.class,'table'))
                load(matname,"siteTbl");
                finfo = dir_canonical(ncName);
                if (isfield(siteTbl.Properties.UserData,'table_source') && isequal(finfo, siteTbl.Properties.UserData.table_source))
                    matok = true;
                else
                    matok = false;
                end
                return;
                if (tryMatfile==1 && strcmpi(ext,".nc"))
                    fprintf("NOTE:  timestamp/size mismatch between matfile and netcdf file\n\t%s\n\t%s\nUsing nc file\n", matname, ncName);
                    siteTbl = [];
                end
            end
        catch
            siteTbl = [];
            matok = false;
            error("error reading %s", matname);
        end
    end     
end
