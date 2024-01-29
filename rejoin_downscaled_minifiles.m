function rejoin_downscaled_minifiles(region, model, ensemble, varnames,scenario, grgn, model_set, do_overwrite, obs_src, gridsize, mdl_yrs, include_zvals, pdf_map_method)
% function rejoin_downscaled_minifiles(region, model, varnames, ensemble,scenario, do_overwrite, is_livneh_run, gridsize, quad_or_single, mdl_yrs, pdf_map_method)
%
%   usage:  rejoin_downscaled_minifiles(region, model, varnames, ensemble,scenario, do_overwrite + optional:  is_livneh_run, gridsize, quad_or_single, mdl_yrs, pdf_map_method)
%
%   currently set up to join downscaled livneh data.  Will need some modifications to work on other gridded obs data.
%
%   varnames:  can be a single variable or multiple variables.  If single, assumes the primary variable is the first
%   (for folder and file names), but will also copy over any additional variables.  For example, ["tasmax", "packed_zvals"]  (or "zvals", 12/30/31, changed name to zvals from packed_zvals) 
%        if varnames are a single row, they are all read from the same set of files and written to the same output fule.  
%        if they are on separate rows, then they are read from separate intermediate files and written to separate output files.
%            This is to handle the change where extended-run                zvals are now written to separate files. 
%
%   12/2/21     changing output to create TLLC output instead of LLT output.
%   12/2/21     changed output for gridded data to TLLC_time, chunked [1,1,nyrs*365]
%   2/4/23      modified to handle separate z_vals files.  (zvals are now written to a 2nd intermdiate output file, and need to be combined into a separate large file.) 
%
%   

    if (nargin < 11), help(mfilename); end
    if (islogical(obs_src))
        if (obs_src)
            obs_src = "livneh";
        else
            obs_src = "nrcan";
        end
    end
    if (~exist('pdf_map_method','var')), pdf_map_method = []; end
    if (~exist("include_zvals","var")), include_zvals = false; end
    if (~isstring(varnames)), error("error:  varnames must be string, not char"); end
    n_var_rows = size(varnames,1);
    varname=varnames(1,1);

    add_src_to_dirs = true;     % kludge.  s/b an input parameter.  for changes made 4/23-5/23 to separate runs using different gridded obs sets.
                                %  earlier, we always worked with a single gridded obs set.
    
    if (isfolder("/Volumes/lacie_1/data"))
        basedir = "/Volumes/lacie_1/data/downscaled/arrm_v2";
    elseif (isfolder("/lustre/"))
        basedir = "/lustre/scratch/iscottfl/downscaled";
    end
    
    if (isempty(pdf_map_method))
        if (strlength(grgn)==0)
%           subdir = sprintf("%s_%s_%s_%s_%s", model, ensemble, varname, scenario, quad_or_single);
            subdir = sprintf("%s_%s_%s_%s",    model, ensemble, varname, scenario);
        else
%           subdir = sprintf("%s_%s_%s_%s_%s_%s", model, ensemble, varname, scenario, grgn, quad_or_single);
            subdir = sprintf("%s_%s_%s_%s_%s",    model, ensemble, varname, scenario, grgn);
        end
    else
        if (strlength(grgn)==0)
%           subdir = sprintf("%s_%s_%s_%s_%s_%s", model, ensemble, varname, scenario, quad_or_single, pdf_map_method);
            subdir = sprintf("%s_%s_%s_%s_%s",    model, ensemble, varname, scenario,                 pdf_map_method);
        else
%           subdir = sprintf("%s_%s_%s_%s_%s_%s_%s", model, ensemble, varname, scenario, grgn, quad_or_single, pdf_map_method);
            subdir = sprintf("%s_%s_%s_%s_%s_%s",    model, ensemble, varname, scenario, grgn,                 pdf_map_method);
        end
    end
    
    basedir = fullfile(basedir, model_set, region);
    if (add_src_to_dirs)
        basedir = fullfile(basedir, obs_src);
        subdir = sprintf("%s_%s",subdir, obs_src);
    end
    if (strncmp(obs_src,"stations",8))
        subdir = sprintf("all_%s", subdir);
    end
    
    minifiledir = fullfile(basedir,scenario, subdir);
    destdir = basedir;

    if (~isfolder(destdir)),     error("error:  dest folder %s doesn't exist", destdir); end
    if (~isfolder(minifiledir)), error("error:  minifile folder %s doesn't exist", minifiledir); end

    if (~strcmp(obs_src,"stations"))
        chunksize = [inf,1,1];
        do_rotate = true;
    else
        chunksize = [];
        do_rotate = false;
    end

    for vn = 1:n_var_rows
        vnames = varnames(vn,:);
        vname = vnames(1);
    
        if (strcmp(obs_src, "nclimgrid"))
            if (strlength(grgn)==0)
                tmplnm      = sprintf("downscaled.%s.%s.%s.%s.nclimgrid.%s.%d.%d.grid_*_%02dx%03d.nc", model, ensemble, vname, scenario, region, mdl_yrs, gridsize);
                outnm       = sprintf("downscaled.%s.%s.%s.%s.nclimgrid.%s.%d.%d.tllc_time.nc",        model, ensemble, vname, scenario, region, mdl_yrs);
            else
                tmplnm      = sprintf("downscaled.%s.%s.%s.%s.%s.nclimgrid.%s.%d.%d.grid_*_%02dx%03d.nc", model, ensemble, vname, scenario, grgn, region, mdl_yrs, gridsize);
                outnm       = sprintf("downscaled.%s.%s.%s.%s.%s.nclimgrid.%s.%d.%d.tllc_time.nc",        model, ensemble, vname, scenario, grgn, region, mdl_yrs);
            end
        elseif (strcmp(obs_src, "livneh"))
            if (strlength(grgn)==0)
                tmplnm      = sprintf("downscaled.%s.%s.%s.%s.livneh_1_16th.%s.%d.%d.grid_*_%02dx%03d.nc", model, ensemble, vname, scenario, region, mdl_yrs, gridsize);
                outnm       = sprintf("downscaled.%s.%s.%s.%s.livneh_1_16th.%s.%d.%d.tllc_time.nc",        model, ensemble, vname, scenario, region, mdl_yrs);
            else
                tmplnm      = sprintf("downscaled.%s.%s.%s.%s.%s.livneh_1_16th.%s.%d.%d.grid_*_%02dx%03d.nc", model, ensemble, vname, scenario, grgn, region, mdl_yrs, gridsize);
                outnm       = sprintf("downscaled.%s.%s.%s.%s.%s.livneh_1_16th.%s.%d.%d.tllc_time.nc",        model, ensemble, vname, scenario, grgn, region, mdl_yrs);
            end
        elseif (strcmp(obs_src,"nrcan"))
            if (isempty(pdf_map_method))
                if (strlength(grgn)==0)
                    tmplnm      = sprintf("downscaled.%s.%s.%s.%s.%s.%d.%d.grid_*_%02dx%03d.nc", model, ensemble, vname, scenario, region, mdl_yrs, gridsize);
                    outnm       = sprintf("downscaled.%s.%s.%s.%s.%s.%d.%d.tllc_time.nc",        model, ensemble, vname, scenario, region, mdl_yrs);
                else
                    tmplnm      = sprintf("downscaled.%s.%s.%s.%s.%s.%s.%d.%d.grid_*_%02dx%03d.nc", model, ensemble, vname, scenario, grgn, region, mdl_yrs, gridsize);
                    outnm       = sprintf("downscaled.%s.%s.%s.%s.%s.%s.%d.%d.tllc_time.nc",        model, ensemble, vname, scenario, grgn, region, mdl_yrs);
                end
            else
                if(strlength(grgn)==0)
                    tmplnm      = sprintf("downscaled.%s.%s.%s.%s.%s.%s.%d.%d.grid_*_%02dx%03d.nc", model, ensemble, vname, scenario, region, pdf_map_method, mdl_yrs, gridsize);
                    outnm       = sprintf("downscaled.%s.%s.%s.%s.%s.%d.%s.%d.tllc_time.nc",        model, ensemble, vname, scenario, region, pdf_map_method, mdl_yrs);
                else
                    tmplnm      = sprintf("downscaled.%s.%s.%s.%s.%s.%s.%s.%d.%d.grid_*_%02dx%03d.nc", model, ensemble, vname, scenario, grgn, region, pdf_map_method, mdl_yrs, gridsize);
                    outnm       = sprintf("downscaled.%s.%s.%s.%s.%s.%d.%s.%s.%d.tllc_time.nc",        model, ensemble, vname, scenario, grgn, region, pdf_map_method, mdl_yrs);
                end
            end
        elseif (strncmp(obs_src, "stations",8))
            if(strlength(grgn)==0)
                tmplnm      = sprintf("downscaled_all.%s.%s.%s.%s.%d.%d.grid_*.nc", model, ensemble, vname, scenario, mdl_yrs);
                outnm       = sprintf("downscaled_all.%s.%s.%s.%s.%d.%d.nc",         model, ensemble, vname, scenario, mdl_yrs);
            else
                tmplnm      = sprintf("downscaled_all.%s.%s.%s.%s.%s.%d.%d.grid_*.nc", model, ensemble, vname, scenario, grgn, mdl_yrs); % downscaled_all.ACCESS-CM2.r1i1p1f1.tasmax.ssp585.gn.1950.2100.grid_N62W162.nc        
                outnm       = sprintf("downscaled_all.%s.%s.%s.%s.%s.%d.%d.nc",   model, ensemble, vname, scenario, grgn, mdl_yrs);
            end
        end
        
        template    = fullfile(minifiledir, tmplnm);    % "/Volumes/lacie_1/data/downscaled/arrm_v2/GFDL_perfect_model/conus/GFDL-HIRAM-C360_cm3/r1i1p1/tasmax/rcp85-X1X3X4/downscaled.GFDL-HIRAM-C360_cm3.r1i1p1.tasmax.rcp85-X1X3X4.gfdl_PM.conus.2086.2115.grid_*.nc"
        outname     = fullfile(destdir,     outnm);     % "/Volumes/lacie_1/data/downscaled/arrm_v2/GFDL_perfect_model/conus/downscaled.GFDL-HIRAM-C360_cm3.r1i1p1.tasmax.rcp85-X1X3X4.gfdl_PM.conus.2086.2115.nc";
    
        fprintf("template:    %s\n", template);
        fprintf("varname:     %s\n", vname);
        fprintf("outname:     %s\n", outname);
        if (length(vnames)>1)
            fprintf("addtl vars:  %s\n", vec2string(vnames(2:end)));
        end
        
        if (~isfolder(minifiledir)), error("error:  minifile folder %s does not exist", minifiledir); end
        if (~isfolder(destdir)), error("error:  destination folder %s does not exist", destdir); end
        finfo = dir(template);
        nfiles = length(finfo);
    %   [status, nfiless] = system(sprintf("ls %s 2>/dev/null | wc -l", template));
    %   if (status ~= 0), error("error checking for input files: %d", status); end
    %   nfiles = str2double(nfiless);
        if (~(nfiles > 0)), error("error:  no files matching template %s", template); end
        
    %   join_minifiles(template,varnames,outname,"verbose",true, "overwrite", do_overwrite,'format','netcdf4','rotate',false,"clean_metadata",true);
        join_minifiles(template,vnames,outname,"verbose",true, "overwrite", do_overwrite,'format','netcdf4','rotate',do_rotate,"clean_metadata",true,"chunksize",chunksize, "unlimited",true, "include_zvals", include_zvals);
    end
end
