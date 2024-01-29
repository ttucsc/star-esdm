function outname = compress_netcdf( inname, outname, varargin )
%
%   kwd/value pairs:
%       outdir          string          name of directory to write output.  If empty, writes to same folder as inname.
%       shuffle         true/false      turns on or off byte shuffle.  defaults: true for integer, false for float types
%       chunksize       [x,y,z]         set chunking.  default:  [ndays,1,1] for LLT files, [1,1,ndays] for std files
%                                           if more than 1 climate variable, use 1 row for each climate variable's
%                                           chunksize
%       chunktype       string          "time" or "latlon".  Sets whether to chunk along time or lat/lon dimension.
%                                           used only if chunksize is empty.
%       deflatelevel    [0..9]          deflate level [default:  [] ].  Set to [] to turn off compression.  
%                                                                           using 0 does no compression, but adds some
%                                                                           overhead.  See NetCDF documentation.
%       int_scaling     true/false or [add_offset,scale_factor]  write output using 16-bit scaled integers "packed format. 
%                                           where unpacked_value = packed_value * scale_factor + add_offset
%                                                     add_offset = dataMin
%                                                   scale_factor = (dataMax - dataMin) / (2^n - 1)   n==# of bits: 16
%                                           default (if int_scaling == true):
%                                               1.  Assumes only 1 climate variable in file
%                                               2.  Reads data, finds range (min, max), and calculates add_offset
%                                                       and scale_factor
%                                               3.  Aborts if variable is larger than maxmem bytes
%       varnames        string array of variables to compress
%                           if empty, will compress any climate variable found in the file
%       overwrite       true/[false]  if false, will not overwrite existing output file
%
%       maxmem          maximum memory to use for reading a variable.  S/b less than 1/3 total system memory.
%                           default:    on Neys:  40 GB 
%                                       on HPCC:   2 GB     (about 6 GB mem available on standard HPCC qlogin.)
%                                                           BUT you can request up to about 256 GB on a qlogin.
%                                                           See HPCC documentation for how to request specific memory
%                                                           amount.
%       verbose         true/false or 0,1 or 2.  0 (false) don't report info on variables  0default]
%                                                1 (true)  report info on variables
%                                                2         report info and attributes on variables
%
%   NOTE:  This won't work properly on station netcdf files!  Still need to add that capability.
%
%   To Do:
%       copy_by_parts needs to be smarter and work by chunks, instead of by lat/lon or time.


%   addpath("Users/iscottfl/Documents/MATLAB",	"/Users/iscottfl/Desktop/atmos2/ARRM_V2".	"/Users/iscottfl/Desktop/atmos2/ARRM_V2/util_general",	"/Users/iscottfl/Desktop/atmos2/ARRM_V2/util_nc",	"/Users/iscottfl/Desktop/atmos2/ARRM_V2/ncdf",	"/Users/iscottfl/Desktop/atmos2/ARRM_V2/ARRM_V2_subs",	"/Users/iscottfl/Desktop/atmos2/ARRM_V2/utils"
    addpath("./util_general",	"./util_nc",	"./ncdf",	"./ARRM_V2_subs",	"./utils");

    t1 = tic();
    [ncin, outname, do_overwrite, varnames, shuffle, chunksize, deflatelevel, calc_scaling, add_offset, scale_factor, verbose, maxmem_orig] = init_params(inname, outname, varargin{:});
    
    fprintf("\n-------------compress_netcdf.m %s---------------\n", datestr(now));
    fprintf("reading file:  %s\n", basename(inname));
    fprintf("from folder:   %s\n", dirname(inname));
    fprintf("writing file:  %s\n", basename(outname));
    fprintf("to folder:     %s\n", dirname(outname));
    
    fprintf("\nCompression settings:\n");
    if (isempty(add_offset) && ~calc_scaling)
        fprintf("output datatype:  float\n");
    elseif (~calc_scaling)
        fprintf("output datatype:  uint16\n");
        fprintf("add_offset:       %.4f\n", add_offset);
        fprintf("scale_factor:     %.8f  (1/ = %.1f)\n", scale_factor, 1/scale_factor);
    else
        fprintf("output datatype:  uint16\n");
        fprintf("add_offset:       (calculate)\n");
        fprintf("scale_factor:     (calculate)\n");
    end
    if (isempty(deflatelevel))
        fprintf("deflate level:    none\n");
    else
        fprintf("deflate level:    %d\n", deflatelevel);
    end
    
    if (isempty(chunksize))
        fprintf("chunksize:        none\n");
    else
        fprintf("chunksize:        [%d, %d, %d]\n", chunksize);
    end
    fprintf("shuffle:          %d\n", shuffle);
    
    fprintf("Variables:\n");
    for i=1:length(varnames)
        fprintf("\t%s\n", varnames(i));
    end
    
    if (calc_scaling)
        [add_offset, scale_factor, vdata] = calc_scale_factors(ncin, varnames);
    else
        vdata = [];
    end

    ncout = update_compression_settings(ncin, outname, varnames, shuffle, chunksize, deflatelevel, add_offset, scale_factor);

    if (~strcmp(ncout.Format, ncin.Format))
        fprintf(2, "NOTE:  format changed from %s to %s\n", ncin.Format, ncout.Format);
    end
    
    netcdf_init(ncout, do_overwrite); 

    nvars = length(varnames);

    for i=1:nvars
        vname = varnames(i);
        vin  = ncin.getvar(vname);
        
        if (vin.total_nbytes < ncObj.maxmem)
            if (is_climate_variable(vname))
                copy_var_whole(   ncin, ncout, vin, vname, i, verbose, vdata);
            else
                copy_var_whole(   ncin, ncout, vin, vname, i, verbose);
            end
        else
            copy_var_by_parts(ncin, ncout, vin, vname, i, verbose);
        end
    end

    elapsed = toc(t1); 
    fprintf("done rewriting.  Output is in %s\n", outname);
    fprintf("elapsed:  %s\n\n", timestr(elapsed));
    
    ncObj.maxmem(maxmem_orig);
    
end

function [add_offset, scale_factor, vdata] = calc_scale_factors(ncin, varnames)

    varname = find_climate_varname(varnames,"all");
    
    if (length(varname) > 1), error("error:  too many climate variables to calculate scale factors: %s", join(varname," ")); end
    
    vv = ncin.getvar(varname);
    if (vv.total_nbytes > ncObj.maxmem)
        fprintf(2,"Variable too large to fit in memory:");
        disp(vv);
        error("error:  can't calculate scale factors for variable that won't fit in memory");
    end

    vdata = nc.readvar(varname);
    
    vmin = min(vdata(:));
    vmax = max(vdata(:));
    add_offset = vmin;
    scale_factor = (vmax - vmin)/65534;
    
end

function copy_var_whole(nc, ncout, vin, varname, i, verbose, vdata)

    vout = ncout.getvar(varname);

    fprintf("reading variable %2d %-20s:  %12d bytes, size: %s\n",  i, vin.Name, vin.total_nbytes(), vec2string(vin.Size,"brackets",'[]'));
    
    if (~exist("vdata","var") || isempty(vdata))
        vout.vdata = nc.readvar(varname);
    end
%     nc.loadvar(vname);

    if (verbose > 0)
        disp(vout);
        if (verbose > 1)
            for j=1:length(vout.Attributes)
                disp(vout.Attributes(j))
            end
        end
    end
            

    fprintf("writing variable %2d %-20s:  %12d bytes, size: %s\n",  i, vin.Name, vin.total_nbytes(), vec2string(vin.Size,"brackets",'[]'));

    ncout.writevar(varname);
    pause(.025);
end


function copy_var_by_parts(ncin, ncout, vin, varname, i, verbose)

    fprintf("copying variable %2d %-20s by parts:  %12d bytes, size: %s\n",  i, vin.Name, vin.total_nbytes(), vec2string(vin.Size,"brackets",'[]'));
    
    vout = ncout.getvar(varname);
    if (verbose > 0)
        disp(vout);
        if (verbose > 1)
            for j=1:length(vout.Attributes)
                disp(vout.Attributes(j))
            end
        end
    end
            

    [~, ~, ~,~, ~, ~, timix] = ncdf_get_llt_dimnames(vin);
    is_llt = timix == 1;
    nx = vin.Size(1);
    ny = vin.Size(2);
    nz = vin.Size(3);
    if (is_llt)
        ntimes = nx;
        for iz=1:nz
            for iy=1:ny
                start=[1, iy,iz];
                vdata = ncin.readvar(varname, start, [ntimes,1,1]);
                ncout.writevar(varname, vdata, start);
            end
            show_progress(iz, nz);
        end
    else
        ntimes = nz;
        for iz=1:ntimes
            start = [1,1,iz];
            vdata = ncin.readvar(varname, start, [nx,ny,1]);
            ncout.writevar(varname, vdata, start);
            show_progress(iz,nz);
        end
    end
    
end

function netcdf_init(ncout, do_overwrite)

%     outvars = setdiff(ncout.varlist, varnames);
%   ncout.writeschema(do_overwrite, "Dimensions", true, "Variables", false, "nofillmode", true);
    ncout.writeschema(do_overwrite, "Dimensions", true, "Variables", false, "nofillmode", false);   % Ouch.  matlab changes _FillValue to "disable" if nofillmode is set to true.
    
end

function ncout = update_compression_settings(ncin, outname, varnames, shuffle, chunksize, deflatelevel, add_offset, scale_factor)

    ncout = ncin.clone;
    ncout.Filename = outname; 
    ncout.Format = "netcdf4";
    
        % make sure no dimensions are unlimited
        
    dimnames = ncout.dimlist;
    for i=1:length(dimnames)
        d = ncout.getdim(dimnames(i));
        d.Unlimited = false;        
    end
        
        % now update the variable's compression settings
        
    nvars = length(varnames);
    
    do_compress = ~isempty(deflatelevel);
    
    for i=1:nvars

        v=ncout.getvar(varnames(i));
        
                % set or clear add_offset and scale_factor
        if (~isempty(add_offset) && ~isempty(add_offset(i)))
            v.putatt("add_offset", single(add_offset(i)));
            v.putatt("scale_factor", single(scale_factor(i)));
            v.Datatype = "uint16";
            v.FillValue = uint16(65535);
            fill_ix = find(strcmp("_FillValue",v.attlist), 1);
            if (~isempty(fill_ix))
%               v.Attributes(fill_ix)=[]; 
                v.putatt("_FillValue", v.FillValue);
            end
        else
            offset_ix = find(strcmp("add_offset",{v.Attributes.Name}), 1);
            if (~isempty(offset_ix))
                v.Attributes(offset_ix) = [];
                v.Datatype="single";
                v.FillValue=float(1e20);
            end
            scale_ix = find(strcmp("scale_factor",{v.Attributes.Name}), 1);
            if (~isempty(scale_ix))
                v.Attributes(scale_ix) = [];
            end                            
%             fill_ix = find(strcmp("_FillValue",v.attlist), 1);
%             if (~isempty(fill_ix))
%                 v.Attributes(fill_ix)=[]; 
%             end
        end
                
        if (~do_compress)
            v.Shuffle = false;
            v.ChunkSize = [];
            v.DeflateLevel = [];
        else
            v.Shuffle = shuffle(i);
            v.ChunkSize = chunksize(i,:); 
            v.DeflateLevel = deflatelevel(i);
        end
    end
end



function  [ncin, outname, overwrite, varnames, shuffle, chunksize, deflatelevel, calc_scaling, add_offset, scale_factor, verbose, maxmem_orig] = init_params(inname, outname, varargin)

    inname   = string(inname);
    outname  = string(outname);
    
                    % parse input
    p = inputParser;
    addParameter(p,"inname",               "",  @(s) isnetcdf(s));
    addParameter(p,"outname",      strings(0),  @(s) isstring(s));
    addParameter(p,"outdir",       strings(0),  @(s) isstring(s));
    addParameter(p,"overwrite",         false,  @(s) islogical(s) || s==1 || s==0);
    addParameter(p,"varnames",    strings(0,0), @(s) ischar(s) || isstring(s));
    addParameter(p,"shuffle",            false, @(s) islogical(s) || s==1 || s==0);
    addParameter(p,"chunksize",             [], @(s) isnumeric(s) && (isempty(s) || length(s)==3));
    addParameter(p,"chunktype",         "time", @(s) isstring(s) && any(strcmpi(s, ["time","latlon"])));
    addParameter(p,"deflatelevel",          6, @(s) (isstring(s) && strcmpi(s,"none")) || (isnumeric(s) && s >= 0 && s <=9));
    addParameter(p,"int_scaling",        false); %,        @(s) isempty(s) || (length(s)==1 && islogical(s) || s==1 || s==0) || (isnumeric(s) && length(s)==2));
    addParameter(p,"verbose",                0, @(s) isnumeric(s) || islogical(s));
    addParameter(p,"maxmem",                [], @(s) isnumeric(s) && s >= 0 && s <=250e9);
                           
    parse(p, "inname", inname, "outname", outname, varargin{:});
    
    inname      = p.Results.inname;
    outname     = p.Results.outname;
    outdir      = p.Results.outdir;
    overwrite   = p.Results.overwrite;
    varnames    = p.Results.varnames;
    shuffle     = p.Results.shuffle;
    chunksize   = p.Results.chunksize;
    chunktype   = p.Results.chunktype;
    deflatelevel= p.Results.deflatelevel;
    int_scaling = p.Results.int_scaling;
    verbose     = p.Results.verbose;
    maxmem      = p.Results.maxmem;
    
    ncin = ncdf(inname);
    ncin.loadvars(true);
    
    [mydir,fname, fext] = fileparts(inname);
    if (isempty(outname))
        outname = sprintf("%s.tllc%s", fname,fext);   % tllc extension for changing chunksize in time,lat,lon format.
    end
    
    ss=extractBefore(outname,2);
    if (~any(strcmp(ss,["/","."])))
        if (isempty(outdir))
            outdir = mydir; 
        end
        outname = fullfile(outdir, outname);
    end
    
    if (isfolder(outname))
        error("error:  outname %s is a folder.", outname);
    elseif (isfile(outname) && ~overwrite)
        error("error:  outname %s is a file.  To overwrite, use ""overwrite"", true", outname);
    end
    
            % get list of climate variables, or check if all varnames are climate variables.
    if (isempty(varnames))
        varnames = intersect(ncin.varlist(), climate_varnames());
    else
        nonvars = setdiff(varnames, std_varnames());
        if (~isempty(nonvars))
            error("error:  varname is not a climate variable: %s\n", nonvars);
        end
    end
    
    calc_scaling = false;
    add_offset = [];
    scale_factor = [];
    if (islogical(int_scaling))
        if(int_scaling)
            calc_scaling = true;
        end
    else
        add_offset = int_scaling(1);
        scale_factor = int_scaling(2);
    end
    
        % get chunking if not specified in the inputs.
        
    v=ncin.getvar(varnames(1));

    [latname, lonname, timename,~, latix, lonix, timix] = ncdf_get_llt_dimnames(v);    

    if (isstring(deflatelevel) && strcmpi(deflatelevel,"none"))
        deflatelevel = [];
        chunksize = [];
        shuffle = false;
    elseif (deflatelevel == 0)
        if (shuffle)
            fprintf("\n\tNOTE:  deflatelevel is zero, so setting shuffle to false\n");
            shuffle = false;
        end
    elseif (isempty(chunksize))
        
        vtime = ncin.getvar(timename);
        vlats = ncin.getvar(latname);
        vlons = ncin.getvar(lonname);

        ntimes = max(vtime.Size);
        nlats  = max(vlats.Size);
        nlons  = max(vlons.Size);
        
        if (strcmpi(chunktype,"time"))
            timechunk = ntimes;
            latchunk  = 1;
            lonchunk  = 1;
        else
            timechunk = 1;
            latchunk  = nlats;
            lonchunk  = nlons;
        end        
        chunksize([latix,lonix,timix]) = [latchunk, lonchunk, timechunk];
    end
    
    [~, on] = get_hostname();
    if (isempty(maxmem))
        if (on.neys)
            maxmem = 40*1024*1024*1024;
        elseif (on.icsf_kmac)
            maxmem = 40*1024*1024*1024;
        elseif (on.hpcc_system)
            maxmem = 2*1024*1024*1024;
        else
            maxmem = 4*1024*1024*1024;
        end
    elseif (on.neys)
        maxmem = min(maxmem, 120e9);
    elseif (on.icsf_kmac)
        maxmem = min(maxmem, 90e9);
    elseif (on.icsf_lmac)
        maxmem = min(maxmem, 10e9);
    end
    
    maxmem_orig = ncObj.maxmem(maxmem);
        
end
        


% function update_history(ncout, nncs, inname, varnames, out_calendar, in_calendar, new_history)
% 
%     histstr = new_history + sprintf("%s:  joining %d minifiles with names like ''%s''; variables: %s", mfilename, nncs, inname, join(varnames," "));
%     
%     if (calendar_length(out_calendar) ~= calendar_length(in_calendar))
%         histstr = histstr + sprintf("; calendar adjusted from %s to %s", in_calendar, out_calendar);
%     end
%     
%     ncdf_update_history(ncout, histstr, now());
% end
% 
%    
    
