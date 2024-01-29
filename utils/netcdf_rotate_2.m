function netcdf_rotate_2(fname, varargin)
% function netcdf_rotate(fname, outname, "varname", varnames, "exclude", excludevars, "order", new_orders, "overwrite", do_overwrite, "format", ncdf_format, "ChunkSize", ChunkSize, "FillValue", FillValue)
%
%   Required Parameters:
%       fname
%   Required, but can be []:  
%       outname       if missing, will insert "llt" or "tll" into input filename for the output filename
%   Optional kwd/value parameters:
%       "varname", varnames     Specify variables to be rotated -- generally, Tmax, tasmax, tmin, ..., pr, precip...
%                               If more than 1 variable, use an array of strings or cell array.
%                               If varnames is empty or not specified, no variables will be rotated, but other parameters, such as chunking and
%                               deflation will be set, as specified by ChunkSize and DeflateLevel (see below).
%               NOTE:  for cmip6 files, must specify varname.  Currently can't parse filenames for variable names of cmip6 files.
%
%       "exclude", excludevars  if present, is list of variables to exclude from output.  Otherwise, all variables are
%                                   copied to the output file.
%       "order", new_orders     "llt" (or ["time","lat","lon"]) or "tll" (or ["lon","lat","time"]), or a string array.  If
%                               string array, 1st dimension should be the one you want to vary the fastest (contiguous)
%                               default output ordering is "llt".  (extension labelling is netcdf-terminology for the ordering (depth, column row), 
%                               whereas specifying order as an array, as in ["time","lat","lon"], is the order used by matlab
%                               and fortran (row, column) or (row, column, depth),
%                               You can specify 1 order to use for all varnames, or 1 order for each varnames that are
%                               being rotated.  If so, use a cell array of orders.
%                               If different orders for each variable being rotated, use cell array, 1 for each variable
%                                   NOTE: ordering of initials for "llt" or "tll" is based on netcdf naming convention,
%                                         which is reverse of Matlab's ordering.
%                                         Netcdf uses C ordering, (last dimension varies most rapidly
%                                         Matlab uses Fortran ordering (first dimension varies most rapidly)
%       overwrite           true/false.  Default is false
%       format             'classic' '64bit', 'netcdf4_classic' or 'netcdf4'.  
%                               Use netcdf4 if original file is classic (netcdf3) and you want to change it to netcdf4.
%                               NOTE:  netcdf4 assigns a default values for several parameters
%                               that are not automatically set in netcdf3 files:
%                                   ChunkSize  (some not-obvious function in the netcdf4 library...)
%                                   _FillValue (9.9...e36)
%                                   (and possibly some others)
%       FillValue          If rotated variables do not have a FillValue defined, I recommend you set this to 1e20 
%                               (which will be converted to single as needed).  This is the most common value used for
%                               double (or single) precision for climate netcdf files.
%                               If different FillValue for each variable to rotate, use cell array, one for each variable
%       maxmem             maxmem size for reading variables as a whole.  If variable size > maxread, then variable
%                               will be read and written in parts.  This will be much slower if the variable is being
%                               rotated.
%                                   
%                                   The following apply to all climate variables (tasmax, tasmin, pr, etc.)
%                                   If more than 1 variable, single entry applies to all climate variables
%                                   If separate chunking or deflation needed on each variable, specify as cell arrays,
%                                       1 cell for each variable. 
%
%       ChunkSize            ChunkSize to use for variables being rotated.  If blank, and netcdf format is netcdf4, the 
%                               netcdf4 library will assign default values, which may not be good for your purposes.
%                               can also be "CONTIGUOUS",  "DEFAULT"  or "PASS_THROUGH"
%                                       "PASS_THROUGH" = "rotate any existing ChunkSize, use defaults for any without
%                                                         ChunkSize"
%                                   NOTE:  CONTIGUOUS doesn't seem to be working in Matlab...
%                               or can be vector of numeric values, [5,10,20];
%                               Use inf to chunk along the full extent of a single dimension.
%                               [inf, 1 1] = chunk along time dimension.  [1, inf, inf] = chunk along lat/lon slices.
%                               For downscaling, set the chunksize to #days (or inf) for time axis, and 1 for lat & lon.
%                               For sharing with outside groups, verify whether they want it chunked along the time axis
%                               or along lat & lon (1 for time dimension, inf for lat & lon), which is what the 
%                               raw downloaded models are using in CMIP6.
%                               If different ChunkSize for each variable to rotate, use cell array, one for each variable
%       DeflateLevel       DeflateLevel to use for variables being rotated.  Use [] or an integer value between 0 & 9. 
%                               To use existing DeflateLevel, set DeflateLevel to []; 
%                           it by parts.
%       verbose             true/false/2
%
%   NOTE:  specify new_orders and ChunkSize in the order that the data should be written.
%           this is reverse of how netcdf displays the dimensions for the variable.
%           for example:  if data should be contiguous by time first, then lat, and finally by lon, then
%           use ["time","lat","lon"] for new_order.  (This is matlab ordering)
%           An ncdump of the output file will show the dimensions as "lon","lat","time".
%           Default output file name uses netcdf ordering.  It will append a pre-extension of "llt" for this.
%
%   NEEDED:
%       GROUPS:  doesn't handle variables in subgroups at present.  This will need some thought.
%       MEMORY:  this uses up to 1.25x the size of the input netcdf file in memory.  Should have a maxmem parameter and
%                break up variables if needed.
%                Entire file should fit comfortably in physical memory, or this program will run VERY SLOWLY.
%                local ncdf classes use low-level netcdf routines to read and write variables in pieces to avoid
%                Matlab's own high-level netcdf routines which can use up to 3x the memory when reading/writing data.
%
%   ChunkSize:
%   This should have an option for setting the ChunkSize for each rotated variable, rather than just DEFAULT, CONTIGUOUS
%   or PASS_THROUGH.  (done?)
%
%   finish and make use of:
%       update_variable_history
%       update_schema
%
%
%   NOTE:  chunksize and deflatelevel stuff removed for now.  Need to add that back in and test...
%
%   12/2/2021:  changed fix_unlimited_dimensions(...) to clear_unlimited_dimensions(...)
%
    t1 = tic();

%   [fname, nc, outname, varnames, excludevars, climvars, new_orders, do_overwrite, ChunkSize, DeflateLevel, ncformat, FillValue, verbose, maxmem_orig] = init_params(fname, varargin{:});
    [fname, nc, outname, varnames, excludevars,        ~, new_orders, do_overwrite,         ~,            ~, ncformat,         ~, verbose, maxmem_orig] = init_params(fname, varargin{:});
    
    fprintf("reading from:       %s\n", fname);
    fprintf("writing to:         %s\n", outname);
    fprintf("rotating variables: %s\n", vec2string(varnames,'brackets','[]'));

%   nc=ncdf(fname,'create',false);
    nc.loadvars([],true);   % load dimension variables only.
    
    if (~isempty_s(ncformat))      % set format to classic or netcdf4. 
        oldformat = nc.Format;
        nc.Format = ncformat; 
    else
        oldformat=[];
    end   
    nvars = length(varnames);
    
    
        % make sure rotation variables are found in input file
    ncvarnames = string({nc.Variables.Name});
    for i=1:nvars
        if (~any(strcmp(ncvarnames, varnames(i))))
            error("error:  variable %s not found in input file %s", varnames(i), fname);
        end
    end

        % create output ncdf
    ncout = nc.clone();
    ncout.Filename = outname;
    update_history(ncout, varnames, oldformat, ncformat)
    
    for j=1:nvars
        v=ncout.get(varnames(j));
        v.rotate_dims(new_orders{j},false);  % change dimension ordering.  false -> don't rotate actual variable, but we haven't read the actual variable data yet anyway.
%         for k=1:length(v.Dimensions)-1
%             v.Dimensions(k).Unlimited = false;
%         end
%         v.Size(end)=0;                       % set up for writing variable by parts later.
%         v.Dimensions(end).Unlimited = true;
    end
    
        % turn off Unlimited on all dimensions except new final dimension.
  
    ncout.clear_unlimited_dimensions();     % need to add back code to set selected dimension(s) to unlimited.  For most climate files, time should be unlimited.
%     for i=1:length(ncout.Dimensions)
%         if (any(strcmp(ncout.Dimensions(i).Name, unlim_dims)))
%             ncout.Dimensions(i).Unlimited = 1;
%         else
%             ncout.Dimensions(i).Unlimited = 0;
%         end
%     end
    
        % remove any excluded variables
    for j=1:length(excludevars)
        vix = ncout.varix(excludevars(j));
        if (isempty(vix))
            fprintf("warning:  exclude-var %s not found in file\n", excludevars(j));
        else
            ncout.Variables(vix) = [];
        end
    end    
    
    fprintf("writing schema...");
    ncout.writeschema(do_overwrite, "Dimensions",true);   % create file, but only write out the dimensions.
    fprintf("done\n");
    
    outvarlist = ncout.varlist();
    nncvars = length(outvarlist);
    nrots = 0;
    for i=1:nncvars
        vname = outvarlist(i);
        vin  = nc.getvar(vname);
        vndims = sum(vin.Size>1);
        
        if (vndims == 1 || vin.total_nbytes < ncObj.maxmem)
%           rotated = copy_var_whole(   nc, ncout, vin, vname, varnames, i, ChunkSize, DeflateLevel, FillValue, new_orders, verbose);
            rotated = copy_var_whole(   nc, ncout, vin, vname, varnames, i,                                     new_orders, verbose);
        else
%           rotated = copy_var_by_parts(nc, ncout, vin, vname, varnames, i, ChunkSize, DeflateLevel, FillValue, new_orders, verbose);
            rotated = copy_var_by_parts(nc, ncout, vin, vname, varnames, i,                                     new_orders, verbose);
        end
        nrots = nrots + rotated;
    end
    elapsed = toc(t1);
    fprintf("\nDone rotating.  Elapsed time:  %s\n", datestr(elapsed/86400,"HH:MM:SS"));
    fprintf("Done rotating:      %s\n", fname);
    fprintf("output written to:  %s\n", outname);
    fprintf("%d variables rotated\n", nrots);
    ncObj.maxmem(maxmem_orig);
end

function rotated = copy_var_whole(nc, ncout, vin, vname, varnames, i, new_orders, verbose)

            vout = ncout.getvar(vname);
%           toType = cast_as(vout.Datatype);
%           vout.vdata = cast(nc.readvar(vout.Name), toType);

            fprintf("reading variable %2d %-20s:  %12d bytes, size: %s\n",  i, vin.Name, vin.total_nbytes(), vec2string(vin.Size,"brackets",'[]'));
            vout.vdata = nc.readvar(vname);
    %       nc.loadvar(vname);

                % look for variable in list of varnames to rotate
            j=find(strcmp(varnames, vname),1);
            if (~isempty(j))
                if (all(strcmp(vin.dimlist(), new_orders{j})))
                    fprintf("NOTE:  variable order not changed!");
                    rotated = false;
                else
                    rotated = true;   % flag variable to be rotated.
                end
                permute_order = find_permute_order(vin, new_orders{j});
            else
                permute_order = [];
                rotated = false;
            end

                % ???? fix ChunkSize if not PASS_THROUGH
            if (verbose > 0)
                disp(vout);
                if (verbose > 1)
                    for j=1:length(vout.Attributes)
                        disp(vout.Attributes(j))
                    end
                end
            end
            
            if (rotated)
                vout.vdata = permute(vout.vdata, permute_order);
            end            
            if (isempty(permute_order))
                fprintf("writing variable %2d %-20s:  %12d bytes, size: %s\n",  i, vin.Name, vin.total_nbytes(), vec2string(vin.Size,"brackets",'[]'));
            else
                fprintf("writing variable %2d %-20s:  %12d bytes, size: %s\n",  i, vin.Name, vin.total_nbytes(), vec2string(vin.Size(permute_order),"brackets",'[]'));
            end
%             ncout.writevar(vname, [],[],[], permute_order);
            ncout.writevar(vname);
            pause(.025);
end

% function update_schema(ncout, FillValue, ChunkSize, DeflateLevel, new_orders)
%         if (~isempty(FillValue{j}))
%             vout.put("/Attributes/_FillValue", FillValue{j});
%         end
%         update_variable_history(vout, new_orders{j}, ChunkSize{j}, FillValue{j});
%         if (strcmp(ChunkSize{j},"DEFAULT"))
%             if (~isempty(vout.ChunkSize)), fprintf("chunksize changed from %s to %s\n", vec2string(vout.ChunkSize,"brackets",'[]'), ChunkSize{j}); end                
%             vout.ChunkSize = "DEFAULT";
%         elseif (strcmp(ChunkSize{j},"CONTIGUOUS"))
%             fprintf("chunksize changed from %s to %s\n", vec2string(vout.ChunkSize,"brackets",'[]'), ChunkSize{j});                
%             vout.ChunkSize = "CONTIGUOUS";
%         end
%         
%         
%                 orig_order = string({vin.Dimensions.Name});
% 
%                fprintf("rotating %s from %s to %s (orig chunk size: %s )\n", vout.Name, vec2string(orig_order,"brackets",'[]'), vec2string(new_orders{j},"brackets",'[]'), vec2string(vout.ChunkSize,"brackets",'[]'));
%                update_variable_history(vout, new_orders{j}, ChunkSize{j}, DeflateLevel{j}, FillValue{j});
%                if (all(strcmp(vin.dimlist(), new_orders{j})))
%                     fprintf("NOTE:  variable order not changed!");
%                     rotated = false;
%                else
%                     rotated = true;   % flag variable to be rotated.
%                end
%                
%                if (~isempty(FillValue{j}))
%                     vout.put("/Attributes/_FillValue", FillValue{j});
%                 end
% 
%                 if (strcmp(ChunkSize{j},"DEFAULT"))
%                     if (~isempty(vout.ChunkSize)), fprintf("chunksize changed from %s to %s\n", vec2string(vout.ChunkSize,"brackets",'[]'), ChunkSize{j}); end                
%                     vout.ChunkSize = "DEFAULT";
%                 elseif (strcmp(ChunkSize{j},"CONTIGUOUS"))
%                     fprintf("chunksize changed from %s to %s\n", vec2string(vout.ChunkSize,"brackets",'[]'), ChunkSize{j});                
%                     vout.ChunkSize = "CONTIGUOUS";
%                 end
% 
% end

function rotated = copy_var_by_parts(nc, ncout, vin, vname, varnames, i, new_orders, verbose)

            vout = ncout.getvar(vname);
%           toType = cast_as(vout.Datatype);
%           vout.vdata = cast(nc.readvar(vout.Name), toType);

                % look for variable in list of varnames to rotate
            j=find(strcmp(varnames, vname),1);
            if (~isempty(j))
                orig_order = string({vin.Dimensions.Name});

                fprintf("rotating %s from %s to %s (orig chunk size: %s )\n", vout.Name, vec2string(orig_order,"brackets",'[]'), vec2string(new_orders{j},"brackets",'[]'), vec2string(vout.ChunkSize,"brackets",'[]'));
                if (all(strcmp(vin.dimlist(), new_orders{j})))
                    fprintf("NOTE:  variable order not changed!");
                    rotated = false;
                else
                    rotated = true;   % flag variable to be rotated.
                end
                permute_order = find_permute_order(vin, new_orders{j});
            else
                permute_order = [];
                rotated = false;
            end

                % fix ChunkSize if not PASS_THROUGH
            if (verbose > 0)
                disp(vout);
                if (verbose > 1)
                    for j=1:length(vout.Attributes)
                        disp(vout.Attributes(j))
                    end
                end
            end
            
            fprintf("reading variable %2d %-20s:  %12d bytes, size: %s\n",  i, vin.Name, vin.total_nbytes(), vec2string(vin.Size,"brackets",'[]'));
            fprintf("writing variable %2d %-20s:  %12d bytes, size: %s\n",  i, vin.Name, vin.total_nbytes(), vec2string(vin.Size(permute_order),"brackets",'[]'));
            
            copy_by_parts(nc, ncout, vin, vout, vname, permute_order)
            pause(.025);
end

function copy_by_parts(ncin, ncout, vin, vout, vname, permute_order)

    if (~isempty(vout.ChunkSize)), error("can't handle chunksize for large variables yet.  sorry"); end
    
%    Report NaN's!!!


%     if (any(strcmp(["single","int32","uint32"], vin.Datatype)))
%         nb=4;
%     elseif (any(strcmp(["int16","uint16"], vin.Datatype)))
%         nb=2;
%     elseif (any(strcmp(["int8","uint8","char"], vin.Datatype)))
%         nb=1;
%     else
%         nb=8;
%     end

    nb=8;       % we have to assume 8 bytes because Matlab's netcdf library reads everything initially as doubles.
        
    isize = vin.Size;
    if (length(isize) == 2)
        nr = isize(1);
        nc = isize(2);
        
        stepsize = floor(maxmem/nb/nr);
        nsteps = ceil(nc/stepsize);
        laststep = nc - (nsteps-1)*stepsize;
        for i=1:nsteps
            istart = 1 + (i-1)*stepsize;
            if (i ~= nsteps)                
                d = ncin.readvar(vname, [1,istart],[nr,stepsize]);
            else
                d = ncin.readvar(vname, [1,istart],[nr,laststep]);
            end
            if (isempty(permute_order))
                ncout.writevar(vname, d, [1,istart]);
            else
                ncout.writevar(vname, d', [istart,1]);
            end
            if (ncobj.verbose), show_progress(i,nsteps); end
        end
    elseif (length(isize) == 3)
        
        nr=isize(1);
        nc=isize(2);
        nz=isize(3);
        
        stepsize = floor(ncObj.maxmem/nb/nr/nc);
        nsteps = ceil(nz/stepsize);
        laststep = nz - (nsteps-1)*stepsize;
        if (ncObj.verbose)
            fprintf("copying %s in %d steps:",vin.Name, nsteps); 
        end
        for i=1:nsteps
            istart = 1 + (i-1)*stepsize;            
            if (i ~= nsteps)
                d = ncin.readvar(vname, [1,1,istart], [nr,nc,stepsize]);
            else
                d = ncin.readvar(vname, [1,1,istart], [nr,nc,laststep]);
            end
            fprintf("r%d.", i);
            if (isempty(permute_order))
                ncout.writevar(d, vname, [1,1,istart]);
            else
                d = permute(d, permute_order);
                if (permute_order(1) == 3)
                    start = [istart,1,1];
                elseif (permute_order(2) == 3)
                    start = [1,istart,1];
                else
                    start = [1,1,istart];
                end
                ncout.writevar(vname, d, start);
            end
            if (ncObj.verbose), fprintf("w%d.", i); end
        end
        fprintf("\n");
    else
        error("can't handle large variables with more than 3 dimensions yet.  sorry");        
    end
end


function permute_order = find_permute_order(v, ordstrings)

    dimnames = string({v.Dimensions.Name});
    nords = length(ordstrings);
    permute_order = zeros(1,nords);
    for i=1:nords
        pord = find(dimnames==ordstrings(i));
        if (isempty(pord) || length(pord) > 1), error("error finding permute order for variable %s", v.Name); end
        permute_order(i)=pord;
    end
end

% function toType = cast_as(dtype)
%     shortstuff = ["single","int8","uint8","int16","uint16","int32","uint32"];
%     if (any(strcmp(shortstuff, dtype)))
%         toType = "single";
%     else
%         toType = dtype;
%     end
% end

function update_history(nc, varnames, oldformat, ncformat)

    hist = sprintf("netcdf_rotate: dimension re-ordering: ");

    hist = hist + sprintf("%s ", varnames);

    if (exist('ncformat','var') && ~isempty(ncformat))
        hist = hist + sprintf("; changed netcdf format from %s to %s", oldformat, ncformat);
    end
    ncdf_update_history(nc, hist); 
end

% function update_variable_history(obj, new_order, ChunkSize, FillValue)
%     hist = sprintf("netcdf_rotate: %s: dimension re-ordering: ", datestr(now,"yyyy-mm-dd HH:MM"));
% 
%     hist = hist + sprintf("from [%s] to [%s]; ", vec2string({obj.Dimensions.Name},"brackets",'[]'), vec2string(new_order));
%     if (~strcmp("ChunkSize","PASS_THROUGH") && ~isempty(obj.ChunkSize))
%         hist = hist + sprintf("%s; ChunkSize changed from %s to %s; ", vec2string(obj.ChunkSize,"brackets",'[]'), vec2string(ChunkSize));
%     end
%     if (exist('newFillValue','var') && ~isempty_s(FillValue))
%         if (isempty(obj.FillValue))
%             oldFillValue = "[]";
%         else
%             oldFillValue = sprintf("%g", obj.FillValue);
%         end
%         hist = hist + sprintf("FillValue changed from %s to %g; ", oldFillValue, FillValue);
%     end
%     try 
%         h = obj.get("/Attributes/history");
%         history = sprintf("%s\n", string(h.Value));
%     catch
%         history = "";
%     end
%     
%     history = history + hist;
%     
%     obj.put("/Attributes/history", history);
% end

function [fname, nc, outname, varnames, excludevars, climvars, new_orders, overwrite, ChunkSize, DeflateLevel, ncformat, FillValue, verbose, maxmem_orig] = init_params(fname, varargin)

    if (isempty(varargin))  % in case
        outname=[]; 
    else
        outname=string(varargin{1});
        varargin = varargin(2:end);
    end
    if (~isempty_s(outname) && any(strcmp(["varnames","new_orders","do_overwrite","ChunkSize","format","FillValue"], outname)))
        error("error:  please specify outname.  Use [] to generate outname from input filename");
    end
                    % parse input
    p = inputParser;

                    % these are the params we want to handle in ARRM_V2_wrapper
    addParameter(p,"varnames", [],  @(s) isempty(s) || ischars(s));
    addParameter(p,"exclude", [],  @(s) isempty(s) || ischars(s));
    addParameter(p,"order", [],  @(s) isempty(s) || ischars(s));
    addParameter(p,"overwrite", false,  @(s) islogical(s) || (isnumeric(s) && (s==0 || s==1)));
    addParameter(p,"ChunkSize", [], @(s) isempty_s(s) || ischar_s(s) || isnumeric(s));     % DEFAULT, CONTIGUOUS, or "PASS_THROUGH" (use original input, rotated)).
    addParameter(p,"DeflateLevel",[], @(s) isempty(s) || (s>=0 && s<=9));
    addParameter(p,"format", [],  @(s) isempty(s) || ischars(s));
    addParameter(p,"FillValue", [],  @(s) isempty(s) || isnumeric(s));
    addParameter(p,"verbose", 0,  @(s) isnumeric(s) || islogical(s));
    addParameter(p,"maxmem", 16*1024*1024*1024,  @(s) isnumeric(s) || islogical(s));      % default to 16 GB reads max.
        
    parse(p, varargin{:});
    
    varnames     = string(p.Results.varnames);
    excludevars  = string(p.Results.exclude);
    new_orders   = string(p.Results.order);
    overwrite    = p.Results.overwrite;
    ChunkSize    = p.Results.ChunkSize;
    DeflateLevel = p.Results.DeflateLevel;
    ncformat     = p.Results.format;
    FillValue    = p.Results.FillValue;
    verbose      = 1*p.Results.verbose;
    maxmem       = p.Results.maxmem;
    
    ncObj.verbose(true);    
    
    maxmem_orig = ncObj.maxmem(maxmem);
                
    nc=ncdf(fname,'create',false);
    
    climvars = intersect(nc.varlist, climate_varnames());        % get list of climate variables in the file (pr, tasmin, tasmax, etc.)
    
    nc.clear_unlimited_dimensions();     % need to add code to set selected dimension(s) to unlimited.  For most climate files, time should be unlimited.
   
        % if varname not specified, extract it from the input filename
%       fix parse_netcdf_filenames for cmip6 files, ian
    if (isempty_s(varnames))
        pinfo = ARRM_V2_parse_netcdf_filenames(fname);
        varnames = string(pinfo.varname);     
        if (isempty(varnames))
            varnames = find_climate_varname(nc.varlist, "all");
        end
        if (isempty_s(varnames)), error("error:  cannot get varname from filename %s", fname); end
    end
    nvars = length(varnames);
    
    if (isempty_s(new_orders)) 
        new_orders = {["time","lat","lon"]}; 
    elseif (~iscell(new_orders))
        new_orders = {new_orders};
    end
    
    nords = length(new_orders);
    for i=1:nvars
        if (i > nords)
            new_orders{i} = new_orders{i-1};
        else
            if (strcmpi(new_orders{i},"llt"))
                new_orders{i} = ["time","lat","lon"];
            elseif (strcmpi(new_orders{i},"tll"))
                new_orders{i} = ["lon","lat","time"];
            end
        end

            % get actual dimension names for new_order
        v=nc.getvar(varnames(i));
        for j=1:length(new_orders{i})
            [~,~,~,dimname] = v.diminfo(new_orders{i}(j));
            if (isempty(dimname)), error("error:  Cannot find matching dimension for %s in variable %s, which has dimensions %s", new_orders{i}(j), varnames(i), vec2string(v.dimlist(),'brackets','[]')); end
            new_orders{i}(j) = string(dimname);
        end
        
    end
    
    if (~iscell(ChunkSize)), ChunkSize = {ChunkSize}; end
    nchunks = length(ChunkSize);
    for i=1:nvars
        if (i > nchunks), ChunkSize{i} = ChunkSize{i-1}; end
        if (ischar_s(ChunkSize{i}))
            ChunkSize{i} = upper(ChunkSize{i}); 
            if (~any(strcmp(ChunkSize{i}, {'DEFAULT','CONTIGUOUS','PASS_THROUGH'})))
                error("error:  bad ChunkSize: %s\n", ChunkSize{i});
            end
        end
    end
    
    if (~iscell(FillValue)), FillValue = {FillValue}; end
    nfills = length(FillValue);
    for i=nfills+1:nvars
        FillValue{i} = FillValue{i-1};
    end
    
        % if outname not specified, create sub-extension like ".llt" and insert it into fname to make the output
        % filename
    if (isempty_s(outname))
        ndims = length(new_orders{1});
        subext = ".";
        for i=ndims:-1:1
            subext = strcat(subext,extractBefore(new_orders{1}(i),2));  % make "llt" or "tll" or whatever.
        end
        [path,fn,fext] = fileparts(fname);
        outname = fullfile(path, sprintf("%s%s%s", fn, subext, fext)); % and assemble outname
    end
    
    
    if (~overwrite && isfile(outname))
        error("error:  file %s already exists.  Set ""overwrite"" to true to overwrite", outname);
    end
    

end