function [outname, logname] = join_minifiles( fns, varnames, outname, varargin )
% function join_minifiles( fn_template, varnames, fnout, latrange, lonrange, ...)
%
%   This was based on my old code, updated to use my ncdf classes.  Joins a set of minifiles into 1 netcdf file.
%   Allows for multiple data variable in the netcdf files.
%   Assumes that lat & lon info is present in the files.
%   Assumes that metadata and is consistent across all files
%   Assumes all dimensions are in ascending order.
%   Assumes sufficient memory to read each minifile's data into memory at once.
%   finds all filenames matching the template, extracts data from minifiles and writes out file fnout
%   If latrange specified (lat1, lon1, lat2, lon2), then only extracts lats & lons in specified range.
%
%   Important assumptions:
%       1.  That data in each file is contiguous (no gaps in time or space), and no gaps across file boundaries
%               e.g., if joining along time, there are no gaps in time, such as one file ending at 11/30/1949 and the
%               next file starting at 1/1/1950 (as in total output missing month of december)
%       2.  That dimension data is identically valued in each file (except along joining dimensions)
%               no rounding differences
%               no difference in increments, such as one file @ 1.5 degree steps, another at 1.0 degree steps.
%               identical time units in all files, as in "days since 1950-01-01"
%       4.  That lat,lon & time dimensions have obvious names (for limiting output ranges)
%               possible values: lat, lats, latitude, latitudes, lon, lons, longitude, longitudes, time
%                   (code ignores case... i.e.,  Lat or LAT also acceptable)
%       5.  netcdf4 Groups contain only metadata, such as run parameters.  
%               NOTE:  metadata from 1st file only will be retained.
%
%   NOTE:  when joining Anne's minifiles, it works fastest if you rotate the output to LLT format.
%          when joining downloaded netcdf files that are broken out by year, it is probably faster to join them in
%          standard order, then rotate the file separately with netcdf_rotate.m.
%
%   Inputs:
%       fns                                     filename template or string array of filenames
%       varnames                                list of varnames to copy
%                                                   Will automatically include all dimensions used by varnames as well.
%       outname                                 output filename
%     optional args:
%     key/value args:
%       'latrange', [lat1,lat2]                 lat range to include in output
%       'lonrange', [lon1,lon2]                 lon range to include in output
%                                                   note:  actual range depends on input.  If entire range is not
%                                                   covered by data, only actual range is output.
%       'daterange',[yyyy,mm,dd; yyyy,mm,dd]    date range to include in output.  can include time: yyyy,mm,dd,hh,mm,ss;
%                                                   NOTE:  include hh,mm,ss if input files have non-integer times (such as  xxx.5 %                                                   
%       'excludes', ["var1","var2"....]         list of variables to exclude
%       'overwrite', tf                         true/false.  true if OK to overwrite existing.  [false];
%       'format', ncdf_fmt                      'classic','64bit','netcdf4_classic','netcdf4'.  Default: "netcdf4" 
%       'keep_dims', tf                         true/false.  If true, retains all dimension variables, even if not used by any retained variable
%       'history','text'                        text to be prepended to existing history from first file
%       'verbose', true/false                   0:  quiet;  1: show progress bars;  2:  output info on each file as read it is processed.
%       'rotate', true/false                    rotate axes (moves time from last dimension to first, or first dimension
%                                                   to last.  Applies to lon/lat/time dimensions and lon/lat/day-of-year
%       'chunksize',[s1,s2,s3]                  set the chunksize to [s1,s2,s3].  Default:  []
%                                                   order of chunksize values should match the order in the input file (independent of rotate setting)
%                                                   if anydimension's  chunksize is set to inf, it will be replaced by the output size for that dimension.
% %       'rename_zvals'                          if true, renames zvals or packed_zvals to varnames(1)_zvals.  default: [true] 
% %                                                   i.e., if varname(1) is tasmax, then zvals -> tasmax_zvals.
%       'unlimited', true/false                 set time's Unlimited flag.  If empty (default), use Unlimited from input files.
%       'checkonly', true/false                 if true, only checks for files and missing dates.
%       'completeonly", true/false              if true, only joins if there are no gaps in the data.
%       'logname', logname                      specify logfile to write info about run.
%                                                   If specified, appends to file if it exists, or creates it if otherwise. 
%                                                   If not specified, logfile is same-name as output nc file, with name of "..._log.txt"
%       "include_zvals", true/false             if true, looks for variable w/ name ending in "zvals", includes it, and
%                                                   renames it to (varname)_zvals.  This handles legacy cases where
%                                                   zvals were named either "zvals" or "packed_zvals".
%                                                   NOTE:   alternately, set to false and add zvals' actual variable name to list of
%                                                           varnames to copy that variable without renaming it.
%       "add_source_metadata", true/false       if true, looks for a subgroup called "Source_metadata".  If that group
%                                                   exists but only contains the source filenames, load the metadata
%                                                   from the files as a subgroup.
%                                               if false, simply copies the Source_metadata group as is (whether or not
%                                               it contains the source file's metadata.)
%       "reorder_vars", [true]/false            if true, sorts variables by size in the output file (helps with access speed) 
%                                                   
%
% % xxxx    'calendar', calendar                    calendar type.  [365-day]  (removed.  now defaults to calendar type of
%                                                   first input file, and warns if calendar type differs for any of the
%                                                   input files.
%
%
%   NOTE:  MOST OF THE NON-DEFAULT OPTIONS HAVEN'T BEEN FULLY TESTED!  
%   NOT TESTED YET:
%       * Combining files with different time slices, rather than different lat/lon regions.
%       * Combining lat/lon regions that are not 1 longitude wide.
%       * limiting the lats or lons or time
%       * changing the calendar --- was buggy, fixed 12/7/2019 icsf
%       * keeping limited list of variables.
%       * excluding variables
%
%
%   12/2/2021:  changed fix_unlimited_dimensions(...) to clear_unlimited_dimensions(...)
%   12/2/2021:  added option to set chunking.
%
%--------------------------

%   addpath("Users/iscottfl/Documents/MATLAB",	"/Users/iscottfl/Desktop/atmos2/ARRM_V2".	"/Users/iscottfl/Desktop/atmos2/ARRM_V2/util_general",	"/Users/iscottfl/Desktop/atmos2/ARRM_V2/util_nc",	"/Users/iscottfl/Desktop/atmos2/ARRM_V2/ncdf",	"/Users/iscottfl/Desktop/atmos2/ARRM_V2/ARRM_V2_subs",	"/Users/iscottfl/Desktop/atmos2/ARRM_V2/utils"
    addpath("./util_general",	"./util_nc",	"./ncdf",	"./ARRM_V2_subs",	"./utils");

        t1 = tic();

%       [ncs, keeper_varnames, unneeded, excludes, outname, latrange, lonrange, daterange, outlats, outlons, calendar, do_overwrite, ncformat, verbose, new_history, do_rotate, chunksize, unlimited, do_clean, fillvalue, checkonly, logname, completeonly, add_source_metadata, reorder_vars  = init_params(fns, varnames, outname, varargin{:});
        [ncs, keeper_varnames, unneeded, excludes, outname,     ~,      ~,      daterange, outlats, outlons, out_calendar, do_overwrite, ncformat, verbose, new_history, do_rotate, chunksize, unlimited, do_clean, fillvalue, checkonly, logname, completeonly, add_source_metadata, reorder_vars] = init_params(fns, varnames, outname, varargin{:});


        [nc_time_info, gaps] = check_daterange(ncs, out_calendar, daterange, outname, logname, verbose);

        if (checkonly || (completeonly && ~isempty(gaps)))
            log_print(logname, 1, "\n\njoin_minifiles:  done checking dates.  Log file is:  %s\n", logname);
            return; 
        end
        
            % scan through all ncs, get complete list of dimensions and their values.
            % create ncout and update the dimension info in ncout to cover the full range of values for each dimension.
            % reduces list of ncs to only those with data in the lat, lon and date ranges.
%       [ncout, ncs, permute_orders] = netcdf_init(ncs, fns, latrange, lonrange, daterange, outlats, outlons, keeper_varnames, unneeded, outname, nc_time_info, do_rotate, chunksize, unlimited, do_clean, ncformat, new_history, fillvalue, add_source_metadata, logname); 
        [ncout, ncs, permute_orders] = netcdf_init(ncs, fns,                                outlats, outlons, keeper_varnames, unneeded, outname, nc_time_info, do_rotate, chunksize, unlimited, do_clean, ncformat, new_history, fillvalue, add_source_metadata, logname); 

        fprintf("done netcdf_init\n"); pause(.01);
        [latname, lonname, timename] = ncdf_get_llt_dimnames(ncout);
        nncs = length(ncs);
        nvars = length(ncout.Variables);
        vnamelen = max(strlength(to_row(keeper_varnames)));
        lats = ncout.getvardata(latname);
        lons = ncout.getvardata(lonname);
        [~,~,~,~,~,tstamps] = ncdf_get_time_info(ncout, timename);
            

            % reorder variables so we put dimensions first, then the remaining variables, smallest to largest.
        if (reorder_vars)
            var_output_order = get_var_output_order(ncout);
            ncout.Variables = ncout.Variables(var_output_order);
        else
            var_output_order = 1:length(ncout.Variables);
        end

                % we'll get the data in chunks so we don't need too much memory.
        log_print(logname, 1, "\nreading %d files, %d lats, %d lons\n",  nncs, length(lats), length(lons));
        log_print(logname, 1, "writing output to:  %s\n", outname);
        log_print(logname, 1, "lat range:  ( %10.4f %10.4f )\n", lats(1), lats(end));
        log_print(logname, 1, "lon range:  ( %10.4f %10.4f )\n", lons(1), lons(end));
        log_print(logname, 1, "time range: ( %10s %10s )\n", datestr_cal(tstamps(1),out_calendar,"yyyy-mm-dd"), datestr_cal(tstamps(end),out_calendar,"yyyy-mm-dd"));
        if (~isempty_s(excludes))
            log_print(logname, 1, "excluding variables:  %s\n", vec2string(excludes),"brackets",'[]');
        end
        log_print(logname, 1, "including variables:  %s\n", vec2string(keeper_varnames(1,:),"brackets",'[]'));
        if (any(~strcmp(keeper_varnames(1,:), keeper_varnames(2,:))))
            log_print(logname, 1, "output as:            %s\n", vec2string(keeper_varnames(2,:),"brackets",'[]'));
        end

        fprintf("creating output file %s\n", outname);
        pause(.25);
        ncout.writeschema(do_overwrite, "filename", outname, "Dimensions",true);
        fprintf("done creating output file\n");    

        dlist = ncout.dimlist();

                % read each file and write it to the output file.
        for j=1:nncs
     %      vdata = get_var_data(ncs{j}, vname, out_calendar);
            for i=1:nvars
                    % load data for all variables that we're keeping.
                if (~any(strcmp(keeper_varnames(1,i), dlist)))       % if variable isn't in list of dimensions, load it from file.
                    ncs{j}.loadvar(keeper_varnames(1,i));
                end
            end
            [~,in_calendar] = ncdf_get_time_info(ncs{j});
            if (calendar_length(in_calendar) ~= calendar_length(out_calendar))
                ncdf_adjust_calendar(ncs{j},out_calendar);      % will adjust all time-dependent variables. 
                        % NOTE:  We shouldn't allow changing calendar when joining downloaded CMIP6 model data.
%               log_error("error:  mismatched calendar length:  %d for 1st file, but %d for %s\n", calendar_length(calendar), calendar_length(in_calendar), ncs{j}.Filename);
            end
            for ii=1:nvars
                ix=var_output_order(ii);
                vname = keeper_varnames(1,ix);
                outvname = keeper_varnames(2,ix);
                if (any(strcmp(vname, dlist))), continue; end     % skip dimensions;  we wrote them already.
    %           log_print(logname, 1, "writing variable %s, %d files\n", vname, nncs)

                v = ncs{j}.getvar(vname);
                if (do_rotate && ~isempty(permute_orders{ix}))
                    v.rotate_dims(permute_orders{ix});
    %                 vdata = permute(vdata, permute_orders{ix});
    %                 ncs{j}.Variables(ix).rotate_dims(permute_orders{ix}, true);
                end

                [start, startpos, sz] = save_results(ncout, ncs{j}, latname, lonname, timename, vname, outvname, logname);
                if (isempty(start) && verbose)
                        log_print(logname, 1, "%5d of %5d: %s: %s  %s\n", j,nncs,  vname, ncs{j}.Filename, "no data in range");
                elseif (verbose)
                    try
                        log_print(logname, 1, "%5d of %5d: %-*s size %s to start %s, %s from file %s\n", j,nncs,  vnamelen, outvname, vec2string(sz, 'format','%5d','brackets','[]'), ...
                                  vec2string(start, 'format','%5d','brackets','[]'), vec2string(startpos,'format',"%s", 'brackets','[]'), basename(ncs{j}.Filename));
                        %
                    catch
                        log_print(logname, 2, "oops!\n");
                    end
                end
            end
            if (~verbose),  show_progress(j,nncs); end
                %   clear out the variable to release memory.
            ncs{j} = [];

        end
        
%             % write NAs to all the gaps.
%             % This shouldn't be necessary, but Matlab doesn't seem to pay attention to the zzzzz parameter, which is set
%             % to fill the file with NAs at creation.
%
        if (~isempty(gaps))
            write_NAs(ncout, gaps, latname, lonname, timename, outvarlist, logname);
        end

        ncout.writeatt("completion_status",int32(1));  % change completion flag to complete from incomplete.   
        elapsed = toc(t1);
        log_print(logname, 1, "\nDone joining files.  Elapsed time:  %s\nOutput written to:  %s\n\n", datestr(elapsed/86400,"HH:MM:SS"), outname);
        
        clear ncs ncout
        pause(1);
       
end

% function vdata = get_var_data(nc, vname, out_calendar)
%     vdata = nc.readvar(vname);
%     v = nc.getvar(vname);   
%         % see if we're time-dependent.  If not, we're done.
%     [timix, ~, ~, timename] = v.diminfo('time');
%     if (isempty(timix)), return; end
% 
%     [~, calendar, svec, dsince] = ncdf_get_time_info(nc, timename);
%         % if we're not changing the calendar, we're done.
% 
%     if (calendar_length(out_calendar) ~= calendar_length(calendar))
%         [~,~,vdata] = ncadjust_calendar(svec, dsince, vdata, calendar, out_calendar, 2, [], timix);
%     end
%     
% end

function [start, startpos, sz] = save_results(ncout, nc, latname, lonname, timename, varname, outvarname, logname)

    % 1.  get ncout's lats, lons
    % 2.  find indexes of nc's lats, lons save
    % 3.  determine if lats, lons are contiguous in output file
    % 4.  if contiguous
    %       write data
    %     else
    %       for each lon
    %           for each lat
    %               write data
    %
    
    start    = [];
    sz       = [];
    startlon = [];
    startlat = [];
    starttim = [];
    
    v = nc.getvar(varname);
    data = v.vdata;
    vdimlist = v.dimlist();
    startpos = strings(1,length(vdimlist));    % (return value for printing in log file) 

    latix = find(strcmp(vdimlist, latname),1);
    lonix = find(strcmp(vdimlist, lonname),1);
    timix = find(strcmp(vdimlist, timename),1);

        % figure out what portion of the data to write out, and where to start writing it.
    if (~isempty(lonix))
        outlons = ncout.getvardata(lonname);    
        inlons = nc.getvardata(lonname);
        lonkeepers  = ismember(inlons, outlons);
        if (sum(lonkeepers) == 0)
            return;
        end
        lon1 = find(lonkeepers,1);
        lon2 = find(lonkeepers,1,'last');
        startlon = find(outlons == inlons(lon1),1);
        nlons = sum(lonkeepers);
        contiguous_lons = (lon2 - lon1 == nlons-1);
        if (~contiguous_lons), log_error(logname, 'can''t do non-contiguous lons write yet'); end
        if (sum(lonkeepers) ~= length(lonkeepers))
            data = extract_data(data, lon1, lon2, lonix, logname); 
        end
        startpos(lonix) = sprintf("%.4f", outlons(startlon));
    end

    if (~isempty(latix))
        outlats = ncout.getvardata(latname);
        inlats = nc.getvardata(latname);
        latkeepers  = ismember(inlats, outlats);
        if (sum(latkeepers)==0)
            return;
        end
        lat1 = find(latkeepers,1);
        lat2 = find(latkeepers,1,'last');
        startlat = find(outlats == inlats(lat1),1);
        nlats = sum(latkeepers);    
        contiguous_lats = (lat2 - lat1 == nlats-1);
        if (~contiguous_lats), log_error(logname, 'can''t do non-contiguous lats write yet'); end
        if (sum(latkeepers) ~= length(latkeepers))
            data = extract_data(data, lat1, lat2, latix, logname); 
        end
        startpos(latix) = sprintf("%.4f", outlats(startlat));
    end

    if (~isempty(timix))
        outstamps = get_tstamps(ncout);
        instamps  = get_tstamps(nc);
        timekeepers = ismember(instamps, outstamps);
        if (sum(timekeepers) == 0)
            return;
        end
        tim1 = find(timekeepers,1);
        tim2 = find(timekeepers,1,'last');
        starttim = find(outstamps == instamps(tim1),1);
        ntims = sum(timekeepers);
        contiguous_times = (tim2 - tim1 == ntims-1);
        if (~contiguous_times), log_error(logname, 'can''t do non-contiguous time write yet'); end
        if (sum(timekeepers) ~= length(timekeepers))
            data = extract_data(data, tim1, tim2, timix, logname); 
        end
        [~,calendar] = ncdf_get_time_info(ncout);        
        startpos(timix) = datestr_cal(outstamps(starttim),calendar, 'yyyy-mm-dd'); 
    end    

        % identify where it needs to go in the output file
    start = make_start(ndims(data), lonix, latix, timix, startlon, startlat, starttim);

        % and write it.
    ncout.writevar(outvarname, data, start);
    sz = size(data);
    if (length(sz) < length(startpos))  % in case last dimension(s) is 1;  matlab drops the last dimension if it is of size 1.
        nd=length(startpos);
        sz(end+1:nd)=1;
    end

        % and fill in any empty startpos fields.
    for i=1:length(startpos)
        if (strlength(startpos(i))==0)
            %startpos(i) = string(data(i,1)); 
            startpos(i) = "*";
        end
    end
end

function write_NAs(ncout, gaps, latname, lonname, timename, varlist, logname)

    ngaps = size(gaps,1);
    if (ngaps == 0), return; end
    [~,calendar] = ncdf_get_time_info(ncout);
    
    log_print(logname, 1, "writing NAs to fill in missing data:\n");
    for j=1:length(varlist)
        vname = varlist(j);
        if (is_climate_variable(vname))
            for i=1:ngaps
                [start, startpos, latix, lonix, timix, sz] = write_NA_gap(ncout, gaps(i,:), latname, lonname, timename, vname);
                log_print(logname, 1, "%5d of %5d: %s size (%d,%d,%d) to start %s, (%8.4f,%8.4f, %s)\n", i,ngaps,  vname, sz, vec2string(start, 'format','%5d','brackets','[]'), startpos(latix), startpos(lonix), datestr_cal(startpos(timix),calendar, 'yyyy-mm-dd'));
            end
        end
    end
end

function [start, startpos, latix, lonix, timix, sz] = write_NA_gap(ncout, gap, latname, lonname, timename, vname)
%
%   writes NAs to the output file for times specified by gap.
%   gap(1) is index (1-based) of beginning of gap
%   gap(2) is index           of end of gap
%
%   this is more complicated than expected because we have to get the dimension ordering from ncout.
        
    outlats = ncout.getvardata(latname);
    outlons = ncout.getvardata(lonname);
        
    outstamps = get_tstamps(ncout);
    
    v = ncout.getvar(vname);
    vdimlist = v.dimlist();
    
    latix = find(strcmp(vdimlist,latname));
    lonix = find(strcmp(vdimlist, lonname));
    timix = find(strcmp(vdimlist, timename));
    
    startpos = nan(1,3);
    
    startpos(lonix) = outlats(1);
    startpos(latix) = outlons(1);
    startpos(timix) = outstamps(gap(1));

    nlats = length(outlats);
    nlons = length(outlons);
    nnas  = gap(2)-gap(1)+1;
    
    sz = nan(1,3);
    sz([lonix,latix,timix]) = [nlons,nlats,nnas];
    
    nas = nan(sz);
    
            % identify where it needs to go in the output file
    start = make_start(ndims(v), lonix, latix, timix, 1, 1, gap(1));

            % and write it.
    ncout.writevar(vname, nas, start);
end

function start = make_start(nd, lonix, latix, timix, startlon, startlat, starttim)
% creates a netcdf 'start' vector, identifying where to start writing the data to the output file.
    start = ones(1,nd);
    if (~isempty(lonix)), start(lonix) = startlon; end
    if (~isempty(latix)), start(latix) = startlat; end
    if (~isempty(timix)), start(timix) = starttim; end
end

% function data = extract_data(data, lon1, lon2, lat1, lat2, tim1, tim2, lonix, latix, timix)
% %   returns subset of data specified by lon, lat and time ranges.
% %   If variable is not dependent on lon (or lat, or time), that range is not used.
%     if (~isempty(lonix) && ~use_all_lons),  data = extract_sub(data, lon1, lon2, lonix); end
%     if (~isempty(latix) && ~use_all_lats),  data = extract_sub(data, lat1, lat2, latix); end
%     if (~isempty(timix) && ~use_all_tims),  data = extract_sub(data, tim1, tim2, timix); end
% end

function data = extract_data(data, ix1, ix2, dimix, logname)
% %   returns subset of data specified by lon, lat and time ranges.
% %   If variable is not dependent on lon (or lat, or time), that range is not used.

    if (dimix == 1)
        data = data(ix1:ix2, :, :, :, :, :, :, :);
    elseif (dimix == 2)
        data = data(:, ix1:ix2, :, :, :, :, :, :);
    elseif (dimix == 3)
        data = data(:, :, ix1:ix2, :, :, :, :, :);
    elseif (dimix == 4)
        data = data(:, :, :, ix1:ix2, :, :, :, :);
    elseif (dimix == 5)
        data = data(:, :, :, :, ix1:ix2, :, :, :);
    elseif (dimix == 6)
        data = data(:, :, :, :, :, ix1:ix2, :, :);
    elseif (dimix == 7)
        data = data(:, :, :, :, :, :, ix1:ix2, :);
    elseif (dimix == 8)
        data = data(:, :, :, :, :, :, :, ix1:ix2);
    else
        log_error(logname, "error: extract_sub:  can't handle %d dimensions", dimix);        
    end
end

function tstamps = get_tstamps(nc)
    [~, cal, start_vec, days_since] = ncdf_get_time_info(nc);   % get calendar info
        tstamps = datenum_cal(start_vec, cal) + round(days_since-.1);       % get datenums in nc's calendar.  Using round(...) here to fix timestamps that are x.5
end

function   [ncs, keeper_varnames, unneeded, excludes, outname, latrange, lonrange, daterange, outlats, outlons, out_calendar, do_overwrite, ncfmt,    verbose, history,     do_rotate, chunksize, unlimited, do_clean, fillvalue, checkonly, logname, completeonly, add_source_metadata, reorder_vars] = init_params(fn_list, varnames, outname, varargin)

    varnames = string(varnames);
    outname  = string(outname);
                    % parse input
    p = inputParser;

    addOptional(p,"latrange",   [-90,90],       @(s) isnumeric(s) && s(1)>= -90 && s(2) >= s(1) && s(2) <= 90);
    addOptional(p,"lonrange",   [0,360],        @(s) isnumeric(s) && s(1)>=   0 && s(2) >= s(1) && s(2) <=360);
    addOptional(p,"daterange",  [],             @(s) isnumeric(s) && isempty(s) || (size(s,2)>=3 && size(s,1)==2 && datenum(s(1,:)) >= datenum(1850,1,1) && datenum(s(2)) <= datenum(2150,12,31) && datenum(s(1)) <= datenum(s(2))));
%   addParameter(p,"calendar",  "365-day",      @(s) ischar_s(s));  %   changing calendar no longer allowed.  This is too complicated to do here. 
    addParameter(p,"overwrite", false,          @(s) islogical(s) || s==1 || s==0);
    addParameter(p,"format",    strings(0),     @(s) ischar_s(s) || isempty_s(s));
    addParameter(p,"keep_all_vars",  false,     @(s) islogical(s));
    addParameter(p,"verbose",   false,          @(s) islogical(s));
    addParameter(p,"history",   "",             @(s) ischar_s(s));
    addParameter(p,"excludes",  strings(0),     @(s) isempty_s(s) || ischars(s));    
    addParameter(p,"rotate",    false,          @(s) islogical(s) || (isnumeric(s) && (s==0 || s==1)));     % rotate axes for lon/lat/time or lon/lat/day-of-year variables
    addParameter(p,"chunksize", [],             @(s) isempty_s(s) || isnumeric(s));    
    addParameter(p,"unlimited", [],             @(s) isempty_s(s) || islogical(s) || (isnumeric(s) && (s==0 || s==1)));    
    addParameter(p,"clean_metadata", true,      @(s) islogical(s) || (isnumeric(s) && (s==0 || s==1)));     % clean ARRM_V2 minifile metadata.
    addParameter(p,"fillvalue", [],             @(s) isnumeric(s) && length(s)==1);
    addParameter(p,"checkonly", false,          @(s) islogical(s) || (isnumeric(s) && (s==0 || s==1)));     % check existence of files and date ranges only.
    addParameter(p,"logname",   "",             @(s) ischar_s(s));
    addParameter(p,"completeonly", false,       @(s) islogical(s) || (isnumeric(s) && (s==0 || s==1)));     % only join files where all the data required is present.
    addParameter(p,"include_zvals", false,      @(s) islogical(s) || s==1 || s==0);     % changed from true.  Now doing separate call to get zvals as separate file.
    addParameter(p,"add_source_metadata", true, @(s) islogical(s) || s==1 || s==0);
    addParameter(p,"reorder_vars", true,        @(s) islogical(s) || s==1 || s==0);
    addParameter(p,"maxmem", 4*1024*1024*1024,  @(s) isnumeric(s));                 % maxmem only used to make sure we can read in each file entirely into memory.
                           
    parse(p, varargin{:});

    latrange            = p.Results.latrange;
    lonrange            = p.Results.lonrange;
    daterange           = p.Results.daterange;
%    calendar           = p.Results.calendar;
    do_overwrite        = p.Results.overwrite;
    ncfmt               = p.Results.format;
    keep_all_vars       = p.Results.keep_all_vars;
    verbose             = p.Results.verbose;
    history             = p.Results.history;
    excludes            = string(p.Results.excludes);
    do_rotate           = p.Results.rotate;
    chunksize           = p.Results.chunksize;
    unlimited           = p.Results.unlimited;
    do_clean            = p.Results.clean_metadata;
    fillvalue           = p.Results.fillvalue;
    checkonly           = p.Results.checkonly;
    logname             = p.Results.logname;
    completeonly        = p.Results.completeonly;
    include_zvals       = p.Results.include_zvals;
    add_source_metadata = p.Results.add_source_metadata;
    reorder_vars        = p.Results.reorder_vars;
    maxmem              = p.Results.maxmem;

    nmiss_calendar = 0;
    ndiff_calendar = 0;
    ndiff_time = 0;
    
    [fdir,fbase,fext] = fileparts(outname);
    if (strlength(fdir)==0), fdir="."; end
    if (~isfolder(fdir)), error("error:  output folder %s doesn't exist or isn't writeable", fdir); end
    if (~strcmpi(fext,".nc")), error("error:  extension for output file '%s' must be '.nc", outname); end
    
    if (strlength(logname)==0)
        logname=fullfile(fdir, sprintf("%s.log", fbase));
        fid = fopen(logname,"w");
        fprintf(fid,"");
        fclose(fid);
    end
    log_print(logname, 1,  "---join_minifiles:  %s ---   %s\n", datestr(now, "yyyy-mm-dd HH:MM:SS"), vec2string(varnames));        
    
    if (isfile(outname) && ~do_overwrite && ~checkonly), log_error(logname, "Error:  output file %s already exists.  Please specify '""overwrite"",true' to overwrite", outname); end
    
    ncObj.verbose(verbose);
    
    [~,fns] = dirmatch(fn_list);        % find matches for all filenames. dirmatch does linux-style name matching.  2nd param is fully-qualified filenames. 
    
    if (isempty(fns))
        log_print(logname, 2, "no files matching: \n");
        log_print(logname, 2, "\t%s\n", fn_list);
        log_error(logname, "error...aborting");
    end
    
     nncs = length(fns);
           % check that no files are larger than maxmem [4 GB].  Will have to
            % figure out how to manage large files...
    for i=1:length(maxmem); if (maxmem(i) <=512), maxmem(i) = maxmem(i) * 1024*1024*1024; end; end
    if(length(maxmem)==2), ncObj.maxmem(maxmem(2)); end
    for i=1:nncs
        s=dir(fns(i));
        mysize = s.bytes;
        if (mysize > maxmem)
            log_error(logname, "error: file %s size too large:  %d bytes\n", basename(fns(i)), mysize);
        end
    end
    
    start_dnums = nan(nncs,1);
    end_dnums   = nan(nncs,1);     
    ncs = cell(nncs,1);
    keepers = true(nncs,1);
    first = true;
    log_print(logname, 1, "reading in ncdf schemas for %d files\n", nncs);
    for i=1:nncs
        show_progress(i, nncs);
        try
            ncs{i} = ncdf(fns(i), "create",false);
        catch me
            log_print(logname, 2, "error opening netcdf file %s\n", fns(i));
            rethrow(me);
        end
        ncs{i}.loadvars([], true);  % load dimensions only.
        [tunits, incal, svec, dsince, tmname] = ncdf_get_time_info(ncs{i});
        if (isempty(incal) || strlength(incal)==0)                              % kludge here to work with sheffield, which doesnt specify what the calendar is.
            v = ncs{i}.getvar(tmname);
            v.putatt("calendar", "standard");
            incal = "standard";
            nmiss_calendar = nmiss_calendar+1; 
        end 
        if (contains(tunits,"minutes"))                                         % also kludge to fix sheffield, which has time in "minutes since"
            v = ncs{i}.getvar(tmname);
            tunits = strrep(tunits,"minutes", "days");
            dsince = dsince/1440;
                % and update the units and days-since info in the nc data structure. 
            v.putatt("units", tunits);
            v.vdata = dsince;
        end
%       dsince = floor(dsince+1e-5);    % truncate days-since to whole days.  
            % make sure time info matches.
        if (i==1)
            time_units = tunits;
            out_calendar = incal;
            start_vec = svec;
%           days_since = dsince;
            time_name = tmname;
            [lat_name, lon_name] = ncdf_get_llt_dimnames(ncs{1});
        else
            if (calendar_length(out_calendar) ~= calendar_length(incal))
                ndiff_calendar = ndiff_calendar + 1;
                log_error(logname, "error:  (fatal): calendar  mismatch on input files:  initial:  %s ;  %s for file %s\n", time_units, tunits, basename(ncs{i}.Filename));
            elseif (~strcmp(time_units,tunits)  || datenum_cal(svec, incal) ~= datenum_cal(start_vec, out_calendar) || ~strcmp(tmname, time_name))
                ndiff_time = ndiff_time + 1;
                if (ndiff_time <= 5)
                    log_print(logname, 2,"warning:  (non-fatal): time info mismatch on input files:  initial:  %s ;  %s for file %s\n", time_units, tunits, basename(ncs{i}.Filename));
                end
            end
        end
        sdnum = datenum_cal(svec, incal);
        start_dnums(i) = sdnum + dsince(1);
        end_dnums(i)   = sdnum + dsince(end);
        
        if (~isempty(daterange))
            sdnum = datenum_cal(daterange(1,:), incal);
            ednum = datenum_cal(daterange(end,:), incal);
            if (start_dnums(i) > ednum || end_dnums(i) < sdnum)
                keepers(i) = false;
                continue
            end
        end        
            % make sure time dimension is of type double.
        ncs{i}.getvar(time_name).cast("double");
        file_lats = ncs{i}.getvardata(lat_name);
        file_lons = ncs{i}.getvardata(lon_name);
        
                % skip this file if data is not in lat/lon range.
        if ((~isempty(latrange) && file_lats(1) > latrange(2) || file_lats(end) < latrange(1)) || (~isempty(lonrange) && file_lons(1)>lonrange(2) || file_lons(end)<lonrange(1)))
            keepers(i) = false;
            continue;
        end
        
        ncs{i}.clear_unlimited_dimensions();
        
        
        if (first)
            outlats = file_lats;
            outlons = file_lons;
            first = false;
        else
            outlats = union(outlats, file_lats);
            outlons = union(outlons, file_lons);
        end        
            
    end
    if (~any(keepers)), error("error:  no files with data in lat/lon range"); end
    ncs = ncs(keepers);
    start_dnums = start_dnums(keepers);
    end_dnums = end_dnums(keepers);
    
        % limit the outlats and outlons to the specified ranges.
        
    if (isempty(latrange))
        latrange = [min(outlats), max(outlats)];
    else    
        keepers = outlats >= latrange(1) & outlats <= latrange(end);
        outlats = outlats(keepers);
    end
    
    if (isempty(lonrange))
        lonrange = [min(outlons),max(outlons)];
    end
    if (any(outlons < 0))
        lonrange = mod(lonrange+180, 360)-180;
    end
    keepers = outlons >= lonrange(1) & outlons <= lonrange(2);
    outlons = outlons(keepers);
    
    if (isempty(daterange))
        start_dvec = datevec_cal(min(start_dnums), out_calendar);
        end_dvec   = datevec_cal(max(end_dnums),   out_calendar);
        daterange = [start_dvec(1:end); end_dvec(1:end)];
    end
    
        % fix date range if 360-day day end-day is day 31 of month.
    if (calendar_length(out_calendar)==360)
        daterange(:,3) = min(daterange(:,3),30);
    end
    
        % trim down the varnames list, excluding any dimension variables and variables in the excludes list.
    if (isempty(varnames)), varnames = string(ncs{1}.varlist()); end
    
    [keeper_varnames, unneeded] = get_varlist(varnames, include_zvals, keep_all_vars, excludes, ncs, logname, verbose);

    if (ndiff_time > 5)
        log_print(logname, 2,"note:  %d files with different start_times\n", ndiff_time);
    end
            
%     if (ndiff_calendar > 0)       % different calendars is now a fatal error.
%         log_print(logname, 2,"note:  %d files with different calendar info\n", ndiff_calendar);
%     end
%             
    if (nmiss_calendar > 0)
        log_print(logname, 2, "note:  %d files with no calendar specification in metadata\n", nmiss_calendar);
        log_print(logname, 2, "       Assuming ""standard"" for calendar type\n");
    end
            
end

function [keeper_varnames, unneeded] = get_varlist(varnames, include_zvals, keep_all_variables, excludes, ncs, logname, verbose)
    % 
    in_varlist = ncs{1}.varlist();
    in_dimlist = ncs{1}.dimlist();
    
%   clim_vars    = varnames(is_climate_variable(varnames));
    
    % first, make sure all the requested varnames are in the varlist from the input files.
    
    if (keep_all_variables)
        varnames = in_varlist;
    else

        missing = setdiff(varnames, in_varlist);
        if (~isempty(missing))
            error("input nd files are missing variables:  %s", join(missing,", ")); 
        end
    end
    
        % include the zvals variable if not already in varnames.
        
                % requested variables, in the same order as they appear in the input file list.
    [~,~,ix] = intersect(varnames, in_varlist, "stable");     
    ordered_varnames = in_varlist(ix); 

        % set up output list of variables, and flags of the ones we want to keep, and flag the dimensions we want to keep.
    keepers = false(1,length(in_varlist));
    dim_keepers = false(1, length(in_dimlist));
    
    keeper_varnames = [to_row(in_varlist); to_row(in_varlist)];  
    
    for i=1:length(ordered_varnames)
        vname = ordered_varnames(i);
        if (any(strcmp(vname, excludes))), continue; end     % skip this variable if its in the excludes list.
            % flag variable as a keeper.
        ix = find(strcmp(vname, in_varlist));
        keepers(ix) = true; %#ok<FNDSB>
        
            % and include its zvals variable
        if (is_climate_variable(vname))
            if (include_zvals)
                std_zvalsname = sprintf("%s_zvals", vname); %clim_vars(i));
                    % see if the std zvals name is in the input list.
                ix = find(strcmp(std_zvalsname, keeper_varnames(2,:)),1);
                if (~isempty(ix))       % found it.  flag it as a keeper.
                    keepers(ix) = true;
                                        % but if it isn't, use the first available old-style zvals variable.
                                        % normally only 1, but in the future, there may be multiple downscaled variables and their zvals, so 
                                        % we'll assume their zvals are in the same order as the climate variables were.
                else
                    ix = find(strcmp("zvals", keeper_varnames(2,:)) | strcmp("packed_zvals", keeper_varnames(2,:)),1);
                    if (isempty(ix))
                        error("error:  cannot find a zvals variable for %s", varnames(i));  % oops.  no std zvals variable found, and no old-style ones left to use either.
                    end
                    keepers(ix) = true;
                    keeper_varnames(2,ix) = std_zvalsname;
                end
            end
                % flag the dimensions as keepers, in both the variables' keepers and in the dim_keepers arrays.
                    % get the variable's list of dimensions
            v = ncs{1}.getvar(vname);
            dims = v.dimlist;
            for j=1:length(dims)
                dimname = dims(j);
                jx = find(strcmp(dimname, keeper_varnames(2,:)),1);
                if (isempty(jx))
                    error("error:  required dimension %s not found in list of variables", dimname);
                else
                    keepers(jx) = true;
                end
                dx = find(strcmp(dimname, in_dimlist),1);
                if (isempty(dx))
                    error("error:  required dimension %s not found in list of dimensions", dimname);
                else
                    dim_keepers(dx) = true;
                end
            end                        
        end
    end
   
        % now create the list of unneeded dimensions and variables.
    if (keep_all_variables)
        unneeded.variables  = strings(0);
        unneeded.dimensions = strings(0);
    else
        keeper_varnames = keeper_varnames(:, keepers);
        unneeded.variables = in_varlist(~keepers);
        unneeded.dimensions = in_dimlist(~dim_keepers);
    end
    
        % now make sure each netcdf file has all the required dimensions and variables.
        %   We do that here so we don't waste time copying data before we discover that we're missing stuff.
        
    nncs = length(ncs);
    keeper_vars = in_varlist(keepers);
    keeper_dims = in_dimlist(dim_keepers);

    log_print(logname, 1, "checking required dimensions and variables in each file\n");

    for i=2:nncs
        difs = setdiff(keeper_vars, ncs{i}.varlist());
        if (~isempty(difs))
            error("error: missing variables (%s) in input file %s", join(difs,", "), ncs{i}.Filename);
        end
        difs = setdiff(keeper_dims, ncs{i}.dimlist());
        if (~isempty(difs))
            error("error: missing dimensions (%s) in input file %s", join(difs,", "), ncs{i}.Filename);
        end
        if (verbose), show_progress(i, nncs); end
    end
end
            
        
%     for i=1:length(clim_vars)
%     
%         if (~isempty(setdiff(varnames, inv_varlist)))
%             if (any(length(matches) < length(varnames)
%     if (include_zvals)
%             % add zvals variable to list of varnames if missing.
%         if (~any(is_zvals_variable(varnames)))
%             invarnames = nc.varlist;
%             ix = find(is_zvals_variable(invarnames));
%             if (isempty(ix))
%                 error("error:  no zvals variable found in input nc files"); 
%             elseif (length(ix) > 1)
%                 error("too many zvals variables found in input.  Please specify zvals variable explicitly in list of varnames to include");
%             end
%             varnames(end+1) = invarnames(ix);
%         elseif (sum(is_zvals_variables(varnames))>1)
%             error("too many zvals variables found in input.  Please specify zvals variable explicitly in list of varnames to include");
%         end
%         
%         zx = find(is_zvals_variable(varnames),1);
%         zvalsname = varnames(zx);
%         
%             % create the proper zvals varname, and assign it to the output variable name
%         jx = find((is_climate_variable(varnames) & ~is_zvals_variable(varnames)),1);  % find first climate variable in list of keeper varnames.
%         if (isempty(jx)), error("error:  no climate variable in list of keeper variables"); end
%         zvalsname(2) = sprintf("%s_zvals",varnames(jx));
%     else
%         zvalsname = strings(0);
%     end
%     
%     for i=1: nncs
%         vlist = intersect(varnames, string(ncs{i}.varlist));  % get list of keeper variables in current nc
%         
%         if (isempty(vlist))
%             log_error(logname, "No matching data for varnames %s in file %s", vec2string(varnames), ncs{i}.Filename);
%         end
% 
%         for j=1:length(vlist)                                 % for each keeper variable:
%             v=ncs{i}.getvar(vlist(j));
%             dlist = v.dimlist();                              % get list of it's dimension variables
% %             no_var = false;
%             for k = 1:length(dlist)                           % for each dimension variable 
% %                  try
%                     vv=ncs{i}.getvar(dlist(k));                   % get nc's dimension variable
% %                 catch
% %                     no_var = true;                            % no variable for this dimension...probably "bnds".
% %                 end
%                 if (~any(strcmp(dlist(k), outdimlist)))       % if it's not in ncout's list of dimensions yet
% %                     if (no_var)
% %                         outdimlist=cat(2, outdimlist,dlist(k));   % and add it to the list of ncout's dimension variables                                      
% %                     else
%                         if (~isempty(fillvalue) && is_climate_variable(vv.Name))
%                             vv.FillValue = cast(fillvalue, vv.Datatype);             % update the fillvalue if necessary.
%                         end
%                         if (~any(strcmp(ncout.varlist(), vv.Name)))
%                         	ncout.putvar(vv.clone());                 % put it there
%                         end
%                         outdimlist=cat(2, outdimlist,dlist(k));   % and add it to the list of ncout's dimension variables                  
% %                     end
%                 else                                          % otherwise
%                     vdimvarout = ncout.getvar(dlist(k));      % add any new values to ncout's copy of the dimension variable                       
%                     vdimvarout.vdata = union(vdimvarout.vdata, vv.vdata, "sorted");  % (default is sorted...just being obvious here.)
%                     vdimvarout.Size = length(vdimvarout.vdata);         % and update all the sizes.
%                     vdimout = ncout.getdim(dlist(k));
%                     vdimout.Length = vdimvarout.Size;
%                 end
%             end
%         end
%         show_progress(i, nncs);
%     end
%     
% end
%     


function permute_orders = rotate_time_dependent_variables(ncout)

    vlist = ncout.varlist();    
    permute_orders = cell(length(vlist),1);
    
    for i=1:length(vlist)
        v = ncout.Variables(i); 
        dlist = v.dimlist();
        if (is_lltvar(v))
            permute_orders{i} = fliplr(1:length(dlist));
            v.rotate_dims(permute_orders{i});     % reverse the dimensions.
            
        else
            permute_orders{i} = [];
        end
    end
end

function tf = is_lltvar(v)

    [latix, ~, ~, ~] = v.diminfo('lat');
    [lonix, ~, ~, ~] = v.diminfo('lon');
    [timix, ~, ~, ~] = v.diminfo('time');
    [doyix, ~, ~, ~] = v.diminfo('doy');
    
    tf = (~isempty(latix) && ~isempty(lonix) && (~isempty(timix) || ~isempty(doyix)));
end



%function [ncout, ncs, permute_orders] = netcdf_init(ncs, fns, latrange, lonrange, daterange, outlats, outlons, keeper_varnames, unneeded, outname, nc_time_info, do_rotate, chunksize, unlimited, do_clean, ncformat, new_history, fillvalue, add_source_metadata, logname)
function  [ncout, ncs, permute_orders] = netcdf_init(ncs, fns,                                outlats, outlons, keeper_varnames, unneeded, outname, nc_time_info, do_rotate, chunksize, unlimited, do_clean, ncformat, new_history, fillvalue, add_source_metadata, logname)
%   Creates output ncdf.
% scan through all ncs, get complete list of dimensions and their values.
% updates the dimension info in ncout to cover the full range of values for each dimension.
%
%  varnames         (1,:) are names of all variables to include from input files, including dimension variables
%                   (2,:) is name to use in output file.
%                   On return, varnames may have been reordered by size for faster file access.

    npts  = length(nc_time_info.days_since);
    nlats = length(outlats);
    nlons = length(outlons);


        % make the output ncdf
    ncout = ncs{1}.clone();
    ncout.Filename = outname;    
    
        % remove unneeded variables
    varlist = ncout.varlist();
    keepers = true(size(varlist));
    for i=1:length(unneeded.variables)
        ix = find(strcmp(unneeded.variables(i), varlist),1);
        keepers(ix) = false;
    end
    ncout.Variables = ncout.Variables(keepers);
    
        % remove unneeded dimensions
    dimlist = ncout.dimlist();
    keepers = true(size(dimlist));
    for i=1:length(unneeded.dimensions)
        ix = find(strcmp(unneeded.dimensions(i), dimlist),1);
        keepers(ix) = false;
    end    
    ncout.Dimensions = ncout.Dimensions(keepers);
    
        % rename any variables (usually only zvals) that need it
    for i=1:size(keeper_varnames, 2)
        varlist =  ncout.varlist();
        if (~strcmp(keeper_varnames(1,i), keeper_varnames(2,i)))
            ix = find(strcmp(keeper_varnames(1,i), varlist),1);
            ncout.Variables(ix).Name = keeper_varnames(2,i);
        end
    end

        % turn off chunking for any non-climate variables

    for i=1:length(ncout.Variables)
        v = ncout.Variables(i);
        if (~is_climate_variable(v.Name) || length(v.Size) ~= length(chunksize))
            v.ChunkSize = [];
        end
    end
        
        % add source files' metadata.  
        % ARRM_V2 output has group "Source_metadata" with (at a minimum) the files used for OBS, HIST and MODEL.
        % Newer output usually includes the actual metadata as well, flagged by presence of Attributre "brief" set to false.
        
    if (add_source_metadata)
        try
            g=ncout.getgrp("Source_metadata");
            OBS_src   = g.getattvalue("OBS_src");
            HIST_src  = g.getattvalue("HIST_src");
            MODEL_src = g.getattvalue("MODEL_src");
        catch
            try
                g=ncout.getgrp("DownscalingParams");
                OBS_src   = g.getattvalue("obsnames");
                HIST_src  = g.getattvalue("histnames");
                MODEL_src = g.getattvalue("mdlnames");
            catch me              
                log_print(logname, 2,  "error:  cannot retrieve list of input source files from %s\n", ncout.Filename);
                throw(me);
            end
        end
        if (isempty_s(HIST_src))      % some early runs had HIST_src empty if it was the first file of MODEL_src.
            HIST_src = Model_src(1); 
        end
        src_grp = make_source_group(OBS_src, HIST_src, MODEL_src, true);
        ncout.putgrp(src_grp);        
    end
   
        % update the group filenames
    update_group_filenames(ncout, outname);
    
        % make time an unlimited variable, and turn off chunking for it. 
    [~,~,~,~,~,~, timix] = ncdf_get_llt_dimnames(ncout);
            
    if (~isempty(unlimited))
        ncout.Dimensions(timix).Unlimited = unlimited;
    end
    
    
        % update chunksize to the length of the time variable if it was passed in as infinite.
    if (any(isinf(chunksize)))
            % chunksize should be in the same order as the input files.
       vnames = ncout.varlist();
       ix=find(is_climate_variable(vnames),1);  % find an output data variable's name
       v=ncs{1}.getvar(vnames(ix));             % get the variable from the input file.
       [~,~,~,~,latix_in, lonix_in] = ncdf_get_llt_dimnames(v);       % find the index positions for lat & lon.
        for i=1:length(chunksize)
            if (isinf(chunksize(i)))
                if     (i==latix_in), chunksize(i) = nlats;
                elseif (i==lonix_in), chunksize(i) = nlons;
                else ,                chunksize(i) = npts;
                end
            end
        end
    end

            % update (or clean up) up some things, like chunking, compression, and fillvalue.
    long_names = strings(0);
    for i=1:length(ncout.varlist)
        vnm = ncout.varlist(i);
        v=ncout.getvar(vnm);
        if (is_climate_variable(vnm))
            if (~is_zvals_variable(vnm))
                if (~isempty(fillvalue))
                    v.FillValue = cast(fillvalue, v.Datatype);      % fix fillvalue
                end
                    % and set missing_value and _Fillvalue if present
                    % (I think these should be removed...)
                try
                    a=v.getatt("missing_value");
                    a.Value=v.FillValue;
                catch
                end
                try
                    a=v.getatt("_FillValue");
                    a.Value = v.FillValue;
                catch
                end
                try
                    lname=v.getattvalue("long_name");
                catch
                    lname = v.Name;
                end
                long_names = join([long_names, lname], ", ");
            end
            v.ChunkSize = chunksize;
        end
        v.DeflateLevel = 0;                             % and turn off compression.
        v.Shuffle = false;

%         if (strcmp(vnm, varnames(1)))
%             
%             if (isempty(long_name) || strlength(long_name)==0)
%                 long_name=varnames(1);
%             end
%         end
    end
    
    
            % get resulting time, lat & lon dat, and update 
                    
    [latname, lonname, timename] = ncdf_get_llt_dimnames(ncout);
    londim  = ncout.getdim(lonname);
    latdim  = ncout.getdim(latname);
    timedim = ncout.getdim(timename);

            % update ncout's lats if necessary
    latvar = ncout.getvar(latname);
    latvar.vdata = outlats;
    latvar.Size = length(outlats);
    latdim.Length = length(outlats);
           
                % update ncout's lons if necessary
    lonvar = ncout.getvar(lonname);
    lonvar.vdata = outlons;
    lonvar.Size = length(outlons);
    londim.Length = length(outlons);
              
                    % next three lines currently unnecessary;  nc_time_info
                    % s/b identical to time info in ncs{1}.
            % update ncout's time info if necessary.
            % update based on output calendar
     
%   units = nc_time_info.time_units;
%   out_calendar = nc_time_info.out_calendar;

    days_since = nc_time_info.days_since;
    
    timevar = ncout.getvar(timename);
%   timevar.putatt("units",units);
    timevar.vdata = days_since;
    timevar.Size = length(days_since);
    timedim.Length = length(days_since);  
         
        % update the netcdf format to netcdf4 (or whatever was specified on input)
    if (~isempty(ncformat))
        ncout.Format = ncformat;
    elseif (strcmp(ncout.Format,'classic'))     % convert to netcdf4 from classic if not specified. 
        ncout.Format = 'netcdf4';
    end
    
        % update the history.
    if (length(fns) == 1)
        template=fns;
    else
        t1=ncs{1};
        template = basename(t1.Filename);
    end
     
    nncs = length(ncs);
%     update_history(ncout, nncs, template, keeper_varnames(2,:), out_calendar, in_calendar, new_history);
    update_history(ncout, nncs, template, keeper_varnames(2,:), new_history);
    
        % update the variables' sizes to match the region and dates we're actually keeping.
    for i=1:size(keeper_varnames,2)
        if (is_climate_variable(keeper_varnames(2,i)))
            v=ncout.getvar(keeper_varnames(2,i));
            update_var(v, londim, latdim, timedim);        
        end
    end
    
        % set up to rotate the climate variablesvariables
    if (do_rotate)
        permute_orders = rotate_time_dependent_variables(ncout);
    else
        permute_orders = cell(length(ncout.varlist()), 1);
    end
        
    if (do_clean)
        clean_up_nc_metadata(ncout, long_names, outlats, outlons);  % remove attributes, etc. that apply to a single gridbox and fix fillvalue and turn off compression stuff.
    end
        
end

function update_group_filenames(grp, outname)
    % makes sure all groups' and subgroups' filenames point to ncout.  
    grp.Filename = outname;
    for i=1:length(grp.Groups)
        update_group_filenames(grp.Groups(i), outname);
    end
end



% function update_history(ncout, nncs, inname, varnames, out_calendar, in_calendar,  new_history)
function update_history(ncout, nncs, inname, varnames, new_history)

    histstr = new_history + sprintf("%s:  joining %d minifiles with names like '%s'; variables: %s", mfilename, nncs, inname, join(varnames,", "));
    
%     if (calendar_length(out_calendar) ~= calendar_length(in_calendar))
%         histstr = histstr + sprintf("; calendar adjusted from %s to %s", in_calendar, out_calendar);
%     end
    
    ncdf_update_history(ncout, histstr, now());
end

function update_var(v, londim, latdim, timedim)
%   updates a variable's size and dimension info

    londimix = find(strcmp(londim.Name, {v.Dimensions.Name}));
    if (~isempty(londimix))
        v.Dimensions(londimix) = londim;
        v.Size(londimix) = londim.Length;
    end
    
    latdimix = find(strcmp(latdim.Name, {v.Dimensions.Name}));
    if (~isempty(latdimix))
        v.Dimensions(latdimix) = latdim;
        v.Size(latdimix) = latdim.Length;
    end
    
    timedimix = find(strcmp(timedim.Name, {v.Dimensions.Name}));
    if (~isempty(timedimix))
        v.Dimensions(timedimix) = timedim;
        v.Size(timedimix) = timedim.Length;
    end
        
end

function output_order = get_var_output_order(ncout)
%   Returns an ordering of the variables by size, with the dimension variables first (in their original order)

    dlist = ncout.dimlist();
    vlist = ncout.varlist();
    
            % put dimension variables at the head of the output list.
    dim_output_order = zeros(length(dlist),1);

            % create output order from dimension variables first.
    j=0;
    for i=1:length(dlist)
        ix = find(strcmp(dlist(i), vlist),1);
        if (~isempty(ix)) % could have a dimension which doesn't have an actual variable array, so only keep track of the ones we found in varlist.
            j=j+1;
            dim_output_order(j) = ix;   
        end
    end
    
    sizes = zeros(length(vlist),1);
    
            % get their sizes, so we can sort them by size and put the biggest variables at the end.
    for i=1:length(vlist)
        if (any(strcmp(vlist(i),dlist)))
            sizes(i)=0;
        else
            sizes(i) = ncout.Variables(i).total_nbytes();
        end
    end
    
    [~,output_order] = sort(sizes);
    
        % put the dimension variables in at the head of the list in the same order they were in in the original file.

    output_order(1:j) =dim_output_order; 
end
    
function clean_up_nc_metadata(nc,long_names, lats,lons)

    [~,fname,ext] = fileparts(nc.Filename);
    fname = strcat(fname, ext);

    nc.put("/Attributes/completion_status",0)
        % remove attributes specific to a single gridbox:
        %   nsites
        %   single_gridbox
        %
    rm_attr(nc, "nsites");
    rm_attr(nc, "single_gridbox");

    try
        cmnts=nc.getattvalue("comments");
        if (contains(cmnts,"downscaling run"))
            nc.putatt("comments",sprintf("Texas Tech Climate Science Center downscaling run, %s, 1950 - 2100,  ( %.4f, %.4f) to ( %.4f, %.4f)",long_names, lats(1),lons(1),lats(end),lons(end))); 
        end
    catch
        % no comments field in sheffield data.
    end
    
    vnames = setdiff(nc.varlist(), nc.dimlist(),"stable");
    nc.putatt("data_variables",join(vnames,"|"));

    
        % these seem to get filled in by matlab, with 9.96920996838687e+36
    [latname, lonname, timename] = ncdf_get_llt_dimnames(nc);
    try       
        v=nc.getvar(timename);
        v.put("FillValue",[]);
        v.put("ChunkSize",[]);
        v.put("DeflateLevel",0);
        v.put("Shuffle",false);
%       v.put("FillValue", []);  % ? s/b/ _FillValue
        v=nc.getvar(latname);
        v.put("FillValue",[]);
        v.put("ChunkSize",[]);
        v.put("DeflateLevel",0);
        v.put("Shuffle",false);
%       v.put("FillValue", []);  % ? s/b/ _FillValue
        v=nc.getvar(lonname);
        v.put("FillValue",[]);
        v.put("ChunkSize",[]);
        v.put("DeflateLevel",0);
        v.put("Shuffle",false);
%       v.put("FillValue", []);  % ? s/b/ _FillValue
    catch
    end

    try
        g = nc.getgrp("RunParams");
        rm_attr(g,"probs");
    catch
    end
    try
        g = nc.getgrp("DownscalingParams");
%       rm_attr(g,"llgrid_size");
        rm_attr(g,"llgrid_lbl");
%       rm_attr(g,"zval_offset");
%       rm_attr(g,"zval_scaling");
%       rm_attr(grp, "gridbox");
        g.put("Attributes/gridbox", [lats(1),lons(1),lats(end),lons(end)]);
        rm_attr(g, "logname");
        rm_attr(g, "ARRM_V2_dir");
        g.put("Attributes/outname", fname);
%         g.put("Attributes/outdir",".");
%         g.put("Attributes/logname","none");
        rm_attr(g, "outdir");
        rm_attr(g, "logname");
    catch
    end

        % other stuff that shouldn't be there:
%   rm_attr(nc,"interpolation_method");        % this was removed from the global atts.  interp_mothod is stored in the DownscalingParams group. 
%   
end

function rm_attr(ncobj, attname)
    vix = ncobj.attix(attname);
    if (~isempty(vix)), ncobj.Attributes(vix) = []; end
end

         
function [nc_time_info, gaps] = check_daterange(ncs, out_calendar, daterange, outname, logname, verbose)

%   creates time variable info for date_range and out_calendar, and checks for gaps in the files.
%   Outputs:
%       time_units      string, s/b something like "days since 1900-01-01 00:00:00"
%       day_one         datevec for 1st day of output.  s/b daterange(1,1:3)
%       days_since      time variable values for output file
%       gaps            nx2 array of start and end indexes for where gaps exist in the input data.
%                           can be used to write NAs to the file after writing all data.
%    

    [tunits, in_calendar, ~, dsince, ~] = ncdf_get_time_info(ncs{1});
%     if (contains(tunits,"minutes"))     % kludge for sheffield, which has time in "minutes since"  % fixed in init_params.
%         tunits = strrep(tunits,"minutes","days");
%         dsince = dsince/1440;
%     end
    [~, timescale, ~] = nc_parse_date_str(tunits);

        % get step size.  Should be 1 day.
    if (timescale ~= 1)
        log_error(logname, "cannot work with time units %s for file %s",  tunits, basename(nc.Filename));
    end
    difdays = diff(dsince);
    step = median(difdays);
    if (abs(step-1) > 1e-5)
        log_error(logname, "median time step is not 1 day:  %.8f for file %s", step, basename(nc.Filename));
    end 

        % make the time variable info

    dnum_one = datenum_cal(daterange(1,:), in_calendar);
    dnum_end = datenum_cal(daterange(2,:), in_calendar);

    ndays = floor(dnum_end - dnum_one + 1);
    day_counts = zeros(ndays,1);
    day_nc = zeros(ndays,1);
    days_since = (0:(ndays-1))';
    dvecs = datevec_cal(dnum_one + days_since, in_calendar);
    day_one = dvecs(1,:);
    nfiles_kept = 0;

    time_units = sprintf('days since %s', datestr_cal(dnum_one, in_calendar, 'yyyy-mm-dd HH:MM:SS'));

%     nc_time_info = struct("time_units", time_units, "day_one", day_one, "days_since", days_since, "out_calendar", out_calendar, "in_calendar", in_calendar);
    nc_time_info = struct("time_units", time_units, "day_one", day_one, "days_since", days_since, "out_calendar", out_calendar);
    
        % now find all the dates present in the input files.

    log_print(logname, 1, "finding date ranges in files\n");
    nncs = length(ncs);
    for i=1:nncs
        [dsince, out_of_range, none_in_range, ~, ~] = get_days_since(dnum_one, dnum_end, out_calendar, ncs{i}, logname);

        if (none_in_range)
            log_print(logname, 2, "warning:  no data in date range for %s\n", basename(ncs{i}.Filename));
            continue;
        end
        nfiles_kept = nfiles_kept + 1;
        if (out_of_range)
            log_print(logname, 2, "warning, %s contains dates outside daterange\n", basename(ncs{i}.Filename));
        end
        if ( verbose > 1)
            date1 = datestr_cal(dnum_one+dsince(1),   out_calendar);
            date2 = datestr_cal(dnum_one+dsince(end), out_calendar);
            log_print(logname, 1, "\t%5d of %5d:  %s - %s  %s\n", i, nncs, date1, date2, basename(ncs{i}.Filename));
        elseif (verbose)
            show_progress(i, nncs);
        end
            % flag days as present, and record which .
        for j=1:length(dsince)
            ix = dsince(j)+1;
            if (ix<=0 || ix> ndays), continue; end
            day_counts(ix) = day_counts(ix)+1;
            if (day_counts(ix) == 1)     % keep track of which file the date first appears in.
                day_nc(ix) = i;
            end
        end  
    end
    
        % now find all the gaps in the data
    day_flags = day_counts > 0;
    gaps = [];
    if (any(~day_flags))
        ix2 = 1;
        ix1 = find(~day_flags(ix2:end),1);      % first missing date
        while (~isempty(ix1))
            ix2 = find(day_flags((ix1+1):end),1);   % first date present after ix1
            if (isempty(ix2))
                ix2 = ndays+1;
            else
                ix2 = ix1+ix2;
            end
                % gaps array has indexes of missing data.  
            gaps(end+1,:)=[ix1,ix2-1]; %#ok<AGROW>        
            ix1 = find(~day_flags((ix2+1):end),1);
            if (~isempty(ix1)), ix1 = ix2 + ix1; end
        end
    
        log_print(logname, 2, "gaps in dates for %s\n", basename(outname));
        for i=1:size(gaps,1)
            date1 = datestr_cal(dnum_one+gaps(i,1)-1, out_calendar, 'yyyy-mm-dd');  % subtract 1 since days_since is zero-based.
            date2 = datestr_cal(dnum_one+gaps(i,2)-1, out_calendar, 'yyyy-mm-dd');
            
            log_print(logname, 2, "\t%s - %s gap for joined file: %s\n", date1, date2, basename(outname));
        end
    else
        log_print(logname, 1, "No date gaps for %s\n", basename(outname));
    end
    
    
    
    % print list of files and date ranges in the files:
    
            % skip this check if all files are using the full date range..

    single_date_range = all(day_counts==nfiles_kept);

    if (~single_date_range && any(day_counts>1))
        log_print(logname, 2, "Warning:  Duplicate dates found in input data.  Files used, by date ranges:\n");
        ix1 = 1;
        nc1 = day_nc(ix1);
        ix2 = 1; 
        while (~isempty(ix2))
            date1 = datestr_cal(dnum_one+ix1, out_calendar, "yyyy-mm-dd");
            ix2 = find((day_nc((ix1+1):end) ~= nc1 & day_nc((ix1+1):end) ~= 0),1);
            if (isempty(ix2))
                date2 = datestr_cal(dnum_end, out_calendar, "yyyy-mm-dd");
                nc2 = nc1;
            else
                ix2 = ix1 + ix2;
                date2 = datestr_cal(dnum_one+ix2, out_calendar, "yyyy-mm-dd");
                nc2 = day_nc(ix2);
            end

            log_print(logname, 1, "\t%s - %s %s\n", date1, date2, ncs{nc1}.Filename);

            ix1 = ix2;
            nc1 = nc2;
        end
    end
end

function [days_since, out_of_range, none_in_range, dvec1, dvec2] = get_days_since(dnum_one, dnum_end, out_calendar, nc, logname)
%   returns days_since for nc, based on dnum_one (which can be different from the datenum of the days_since value in the file's time info);
%   dnum_one and dum_end are the datenums of the output file, based on the
%   data found in all the files.
%   out_of_range is true if any dates are outside of datenums dnum_one, dnum_end
%   days_since are days since dnum1
%
%       Code currently doesn't support a different output calendar than the
%       input data, and aborts if any of the input files have a different
%       calendar than the first input file's calendar
%               Obsolete:
% %  and in 'calendar' units (calendar may be different from nc's calendar (incal)
% %  If calendar is longer than incal (365 vs 360, e.g., there will be gaps in days_since.

        [tunits, incal, svec, days_since, ~] = ncdf_get_time_info(nc);
        
        if (calendar_length(incal)~= calendar_length(out_calendar)), error("error:  calendar mismatch:  incal: %s   calendar: calendar", incal, out_calendar); end

        [~, timescale, ~] = nc_parse_date_str(tunits);
        
            % get step size.  Should be 1 day.
        if (timescale ~= 1)
            log_error(logname, "cannot work with time units %s for file %s",  tunits, basename(nc.Filename));
        end
        difdays = diff(days_since);
        step = median(difdays);
        if (abs(step-1) > 1e-5)
            log_error(logname, "median time step is not 1 day:  %.8f for file %s", step, basename(nc.Filename));
        end
        
        days_since = round(days_since - .1);    % to get rid of timestamps not on midnight (some are x.5).
        
        
        dnum1 = datenum_cal(svec, incal);       % Ian, there's a bug here when incal is not equal to output calendar!!!
        dnums = dnum1 + days_since;
        dvecs = datevec_cal(dnums, incal);
        dnums = unique(datenum_cal(dvecs, out_calendar));
        out_of_range = any(dnums < dnum_one) || any(dnums >= dnum_end+1);
        none_in_range = (dnums(end) < dnum_one || dnums(1) > dnum_end);
        
        dvec1 = dvecs(1,:);
        dvec2 = dvecs(end,:);
        
        days_since = dnums - dnum_one;
end

function log_print(logname, typ, fmt, varargin)
        % if not running parallel, does fprintf of varargins to console and to file
        % else just does fprintf to file.
    fid = fopen(logname,'a');
    fprintf(fid, fmt, varargin{:});
    fclose(fid);
    fprintf(typ, fmt, varargin{:});
end

function log_error(logname, fmt, varargin)
        % if not running parallel, does fprintf of varargins to console and to file
        % else just does fprintf to file.
    fid = fopen(logname,'a');
    fprintf(fid, fmt, varargin{:});
    fclose(fid);
    error(fmt, varargin{:});
end


    
    
