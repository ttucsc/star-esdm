function [ ncdf_out, vdata, nc_in, varUnits, latout, lonout ] = ncdf_read_file(fname, varName, latrange, lonrange, daterange, out_calendar, do_random, out_varName, is_indexes, closure_flags)
% function to read selected data from a netcdf file.
% Returns data in netcdf order (lon, lat, time)
%   which is the matlab orientation of indata(time, lat, lon).
%
%   Inputs:
%           fname           filename to read (or ncdf object with (unmodified) meta-data of a previously opened file
%           varName         variable to read from file
%           latrange,       [min, max] lat & lon ranges. will return all gridcells needed to cover selected area
%           lonrange            specify lonrange as 0-360, not -180 to 180. 
%                               if (lonrange crosses 360-degree boundary, give lonrange like [350, 20].
%                               can also do same for latrange, but would make sense only if you wanted to exclude the
%                               middle of the range...such as looking at arctic and antarctic regions only.
%           daterange      either [yr1, yr2] or [yr1, mo1, day1; yr2, mo2, day2]  specifying start and end date range
%                               if empty, then reads all dates
%                               if yr1 is 1, starts at 1st year of data;  if yr2 is 1, ends at last year of data.
%           do_random       bool.  If true, uses random-date inserts or deletes when adjusting calendar length.
%                                   Can also be an integer, in which case it is used as the seed for the random number
%                                   generator, so the results can be repeated.
%           out_varName     optional.  defaults to varName.  Specify if you need a different varname for the output ncdf
%           out_calendar    optional.  can change calendar if desired.  if empty, calendar is not changed.
%           is_indexes      optional.  boolean.  NOTE:  if true, latrange, lonrange, daterange are indexes of dimensions,
%                               rather than desired ranges.  I.e. latrange,lonrange,daterange are n:m, not [n,m]
%           do_loaddims     optional.  boolean.  If false, assumes fname is a ncdf object with dimensions already loaded
%           closure_flags   optional.  specify closure for lat, lon ranges;  1- or 2-element vector or 2-row matrix. values can be -2, -1, 0, 1, 2
%                               if # cols == 1, applies to left & right bounds.  if length 2, left & right (or top and
%                               if # rows == 1, applies to both lat & lon.  if 2 rows, 1st is for lat, 2nd for lon.
%                               bottom) bounds.  default:  [1,1]:  "at least"
%                                   -2:     absolutely inside  (end point is inside range)
%                                   -1:     at most  (no points outside range
%                                    0:     must match exactly
%                                    1:     at least (smallest range including end point)
%                                    2:     absolutely outside (end point is outside range)

%
%   Returns:
%       ncdf_out            ncdf object containing original schema, but with lats, lons & dates restricted to specified
%                               lats, lons & date ranges.
%                               variable data is returned in matlab ordering (time, lat, lon), with time varying fastest
%                                   i.e., timesteps for a single (lat,lon) location are contiguous in the vdata.
%                                   to reorder to alternate ordering, use data = permute(...)  on variable's vdata.
%                                   NOTE:  netcdf considers this ordering (lon,lat,time)
%       vdata               just the variable data (also contained in ncdf_out
%       ncdf_in             ncdf object with all metadata and dimension data of input file (does not contain the
%                               requested variable's data)

%  icsf 10/12/2017
%   icsf 11/27/2019:  needs to be modified to read in sections so it crosses the meridian if necessary.
%   icsf 10/20/2021:  fixed daterange when ix_indexes is true.  Also added option to return vdata.
%   icsf 3/2/2022:    added varUnits to output.
%   icsf 6/7/2022:    modified to return lats & lons as doubles.
%-------------------------------------------------
    
    if (~exist('daterange','var') || isempty_s(daterange)), daterange = []; end
    if (~exist('out_calendar','var')), out_calendar = []; end
    if (~exist('do_random','var')), do_random = false; end
    if (~exist('out_varName','var') || isempty(out_varName)), out_varName = varName; end
    if (~exist('is_indexes','var') || isempty_s(is_indexes)), is_indexes = false; end
    if (~exist('closure_flags','var') || isempty(closure_flags)), closure_flags = [1,1]; end      % "at least", for both lat & lon.
    
    if (isa(fname,'ncdf'))
        nc_in=fname;
    else
        nc_in = ncdf(fname);
        nc_in.loadvars([], true);
    end
    
%     [latname, lonname, timename] = ncdf_get_llt_dimnames(nc_in);
% 
%     lats = double(ncdf_getvar(nc_in, latname, true));
%     lons = double(ncdf_getvar(nc_in, lonname, true));    
    
    [lats, lons] = ncdf_get_latlons(nc_in);  
    
    if (~isempty(daterange))
        daterange = 1*daterange;
        if (length(daterange) == 2)
            daterange = [daterange(1),1,1; daterange(2),12,31];
        end
                %   get start, end dates from ncobj if either daterange's yr is 1.
                %   this allows using true to specify start at beginning or end at end of data.
        if (daterange(1,1) == 1 || daterange(2,1) == 1)
            [~     , calendar, start_vec, days_since] = ncdf_get_time_info(nc_in, varName);   % don't need timeunits because we have start_vec.
            dnums = days_since + datenum_cal(start_vec, calendar);
            svec = datevec_cal(dnums(1),calendar);
            evec = datevec_cal(dnums(2),calendar);
            if (daterange(1,1) == 1), daterange(1,1) = svec(1); end
            if (daterange(2,1) == 1), daterange(2,1) = evec(1); end
        end
    else    
        timename = ncdf_get_timename(nc_in);
        [~,calendar,start_vec,days_since] = ncdf_get_time_info(nc_in, timename);        
        if (is_indexes)
            date_ix=1:length(days_since);
        end
        start_dnum = datenum_cal(start_vec, calendar);
        dnums = floor(start_dnum + days_since + .00001);        % truncate to start of day.  Adding .00001, which is slightly < 1 second, in case of computer arithmetic.
        daterange = datevec_cal([dnums(1); dnums(end)], calendar);
    end
    
    ncvar = nc_in.get(varName);
    try
        varUnits = ncvar.getattvalue('units');
    catch
        varUnits = "unknown";
    end
    FillValue = ncvar.FillValue;
    if (isempty(FillValue))
        try
            FillValue = ncvar.getattvalue('_FillValue');
            ncvar.FillValue = FillValue;
        catch
        end
    end
    
                % get lat, lon & time indexes & values bracketing latrange and lonrange and daterange.
                
    if (is_indexes)
        if (isempty(latrange)), latrange=1:length(lats); end
        if (isempty(lonrange)), lonrange=1:length(lons); end
        latix = latrange;
        lonix = lonrange;
            % data must be in chronological order, with no duplicates.
        keepers_ix = date_ix;
        ix1 = min(keepers_ix);
        ix2 = max(keepers_ix);
        outnpts = ix2 - ix1 + 1;
        keepers = false(1,outnpts);
        keepers(keepers_ix - ix1 + 1) = true;   
        dest_ix = 1:outnpts;
        
        latout = lats(latix);
        lonout = lons(lonix);
        
        
    else
        if (isempty(latrange)), latrange=[lats(1),lats(end)]; end
        if (isempty(lonrange)), lonrange=[lons(1),lons(end)]; end
        [latix, lonix, ~            , latout, lonout] = latlon_region(latrange, lonrange, lats,lons, closure_flags);      % closure flags defaults [1,1] gives us "at least"... range is potentially closed, but covers entire latrange or lonrange
        
%       [tunits, calendar, start_vec, timevals] = ncdf_get_time_info(ncobj, varName);
        if (~exist('calendar','var'))
            [~     , calendar, start_vec, days_since] = ncdf_get_time_info(nc_in, varName);   % don't need timeunits because we have start_vec.
        end
        
        days_since = floor(days_since);     % in case days_since are midday (.5)                                                                                           

        if (isempty_s(daterange))
            outnpts = length(days_since);
            keepers = true(1,outnpts);
        else
%           [keepers, keepers_ix, dest_ix, outnpts] = find_keepers_in_date_range(daterange, calendar, start_vec, days_since);
            [keepers, ~,          dest_ix, outnpts] = find_keepers_in_date_range(daterange, calendar, start_vec, days_since);
        end
    end
    if (sum(keepers)==0)
        ncdf_out = [];
        return;
    end
    
    if (isempty_s(out_calendar)), out_calendar = calendar; end
            % output ncdf
    ncdf_out = make_ncdf(latout, lonout, daterange, calendar, out_varName, out_varName, varUnits, FillValue, []);
            % create output array
            
            % extract data from the file and replace FillValues with nan's.
    vdata = read_data(nc_in, varName, latix, lonix, keepers, outnpts, dest_ix);
    
        % this should happen automatically.
%     if (~isempty(FillValue))
%         try
%             myfills = vdata==FillValue;
%         catch
%             fprintf("oops:  myfills\n");
%         end
%         vdata(myfills) = nan;
%     end
    
    ncvar = ncdf_out.get(varName);
    ncvar.put('vdata',vdata);
    
    if (calendar_length(out_calendar) ~= calendar_length(calendar)) 
         ncdf_adjust_calendar(ncdf_out, out_calendar, varName, true, start_vec, do_random);
    end
    
    if (nargout > 1)
        vdata = ncvar.vdata;
    end
    
end

function ncdf_out = make_ncdf(lats, lons, daterange, out_calendar, varName, longvarName, varUnits, FillValue, vdata)

    ncdf_out = ncdf();
    nlats = length(lats);
    nlons = length(lons);
    latdim = Dimension('lat',nlats, false);
    ncdf_out.putdim(latdim);
    ncdf_out.putvar('lat',lats,'units','degrees_north','long_name','latitude','Dimensions',{latdim});
    londim = Dimension('lon',nlons, false);
    ncdf_out.putdim(londim);
    ncdf_out.putvar('lon',lons,'units','degrees_east','long_name','longitude','Dimensions',{londim});
    
    timedim = add_time(ncdf_out, daterange, out_calendar, daterange(1,:));
    
    dimensions = {timedim, latdim, londim};
    attributes = {Attribute('longname',longvarName),Attribute('units',varUnits)};

    if (exist('vdata','var'))
        add_var(ncdf_out, varName, attributes, dimensions, FillValue, vdata);
    else
        add_var(ncdf_out, varName, attributes, dimensions, FillValue);
    end
    
end

function timedim = add_time(ncdf_out, daterange, out_calendar, start_date)

    calunits = sprintf('days since %04d-%02d-%02d 00:00:00', start_date);
    cal_len = calendar_length(out_calendar);
    if (cal_len == 360)
            % in case day-of-month is 31...
        daterange(1,3) = min(30, daterange(1,3));        
        daterange(2,3) = min(30, daterange(2,3));
        npts = datenum360(daterange(2,:)) - datenum360(daterange(1,:)) + 1;
        time_vals =  (0:(npts-1)) +  datenum360(daterange(1,:)) - datenum360(start_date);
    elseif (cal_len == 365)
        npts = datenum365(daterange(2,:)) - datenum365(daterange(1,:)) + 1;
        time_vals =  (0:(npts-1)) +  datenum365(daterange(1,:)) - datenum365(start_date);
    else
        npts = datenum(daterange(2,:)) - datenum(daterange(1,:)) + 1;
        time_vals =  (0:(npts-1)) +  datenum(daterange(1,:)) - datenum(start_date);
    end
    
    timedim = Dimension('time',npts,false);
    ncdf_out.putdim(timedim);
    ncdf_out.putvar('time', time_vals, 'units',calunits,'calendar',out_calendar,'Dimensions',{timedim});
end

function add_var(ncdf_out, varName, attributes, dimensions, FillValue, vdata)

    if (strcmp(ncdf_out.Format,'64bit'))
        FillValue_lbl = '_FillValue';
    else
        FillValue_lbl = 'FillValue';
    end
    if (exist('vdata','var'))
        if (isempty_s(vdata)), vdata = make_vdata(dimensions); end
        myvar = Variable(varName, vdata, 'Attributes',attributes, 'Dimensions',dimensions,FillValue_lbl,FillValue) ;
    else       
        myvar = Variable(varName,        'Attributes',attributes, 'Dimensions',dimensions,FillValue_lbl,FillValue, 'Datatype','single') ;
    end
    ncdf_out.putvar(myvar);
end


function vdata = make_vdata(dimensions)
    ndims = length(dimensions);
    siz = zeros(1,ndims);
    for i=1:ndims
        siz(i) = dimensions{i}.Length;
    end
    vdata = single(nan(siz));
end

function    vdata = read_data(ncobj, varName, latix, lonix, keepers, outnpts, dest_ix)
%       note:  FillValues not replaced upon retur.

    nlats = length(latix);
    nlons = length(lonix);
        % allocate an output array of all nans
    vdata = nan(outnpts, nlats, nlons);
        
    ncvar = ncobj.get(varName);
    ix1 = find(keepers,1);
    ix2 = find(keepers,1,'last');
%     timix = ix1:ix2;   check this, Ian!
    
    if (iscontiguous(latix) && iscontiguous(lonix))
                % lat/lon region is contiguous in file.  We can do a single read
        [start, count, ixlat, ixlon, ixtime] = ncdf_make_readinfo_var(ncvar, latix, lonix, ix1:ix2);
        if (~isempty(ncvar.vdata))
            indata = extract(ncvar, start, count);
        else
            indata = ncobj.readvar(varName, start, count, [], ncObj.verbose);
        end
        if (ixtime ~= 1 || ixlat ~= 2 || ixlon ~= 3)
            indata = permute(indata, [ixtime, ixlat, ixlon]);
        end
    
            % insert the valid data into the output array.
            % we index with keepers insteadof just assigning all of indata because the keepers list might not be continguous (there may be gaps)
        vdata(dest_ix,:,:) = indata(keepers(ix1:ix2),:,:);
        
        
    else
            % lat/long region is not contiguous.  Read each gridcell separately.
            in_npts = ix2-ix1+1;
        for lo=1:length(lonix)
            for la = 1:length(latix)
                        % figure out where and how much to read from 
                [start, count, ixlat, ixlon, ixtime] = ncdf_make_readinfo_var(ncvar, latix(la), lonix(lo), ix1:ix2);
                indata = ncobj.readvar(varName, start, count,[], ncObj.verbose);
                indata = reshape(indata, in_npts,1);
                if (ixtime ~= 1 || ixlat ~= 2 || ixlon ~= 3)
                    indata = permute(indata, [ixtime, ixlat, ixlon]);
                end

                    % insert the valid data into the output array.
                    % we index with keepers insteadof just assigning all of indata because the keepers list might not be continguous (there may be gaps)
                vdata(dest_ix,la,lo) = indata(keepers(ix1:ix2));
            end
        end
    end
                
end

    % extracts, equivalent to an ncread, with a start and count.  
    % there's probably a better way to do this...
function indata = extract(ncvar, start, count)
    indata = ncvar.vdata;
    if (all(start==1) && all(count == ncvar.Size)), return; end     
    ndims = length(start);
    iend   = start + count - 1;
    if (ndims == 2)        
        indata = indata(start(1):iend(1), start(2):iend(2));
    elseif (ndims == 3)
        indata = indata(start(1):iend(1), start(2):iend(2), start(3):iend(3));
    elseif (ndims == 4)
        indata = indata(start(1):iend(1), start(2):iend(2), start(3):iend(3), start(4):iend(4));
    elseif (ndims == 5)
        indata = indata(start(1):iend(1), start(2):iend(2), start(3):iend(3), start(4):iend(4), start(5):iend(5));
    elseif (ndims == 6)
        indata = indata(start(1):iend(1), start(2):iend(2), start(3):iend(3), start(4):iend(4), start(5):iend(5), start(6):iend(6));
    elseif (ndims == 7)
        indata = indata(start(1):iend(1), start(2):iend(2), start(3):iend(3), start(4):iend(4), start(5):iend(5), start(6):iend(6), start(7):iend(7));
    else
        error("extract:  can't handle %d dimensions", ndims);
    end
end

% function    vdata = read_data(ncobj, varName, latix, lonix, timix, prenans, postnans)
%     ncvar = ncobj.get(varName);
%     [start, count, ixlat, ixlon, ixtime] = ncdf_make_readinfo_var(ncvar, latix, lonix, timix);
%     vdata = ncobj.readvar(varName, start, count);
%     if (ixtime ~= 1 || ixlat ~= 2 || ixlon ~= 3)
%         vdata = permute(vdata, [ixtime, ixlat, ixlon]);
%     end
%     
%     [~, nlats, nlons] = size(vdata);
%     
%     vdata = [nan(prenans, nlats, nlons); vdata; nan(postnans, nlats, nlons)];
% end

