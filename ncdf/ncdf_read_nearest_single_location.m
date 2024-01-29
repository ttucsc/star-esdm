function [ vdata, latout, lonout, out_datevecs, out_calendar, units ] = ncdf_read_nearest_single_location(fname, varName, lat, lon, daterange, out_calendar, do_random)
% [ vdata, latout, lonout, out_datevecs, out_calendar ] = ncdf_read_nearest_single_location(fname, varName, lat, lon, daterange, out_calendar, do_random)
%
% function to read specified data from a netcdf file for nearest location to the lat,lon.
% Returns data in netcdf order (lon, lat, time)
%   which is the matlab orientation of indata(time, lat, lon).
%
%   Inputs:
%           fname           filename to read (or ncdf object with (unmodified) meta-data of a previously opened file
%           varName         variable to read from file
%           lat,            lat & lon of desired location. will return closest gridcell to lat/lon point.  
%                               NOTE:  NOT INTERPOLATED TO SPECIFIED POINT
%           lon                 specify lonrange as 0-360, not -180 to 180. 
%           daterange      either [yr1, yr2] or [yr1, mo1, day1; yr2, mo2, day2]  specifying start and end date range
%                               if empty, then reads all dates
%                               if yr1 is 1, starts at 1st year of data;  if yr2 is 1, ends at last year of data.
%           do_random       bool.  If true, uses random-date inserts or deletes when adjusting calendar length.
%                                   Can also be an integer, in which case it is used as the seed for the random number
%                                   generator, so the results can be repeated.
%           out_calendar    optional.  can change calendar if desired.  if empty, calendar is not changed.
%           is_indexes      optional.  boolean.  if true, latrange, lonrange, daterange are indexes of dimensions,
%                               rather than desired ranges.
%
%   Returns:
%       vdata               1-D column array covering daterange (even if available data doesn't cover entire period)
%                               variable data is returned in matlab ordering (time, lat, lon), with time varying fastest
%                               DFata matching FillValue of file is replaced with NAs
%       latout, lonout      Actual lat & lon of point returned
%       out_datevecs        datevecs of output dates.

%   icsf 12/20/2019

%-------------------------------------------------
    
    if (~exist('daterange','var')), daterange = []; end
    if (~exist('out_calendar','var')), out_calendar = []; end
    if (~exist('do_random','var')), do_random = false; end
    
    if (isa(fname,'ncdf'))
        ncobj=fname;
    else
        ncobj = ncdf(fname);
        ncobj.loadvars([], true);
    end
    [latname, lonname, timename] = ncdf_get_llt_dimnames(ncobj);

    lats = ncdf_getvar(ncobj, latname, true);
    lons = ncdf_getvar(ncobj, lonname, true);
    
    units = ncdf_getvar_info(ncobj, varName);
        
    if (~isempty(daterange))
        daterange = double(daterange);  % in case it was an integer or logical type...
        if (length(daterange) == 2)
            daterange = [daterange(1),1,1; daterange(2),12,31];      % ok if we don't go to end of the day (23:59:59) because we're truncating all hours to 00:00:00 below.
        end
                %   get start, end dates from ncobj if either daterange's yr is 1.
                %   this allows using true to specify start at beginning or end at end of data.
        if (daterange(1,1) == 1 || daterange(2,1) == 1)
            [~     , calendar, start_vec, days_since] = ncdf_get_time_info(ncobj, varName);   % don't need timeunits because we have start_vec.
            dnums = days_since + datenum_cal(start_vec, calendar);
            svec = datevec_cal(dnums(1),calendar);
            evec = datevec_cal(dnums(2),calendar);
            if (daterange(1,1) == 1), daterange(1,1) = svec(1); end
            if (daterange(2,1) == 1), daterange(2,1) = evec(1); end
        end
    else            
        [~,cal,startvec,days_since] = ncdf_get_time_info(ncobj, timename);        
        start_dnum = datenum_cal(startvec, cal);
        dnums = start_dnum + floor(days_since);
        daterange = floor(datevec_cal([dnums(1); dnums(end)], cal));
        daterange = daterange(:,1:3);    
    end
    
    ncvar = ncobj.get(varName);
    FillValue = ncvar.FillValue;
    if (isempty(FillValue))
        try
            FillValue = ncvar.getattvalue('_FillValue');
            ncvar.FillValue = FillValue;
        catch
        end
    end
    
            % get lat, lon & time indexes & values bracketing latrange and lonrange and daterange.

%   [latix, lonix, ~            , latout, lonout] = latlon_region(latrange, lonrange, lats,lons, closure_flags);      % closure flags defaults [1,1] gives us "at least"... range is potentially closed, but covers entire latrange or lonrange
    [latix, lonix, latout, lonout] = closest(lat, lon, lats, lons, 1, 1);      % closure flags defaults [1,1] gives us "at least"... range is potentially closed, but covers entire latrange or lonrange

%       [tunits, calendar, start_vec, timevals] = ncdf_get_time_info(ncobj, varName);
    if (~exist('calendar','var'))
        [~     , calendar, start_vec, days_since] = ncdf_get_time_info(ncobj, varName);   % don't need timeunits because we have start_vec.
    end

    days_since = floor(days_since);     % in case days_since are midday (.5)                                                                                           

    if (isempty_s(daterange))        
        dnums = datenum_cal(start_vec, calendar) + days_since;
        daterange = datevec_cal([dnums(1);dnums(end)], calendar);
        daterange = daterange(:,1:3); 
        outnpts = length(days_since);
        keepers = true(outnpts,1);
        dest_ix = keepers;
    else
%       [keepers, keepers_ix, dest_ix, outnpts] = find_keepers_in_date_range(daterange, calendar, start_vec, days_since);
        [keepers, ~,          dest_ix, outnpts] = find_keepers_in_date_range(daterange, calendar, start_vec, days_since);
    end
    
    if (sum(keepers)==0)
        vdata = [];
        out_datevecs = [];
        return;
    end
    
    if (isempty_s(out_calendar)), out_calendar = calendar; end
            % output ncdf
            % create output array
            
            % extract data from the file and replace FillValues with nan's.
    vdata = read_data(ncobj, varName, latix, lonix, keepers, outnpts, dest_ix);
    if (~isempty(FillValue))
        myfills = vdata==FillValue;
        vdata(myfills) = nan;
    end
        
    if (calendar_length(out_calendar) ~= calendar_length(calendar))
        my_days_since = 0:(outnpts-1);
        [ ~, out_days_since, vdata ] = ncadjust_calendar( daterange(1,:), my_days_since, vdata, calendar, out_calendar, do_random, daterange(1,:), 1  );
    else
        out_days_since = 0:(outnpts-1);
    end

    if (nargout > 3)
        if (isempty(daterange))
            [~,cal,startvec,days_since] = ncdf_get_time_info(ncobj, timename);        
            start_dnum = datenum_cal(startvec, cal);
            dnums = start_dnum + floor(days_since);
            daterange = floor(datevec_cal([dnums(1); dnums(end)], cal));
        end            
        dnums = datenum_cal(daterange(1,:),out_calendar) + out_days_since;
        out_datevecs = datevec_cal(dnums, out_calendar);
    end
    
end


function    vdata = read_data(ncobj, varName, latix, lonix, keepers, outnpts, dest_ix)
%       note:  FillValues not replaced upon return.

    nlats = length(latix);
    nlons = length(lonix);
        % allocate an output array of all nans
    vdata = nan(outnpts, nlats, nlons);
        
    ncvar = ncobj.get(varName);
    ix1 = find(keepers,1);
    ix2 = find(keepers,1,'last');
    
    [start, count, ixlat, ixlon, ixtime] = ncdf_make_readinfo_var(ncvar, latix, lonix, ix1:ix2);
    indata = ncobj.readvar(varName, start, count);
    if (ixtime ~= 1 || ixlat ~= 2 || ixlon ~= 3)
        indata = permute(indata, [ixtime, ixlat, ixlon]);
    end

        % insert the valid data into the output array.
    vdata(dest_ix,:,:) = indata(keepers(ix1:ix2),:,:);    
                
end
