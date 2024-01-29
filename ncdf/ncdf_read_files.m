function [ncdf_out, vdata, ncdfs_in] = ncdf_read_files(fnames, varNames, latrange, lonrange, daterange, out_calendar, do_random, out_varName, is_indexes, closure_flags)
% function to read selected data from a netcdf file.  
%   data returned brackets the latrange & lonrange specified, so if a single point is specified, will generally return
%   the four surrounding gridpoints from the input data.
%       Use ncdf_read_nearest_single_location(...) to get the data for the single closest data point.
% Returns data in netcdf order (lon, lat, time)
%   which is the matlab orientation of indata(time, lat, lon).
%
%   This program joins data along the time dimension only.  Does not join different lat/lon files.
%
%   Note:   this will be REAAAALLLLY slow if file is not organized with data contiguous in time 
%           i.e., like CSC's *.llt.nc files, or for netcdf4 files, if chunked sensibly 
%
%   Assumes:  lats & lons are identical between files.  
%
%   Required Inputs:
%           fnames           filenames to read (or ncdf objects with (unmodified) meta-data of a previously opened file
%           varNames         variable to read from file
%           latrange,       [min, max] lat & lon ranges. will return all gridcells needed to cover selected area
%           lonrange            specify lonrange as 0-360, not -180 to 180. 
%                               if (lonrange crosses 360-degree boundary, give lonrange like [350, 20].
%                               can also do same for latrange, but would make sense only if you wanted to exclude the
%                               middle of the range...such as looking at arctic and antarctic regions only.
%   Optional Inputs:
%           daterange      either [yr1, yr2] or [yr1, mo1, day1; yr2, mo2, day2]  specifying start and end date range
%                               if empty, then reads all dates
%                               if yr1 is 1, starts at 1st year of data;  if yr2 is 1, ends at last year of data.
%           out_calendar    optional.  can change calendar if desired.  if empty, calendar is not changed.
%           do_random       bool.  If true, uses random-date inserts or deletes when adjusting calendar length.
%                                   Can also be an integer, in which case it is used as the seed for the random number
%                                   generator, so the results can be repeated.
%           out_varName     optional.  defaults to varName.  Specify if you need a different varname for the output ncdf
%           is_indexes      optional.  boolean.  if true, latrange, lonrange, daterange are indexes of dimensions,
%                               rather than desired ranges.
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
%
%   Returns:
%       ncdf_out        ncdf object containing appropriate metadata, including
%                           lat, lon, time, & vardata. 
%       vdata           just the variable's data
%       ncdfs_in        ncdf objects containing metadata (all meta-data plus dimension variables' data) of the input files
%                           NOTE:  if fnames is actually ncdf objects, ncdfs_in will be fnames, including all data
%                                  passed in
%
%   10/21/2021  icsf  modified to return vdata and ncdfs_in.  Also updated header comments.
%   06/07/2022  icsf  modified to return lats & lons as doubles.
%
%-------------------------------------------------
    

%  latrange, lonrange?  what if vector of lats, vector of lons?
    if (~exist('daterange',    'var')), daterange      =[]; end
    if (~exist('do_random',    'var')), do_random      =[]; end
    if (~exist('out_calendar', 'var')), out_calendar   =strings(0); end
    if (~exist('is_indexes',   'var')), is_indexes     =[]; end
    if (~exist('closure_flags','var')), closure_flags  =[1,1]; end

    if (ischar(fnames))
        ncdfs_in = {fnames};
    elseif (isstring(fnames))
        fnames = cellstr(fnames);
    else
        fnames = num2cell(fnames);
    end
    nfiles = length(fnames);
    for i=1:nfiles
        if (ischar(fnames{i}) || isstring(fnames{i}))
            ncdfs_in{i} = ncdf(fnames{i});    % read ncdf metadata
            ncdfs_in{i}.loadvars([],true);    % load all dimension variables
        else
            ncdfs_in{i} = fnames{i};
        end
    end
    if (ischar(varNames)), varNames={varNames}; end
    nfiles = length(ncdfs_in);
    if (nfiles > length(varNames))
        varNames{end+1:nfiles}=varNames{end};
    end
    if (~exist('out_varName','var') || isempty(out_varName)), out_varName =varNames{1}; end
    
    if (isempty(out_calendar) || strlength(out_calendar)==0)
        [~,out_calendar] = ncdf_get_time_info(ncdfs_in{1});
    end
        
    if (isempty(daterange))
        for i=1:length(ncdfs_in)
            [~,calendar,~,~,~,tstamps] = ncdf_get_time_info(ncdfs_in{i});
            if (i==1)
                tstart = tstamps(1);
                tend   = tstamps(end);
            else
                tstart = min(tstart, tstamps(1));
                tend   = max(tend,  tstamps(end));
            end
        end
        daterange = [datevec_cal(tstart, calendar); datevec_cal(tend, calendar)];
    elseif (length(daterange) == 2)        % years given, only.
        daterange = [daterange(1),1,1; daterange(2),12,31];
    end
    

            % get lats & lons from 1st file.
    nc1 = ncdfs_in{1};
%     lons=double(ncdf_getvar(nc1, 'lon',true));
%     lats=double(ncdf_getvar(nc1, 'lat',true));

    [lats, lons] = ncdf_get_latlons(nc1);  

    ncvar = nc1.get(varNames{1});
    FillValue = ncvar.FillValue;
    if (isempty(FillValue))
        try
            FillValue = ncvar.getattvalue('_FillValue');
            ncvar.FillValue = FillValue;
        catch
        end
    end
    try
        units = ncvar.getattvalue('units');
    catch
        units = '';
    end
    
            % get lat & lon indexes of range of points needed  (actually, the actual lats & lons, ian)
            
%   [latix,lonix, meridian_flag, qlats, qlons] = latlon_region(latrange, lonrange, lats,lons);
    [~,    ~,     ~,             qlats, qlons] = latlon_region(latrange, lonrange, lats,lons, closure_flags);
    
    vdata = [];
    ncdf_out = make_ncdf(qlats,qlons, daterange, out_calendar, out_varName, out_varName, units, FillValue, vdata);
    
    for i=1:nfiles
        ncdf_in = ncdf_read_file(ncdfs_in{i}, varNames{i}, latrange, lonrange, daterange, out_calendar, do_random, varNames{i}, is_indexes, closure_flags);
        if (isempty(ncdf_in)), continue; end
        insert(ncdf_out, out_varName, ncdf_in, varNames{i});
    end
    
    if (nargout > 1)
        v = ncdf_out.getvar(out_varName);
        vdata = v.vdata;
    end
    
end

function ncdf_out = make_ncdf(lats, lons, daterange, out_calendar, varName, longname, units, FillValue, vdata)
    %  create an ncdf object using input info.

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
    attributes = {Attribute('longname',longname),Attribute('units',units)};

    if (exist('vdata','var'))
        add_var(ncdf_out, varName, attributes, dimensions, FillValue, vdata);
    else
        add_var(ncdf_out, varName, attributes, dimensions, FillValue);
    end
    
end

function timedim = add_time(ncdf_out, daterange, out_calendar, start_date)
    % add time dimension to ncdf_out.
    % returns the dimension object as well.
    % note:  ncdf_out is a handle object, so is effectively passed by reference.  Modifying it changes calling code's
    % ncdf_out object.

    units = sprintf('days since %04d-%02d-%02d 00:00:00', start_date);
    cal_len = calendar_length(out_calendar);
    if (cal_len == 360)
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
    ncdf_out.putvar('time', time_vals, 'units',units,'calendar',out_calendar,'Dimensions',{timedim});
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


function insert(ncdf_out, out_varName, nc_in, varName )
% inserts all non-NA's from nc_in's varName into ncdf_out's varName.
%   assumes nc_in and nc_out's varName have same time, lats, lons and data order.
%   creates output var's vdata if initially empty.

    in_vdata  = nc_in.getvardata(varName);
    ncvar = ncdf_out.get(out_varName);
    if (isempty(ncvar.vdata)), ncvar.vdata = nan(size(in_vdata)); end
            % insert non-NAs into output data array.
    flags = ~isnan(in_vdata);
    ncvar.vdata(flags) = in_vdata(flags);
end
