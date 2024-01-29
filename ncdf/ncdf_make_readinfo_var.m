function [start, count, ixlat, ixlon, ixtime] = ncdf_make_readinfo_var(ncvar, ilats, ilons, itimes)
%    
%   Returns start and count vectors to read.  Looks at order of dimensions for ncdf variable
%   sets up start and count for call to any of the following for retrieving data from the file:
%           data = netcdf.getVar(ncid, varid, start, count)
%           data = ncread(source, varName, start, count, stride)
%           data = ncobj.readvar(varName, start, count, stride)
%           ncobj = ncobj.readnc(varName, start, count, stride)
%
%   Data can be missing lat or lon (or both) dimensions.  If so, then start and count vectors have 2 (or 1) elements.
%   and ilats or or ilons is ignored.  For example, if data only is for one longitude, data may have 3 dimensions
%   [lon,lat,time], where lon is a singleton dimensions.  Or it may only have 2 dimensions, [lat, time].
%
%   This function does not read any info from the file itself.  ncdf_var should be an ncdf Variable object obtained from
%       my_ncdf = ncdf(fname);
%       ncvar = my_ncdf.get(varName);
%
%   or an ncinfo's variable struct from
%
%       finfo = ncinfo(fname);
%           (find index of variable of interest.

%   Inputs:
%       ncvar           netcdf variable object
%       ilats           start, ending lat index (1-based).  If singleton, then assume retrieve only 1 lat
%       ilons           start, ending lon index (1-based).  If singleton, then assume retrieve only 1 lon
%       itimes          start, ending times index.  if singleton, then assume retrieve only 1 time.
%
%   Outputs:
%       start           start vector for ncobj.getvardata (subtract 1 for netcdf.getVar or ncget_ic)
%       count           count vector for ncobj.getvardata (or netcdf.getVar or ncget_ic)
%       ixlat, ixlon,   dimension index for lat, lon & time.  (matlab 1-based dimensions)
%       ixtime          ixlat tells which dimension (1,2 or 3) is latitude, etc..
%
%                       To get data ordered [lat,lon,time], use 
%                           data = ncget_ic(fname, varName, group, start, count);
%                           mydata = permute(data, ixlat, ixlon, ixtime);
%                       NOTE: netcdf dimension order for the variables is reversed from matlab's, and zero-based.)

    latnames={'lat','latitude','LAT','LATITUDE','Latitude'};
    lonnames={'lon','longitude','LON','LONGITUDE','Longitude'};
    timenames={'time','Time','TIME'};       
    
    if (~isa(ncvar,'Variable')), throw(MEcxception('NCDF:NOTAVAR',sprintf('Error:  %s (%s):  not a ncdf Variable', ncvar.Name, class(ncvar)))); end
    
    ixlat=[];
    ixlon=[];
    ixtime=[];
    
            % find dimension indexes
    dimnames = {ncvar.Dimensions.Name};
    for i=1:length(dimnames)
        if (any(strcmpi(dimnames{i},latnames)))
            ixlat = i;
        elseif (any(strcmpi(dimnames{i},lonnames)))
            ixlon = i;
        elseif (any(strcmpi(dimnames{i},timenames)))
            ixtime = i;
        else
            throw(MException('ICSF:BAD_VARINFO', sprintf('error:  can''t get dimension info. %s', ncvar.Name)));
        end
    end
    
            % get start and count for each index.
    if (~isempty_s(ixlon))
        start(ixlon) = ilons(1);
        if (length(ilons)>1)
            count(ixlon) = ilons(end)-ilons(1)+1;
        else
            count(ixlon) = 1;
        end
    end
    
    if (~isempty_s(ixlat))
        start(ixlat) = ilats(1);
        if (length(ilats)>1)
            count(ixlat) = ilats(end)-ilats(1)+1;
        else
            count(ixlat) = 1;
        end
    end

    start(ixtime) = itimes(1);
    count(ixtime) = itimes(end)-itimes(1)+1;    % will be 1 if itimes is only 1 element...i.e., read just 1 value

    
end
