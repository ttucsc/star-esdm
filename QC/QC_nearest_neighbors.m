function [neighbors, siteTbl] = QC_nearest_neighbors(stn, siteTbl, nnbrs, maxdist, minpct, maxdegrees, start_date, end_date, show_dates, csvname)
%function [neighbors, siteTbl] = QC_nearest_neighbors(stn, siteTbl, nnbrs, maxdist, minpct, maxdegrees, start_date, end_date, show_dates, csvname)
% locate nearest neighbors to a station or lat/lon
%
%   stn         station.  Can be:
%                   ix          station index in siteTbl
%                   stnID       station ID (char or string)
%                   [lat,lon]   lat, lon of point
%   siteTbl     siteTbl from netcdf file.  can be:
%                   ncName      name of netcdf file
%                   siteTbl     table extracted via cal to QC_get_site_tbl
%   nnbrs       # of neighbors to retrieve (max...may retrieve fewer)
%                   if 0, empty or missing, retrieves all within maxdist
%   maxdist     max dist to neighbors, in km
%                   if 0, empty or missing, maxdist is set by maxdegrees
%   minpct      min percentage of valid points
%                   excludes sites with pct_valid < minpct
%                   if 0, empty or missing, does not limit by pct_valid
%   maxdegrees  max degrees in lat/lon to start search.  defaults to lat=5 degrees, lon = 5/cos(lat)
%

    if (~istable(siteTbl)), siteTbl = QC_get_site_table(siteTbl); end
    if (ischar_s(stn))
        stn = QC_get_site_table(siteTbl, "stnID", stn);
        latpt = stn.lat(1);
        lonpt = stn.lon(1);
    elseif (istable(stn))
        latpt = stn.lat(1);
        lonpt = stn.lon(1);
    elseif (isnumeric(stn) && length(stn) == 2)
        latpt = stn(1);
        lonpt = stn(2);
    elseif (isnumeric(stn) && length(stn) == 1)
        stn = siteTbl(stn,:);
        latpt = stn.lat(1);
        lonpt = stn.lon(1);
    else
        error('QC_nearest_neighbors:  bad stn info');
    end
    
    if (~exist('show_dates','var') || isempty(show_dates))
        show_dates = (nargout == 0);
    end
    
    neighbors = siteTbl;
    neighbors.index = (1:size(neighbors,1))';

    if (nargin == 2)
        nnbrs = 1;
    else
        if (~exist('nnbrs','var') || isempty_s(nnbrs) || nnbrs == 0),                 nnbrs = 0;          end
        if (~exist('maxdist','var') || isempty_s(maxdist)),                           maxdist = 0;        end
        if (~exist('minpct','var') || isempty_s(minpct)),                             minpct = 0;         end
        if (~exist('maxdegrees','var') || isempty_s(maxdegrees) || maxdegrees == 0),  maxdegrees = 5.0;   end
    end
    
    maxlatdegrees = maxdegrees;
    maxlondegrees = maxdegrees / cos(latpt*pi/180);
    
    if (lonpt > 180), lonpt = -360 + lonpt; end  % make sure lonpt is in range -180 to 180

    starters = (abs(neighbors.lat-latpt) < maxlatdegrees) & (abs(neighbors.lon - lonpt) < maxlondegrees);
    
    neighbors = neighbors(starters,:);
        
    if (~isempty_s(minpct) && minpct > 0)
        valid = (neighbors.pctValid >= minpct);
        neighbors = neighbors(valid,:);
    end
    nstations = size(neighbors,1);

    latpt = repmat(latpt,nstations,1);
    lonpt = repmat(lonpt,nstations,1);
    
    
    scaling = 2*pi*6371/360;
    [dist, az] = distance(latpt, lonpt, neighbors.lat, neighbors.lon); 
    neighbors.dist = dist * scaling;
    neighbors.az = az;
    neighbors = sortrows(neighbors,'dist');
    
    if (~isempty_s(maxdist) && maxdist > 0 && maxdist < neighbors.dist(end))
        lastix = find (neighbors.dist > maxdist, 1) - 1;
    else
        lastix = nstations;
    end
    
    if (~isempty_s(nnbrs) && nnbrs > 0)
        lastix = min(lastix, nnbrs);
    end
    
    neighbors = neighbors(1:lastix, :);
    
    if (exist('csvname','var') && ~isempty(csvname))
        t = neighbors;
        if (~ischar_s(t.startDate))
            try
                z = datestr(t.startDate,'yyyy-mm-dd');
                t.startDate = z;
                z = datestr(t.endDate,'yyyy-mm-dd');
                t.endDate =z;
            catch
            end
        end
        writetable(t,csvname,'QuoteStrings',true);
    end

    
    if (exist('start_date','var') && ~isempty_s(start_date))
        if (islogical(start_date) && ~start_date)
            return;
        else
            if (~exist('end_date','var'))
                neighbors = QC_get_data(neighbors, [], start_date);         % might be able to replace [] with nc if we read it earlier anywhere, Ian...
            else
                neighbors = QC_get_data(neighbors, [], start_date, end_date);
            end
        end
    end
    
        % make start_date and end_date understandable date strings if output not saved in an output variable
    if (show_dates)
        calendar = neighbors.Properties.UserData.calendar;
        try
            if (calendar_length(calendar) == 360)
                neighbors.startDate = datestr360(neighbors.startDate,'yyyy-mm-dd');
                neighbors.endDate   = datestr360(neighbors.endDate,  'yyyy-mm-dd');
            elseif (calendar_length(calendar) == 365)
                neighbors.startDate = datestr365(neighbors.startDate,'yyyy-mm-dd');
                neighbors.endDate   = datestr365(neighbors.endDate,  'yyyy-mm-dd');
            else
                neighbors.startDate = datestr(neighbors.startDate,'yyyy-mm-dd');
                neighbors.endDate   = datestr(neighbors.endDate,  'yyyy-mm-dd');
            end
        catch
        end
    end
end

