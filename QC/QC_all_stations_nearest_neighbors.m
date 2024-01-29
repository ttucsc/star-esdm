function neighbors = QC_all_stations_nearest_neighbors(stnID, siteTbl, varargin)
% function neighbors = QC_all_stations_nearest_neighbors(stn, sites, nnbrs, maxdist, + kwd/value pairs)
%                                       minpct, dateRange, fullYears, minYears, onlyAll, keepSameLoc, matchID
%
% Locate nearest neighbors to a station or lat/lon
%
%
%       earlier version:  function call was:  function neighbors = QC_all_stations_nearest_neighbors(stn, sites, nnbrs, maxdist, minpct, dateRange, full_or_minyrs, boundRange, onlyAll, keepSameLoc, matchID)

%   stn             station.  Can be:
%                       ix          station index in siteTbl
%                       stnID       station ID (char or string)
%                       [lat,lon]   lat, lon of point
%   sites           siteTbl from netcdf file.  can be:
%                       ncName      name of netcdf file
%                       csvname     name of csv file with fields stnID, stnName, lat,lon [, elev]
%                       siteTbl     table extracted via cal to QC_get_site_tbl
%                                   from all_Locations.csv, 
%                                   or cell array of the above
%   nnbrs           # of neighbors to retrieve (max...may retrieve fewer)
%                       if 0, empty or missing, retrieves all within maxdist
%   maxdist         max dist to neighbors, in km
%                       if 0, empty or missing, maxdist is set 5 degrees (translated to km, using latitude)
%
%   Optional keyword/value pairs:
%       , maxdist, minpct, dateRange, fullYears, minYears, searchType, onlyAll, keepSameLoc, matchID
%
%   minpct          min percentage of valid points
%                       excludes sites with pct_valid < minpct
%                       if 0, empty or missing, does not limit by pct_valid
%                   NOTE:  **** minpct is currently for all records for station, not limited by dateRange ****
%   dateRange       range of dates -- exclude any stations whose dates do not overlap date range
%                       can be [yr1, yr2] or [yr1,mo1,day1;  yr2, mo2, yr2]
%   fullYears       boolean.  if true, include only sites that cover the netire date Range
%   minYears        numeric.   include only stations with at least minYears of data.
%       "fullYears", t/f    if true, include only stations that cover entire date range specified.
%                           if false, include any station with data anywhere in the date range specified.
%       "minYears", #       if > 0, include only stations with at least minYears of data (in date range, if specified)
%       "searchType", typ   "stnID","stnName","any"  
%                           default:  "stnID"
%                           or "like_stnID", "like_stnName", "like_any"
%                               add "like" to match in searchType on any partial match.
%                               otherwise, match is exact match only.
%       "showDates", t/f     if true, dates returned as datestr(...,'yyyy-mm-dd');  else dates returned as datenums
%   onlyAll         0,1,2.    0 include all stations
%                                   1 include only stations with valid data in all siteTables
%                                   2 like 1, but insert empty record for sites with no matches.
%   keepSameLoc     (optional)  if true, keep all records with same lat/lon as last, even if it exceeds the count nnbrs
%                                   This covers cases where there are multiple stations for the same lat/lon.
%   matchID         (optional)  if true, keeps station that matches ID, even if not the closest to the given lat & lon.
%
%       One thing this lacks is selecting nearest station(s) with more than X valid data points.
%       minpct doesn't really give us that, though if pctValid is high enough, this gives us a good proxy.
%       minYears is another good proxy.
%

    siteFiles = ["tmin","tmax","prec","global"];
        
    [stnID, siteTbl, nnbrs, maxdist, searchType, minpct, dateRange, fullYears, minYears, onlyAll, keepSameLoc, matchID] = initParams(stnID, siteTbl, varargin{:})
            
    km_per_deg = 2*pi*6371/360;        % km per degree lat or per degree lon @ equator
    
    if (~isempty_s(minpct) && minpct > 0)
        valid = (siteTbls{itbl}.pctValid >= minpct);
        if (sum(valid) == 0), continue; end
        siteTbls{itbl} = siteTbls{itbl}(valid,:);
    end

    if (exist('dateRange','var') && ~isempty(dateRange) && any(strcmp(siteTbls{itbl}.Properties.VariableNames,"startDate")))
        valid = (datenum(siteTbls{itbl}.startDate) <= datenum(dateRange(2,:)) & datenum(siteTbls{itbl}.endDate) >= datenum(dateRange(1,:)));
        if (sum(valid) == 0), continue; end
        siteTbls{itbl} = siteTbls{itbl}(valid,:);                
    end

    
        % for extending the table if we put more than 200 points in it.
    extenders = create_table({"",0,0}, 200, siteTbls);
    
        % in case they gave us a row of stnIDs or indexes
    if ((isstring(stnID) && isrow(stnID)) || (isnumeric(stnID) && isrow(stnID) && length(stnID)>2)), stnID=stnID'; end
        
    nstns = size(stnID,1);
    neighbors=[];
    
    jj=0;
    for j=1:nstns
        
        jj=jj+1;            % some feedback for user if there are lots of stations.
        if (nstns >=100)
            pct10 = round(nstns/10);
            pct5 = round(nstns/20);
            pct = round(nstns/100);
            if (mod(j,pct10)==0)
                fprintf('%d', round(100*j/nstns));
                jj=0;
            elseif (mod(jj,pct5)==0)
                fprintf('+');
            elseif (mod(jj,pct)==0)
                fprintf('.');
            end
        end
            
        if (isstring(stnID) || ischar(stnID))                                  % stn is siteID to find neighbors for
            if (isstring(stnID))
                stID=stnID(j);
            else
                stID=string(stnID(j,:));
            end
            [stnID, latpt, lonpt] = find_station_by_name(stID, siteTbls);
            if (isempty(stnID))
            % find next open slot in neighbors and insert nbrs into neighbors at that point.
                if (isempty(neighbors))
                    neighbors = create_table({stID, 0, 0}, 1, siteTbls);
                else
                    tn=neighbors;
                    neighbors=[tn;create_table({stID, 0, 0}, 1, siteTbls)];
                end
                
                continue; 
            end
        elseif (istable(stnID))                               % stn is a line from a table w/ lat & lon
            stnID = stnID.stnID(j);
            latpt = stnID.lat(j);
            lonpt = stnID.lon(j);
        elseif (isnumeric(stnID) && size(stnID,2) == 1)         % stn is index of station to find neighbors for
            stnID = siteTbl.stnID(stnID(j));
            latpt = siteTbl.lat(stnID(j));
            lonpt = siteTbl.lon(stnID(j));
        elseif (isnumeric(stnID) && size(stnID,2) == 2)         % stn is lat & lon 
            latpt = stnID(j,1);
            lonpt = stnID(j,2);
            stnID = sprintf('(%.4f,%.4f)', latpt, lonpt);
        else
            error('QC_all_stations_nearest_neighbors:  bad stn info');
        end
        
        nbrs = create_table({stnID, latpt, lonpt}, 20, siteTbls);
        
        for itbl=1:ntbls
            siteTbl = siteTbls{itbl};

            siteTbl.lon = mod(siteTbl.lon+180,360) - 180;
        
            tbl = siteTbl;

            if (maxdist > 0)
                maxlatdegrees = 5*maxdist/km_per_deg;
                maxlondegrees = maxlatdegrees / cos(latpt*pi/180);
            else
                maxlatdegrees = 5;
                maxlondegrees = 5 / cos(latpt*pi/180);
            end

            if (lonpt > 180), lonpt = -360 + lonpt; end  % make sure lonpt is in range -180 to 180

                    % first, shrink the list to something manageable by excluding all the points far from the requested lat/lon.
                    %  bug here for stations close to 180 degrees longitude, ian!
            starters = (abs(tbl.lat-latpt) < maxlatdegrees) & (abs(tbl.lon - lonpt) < maxlondegrees);

            tbl = tbl(starters,:);

            nstations = size(tbl,1);

            latpts = repmat(latpt,nstations,1);
            lonpts = repmat(lonpt,nstations,1);


            [dist, az] = distance(latpts, lonpts, tbl.lat, tbl.lon); 
            tbl.dist = round(dist * km_per_deg, 3);        % distance in km from point, rounded to nearest meter.
            tbl.az = az;

            if (isempty(tbl)), continue; end

            if (~isempty_s(maxdist) && maxdist > 0)
                keepers = tbl.dist <= maxdist;
                if (sum(keepers) == 0), continue; end
                tbl = tbl(keepers,:);
            end
            
            if (onlyAll)
                keepers = tbl.pctValid > 0;
                if (sum(keepers) == 0), continue; end
                tbl = tbl(keepers,:);
            end
    
            nbrs = insert_neighbor(stnID, j, nbrs, tbl, itbl, extenders); 
        end
        
            % trim nbrs down to only the used part.
        next_open = find(strlength(nbrs.nbrID)==0, 1);
        if (~isempty(next_open)), nbrs(next_open:end,:) = []; end
        nbrs = sortrows(nbrs,{'dist','nbrID'});
        
        if (onlyAll>0)
            keepers = onlyAlls(nbrs);
            if (sum(keepers) > 0)
                nbrs = nbrs(keepers,:);
            elseif (onlyAll == 1)
                continue;
            else
                nbrs = create_table({stnID, 0, 0}, 1, siteTbls);
                nbrs{1,1}=j;
            end
        end

        
        if (~isempty_s(nnbrs) && nnbrs > 0 && size(nbrs,1) > nnbrs)
            if (keepSameLoc)
                mnbrs = find(nbrs.dist <= nbrs.dist(nnbrs)+.012,1,'last');      % .012km  is slightly over .0001 degrees, so will catch stations where location matches to nearest .0001 degrees.
            else
                mnbrs = nnbrs;
            end
            nbrs = nbrs(1:mnbrs,:); % nbrs(nnbrs+1:end,:)=[];
        end
        
        if (matchID)        % if matchID, keep only stations with matching IDs (there should be at most 1)
            ix = find(nbrs.nbrID == stnID);
            if (~isempty(ix))
                if (islogical(matchID) || matchID == 1)
                    nbrs = nbrs(ix,:);
                else
                    nbrs = nbrs(1:ix,:);
                end                    
            end
                            % if not match, leave all stations in the list.
        end
        
            % find next open slot in neighbors and insert nbrs into neighbors at that point.
        if (isempty(neighbors))
            neighbors = nbrs;
        else
            tn=neighbors;
            neighbors=[tn;nbrs];
        end
    end
    
        % add some flags for problem records
    nrecs = size(neighbors,1);
    flags = strings(nrecs, 1);
    
        % check for dups.
            % multiple matches
    dups = neighbors.n(1:end-1) - neighbors.n(2:end);
    fd = find(dups==0);
    dups(fd+1)=0;
    flags(dups == 0) = "m";
            % different stnID than given
    difs = neighbors.stnID ~= neighbors.nbrID;
    flags(difs) = flags(difs)+"D";
            % distance > 1 km
    dist = neighbors.dist > 1;
    flags(dist) = flags(dist)+"+";
    
    t=table(flags);
    neighbors = [t, neighbors];
    
    
    
        % and display it if output not being captured to a variable.
    if (nargout==0)
        disp(neighbors);
    end

end

function [stnID, latpt, lonpt] = find_station_by_name(stn, siteTbls)

    for i=1:length(siteTbls)
        stnix = find(siteTbls{i}.stnID == string(stn));

        if (~isempty(stnix)), break; end

        stnix = find(siteTbls{i}.stnName == string(stn));
        if (~isempty(stnix)), break; end
    end
    if (~isempty(stnix))
        stnID = siteTbls{i}.stnID(stnix);
        latpt = siteTbls{i}.lat(stnix);
        lonpt = siteTbls{i}.lon(stnix);
        return;
    end

        % no exact match.  look for substring in station names
    for i=1:length(siteTbls)
        stnix = find(contains(siteTbls{i}.stnName,string(stn),'IgnoreCase',true),1);
        if (~isempty(stnix)), break; end
    end
    if (~isempty(stnix))
        stnID = siteTbls{i}.stnID(stnix);
        latpt = siteTbls{i}.lat(stnix);
        lonpt = siteTbls{i}.lon(stnix);
    else
        stnID=[];
        latpt = nan;
        lonpt = nan;
    end
    
    
    
end
    

function  neighbors = create_table(stn, nstns, siteTbls)
    % make room for 20 matches.  Matlab will expand beyond 20 if needed.
    % we'll truncate to the actual number of matches later,
    
    ntbls   = length(siteTbls);
    n       = zeros(nstns,1);
    stnID   = strings(nstns,1);
    nbrID   = strings(nstns,1);
    stnName = strings(nstns,1);
    lat     = nan(nstns,1);
    lon     = nan(nstns,1);
    elev    = nan(nstns,1);
    state   = strings(nstns,1);
    dist    = nan(nstns,1);
    az      = nan(nstns,1);
    xtra    =cell2table(cell(nstns,5*ntbls));
        
    neighbors = [table(n, stnID, nbrID, stnName, lat, lon, elev, state, dist, az),xtra];

    if (~isempty(stn))
        neighbors.stnID(1)=stn{1};
        neighbors.lat(1)=stn{2};
        neighbors.lon(1)=stn{3};
    end
    
        % there's probably a better way of adding columns for each table without getting the warning about
        % inefficiency of table growing on each loop.  But this works.
    startDate = repmat("NA", nstns, 1); 
    endDate   = repmat("NA", nstns, 1); 
    pctValid  = nan(nstns, 1);
    years     = nan(nstns, 1);
    ixcol     = zeros(nstns, 1);
    t=table(startDate,endDate,pctValid,years,ixcol);
    vnt=t.Properties.VariableNames;
    for i=1:length(siteTbls)
        for j=1:5
            ixn=5*(i-1)+j;
            vnn=sprintf('Var%d',ixn);
            neighbors.(vnn) = t.(vnt{j});
        end
        if (isfield(siteTbls{i}.Properties.UserData,'varName'))
            vname=siteTbls{i}.Properties.UserData.varName;
            vnames={sprintf('%sStart',vname),sprintf('%sEnd',vname),sprintf('%sPctVld',vname),sprintf('%sYrs',vname),sprintf('%sIx',vname)};
        else
            vnames={sprintf('startDate%d',i),sprintf('endDate%d',i),sprintf('pctValid%d',i),sprintf('Years%d',i),sprintf('ix%d',i)};
        end
        neighbors.Properties.VariableNames(10+5*(i-1)+(1:5))=vnames;

    end
end

function neighbors = insert_neighbor(stnID, stnix, neighbors, tbl, itbl, extenders)
    % inserts tbl into neighbors, either matching existing or at next available slot.

    if (strlength(stnID) == 0), stnID = neighbors.stnID(1); end
    nbrlen = size(neighbors,1);
    vns=neighbors.Properties.VariableNames;
    tbllen = size(tbl,1);
    
    ix=find(strlength(neighbors.nbrID)==0,1);
    no_nbrs = ~isempty(ix) && ix==1;
    
    for i=1:tbllen
        if (no_nbrs)
            ix = 1;
            no_nbrs = false;
        else
            ix = find(neighbors.stnID == stnID & neighbors.nbrID == tbl.stnID(i));
        end
        if (isempty(ix))
            next_open = find(neighbors.n==0, 1);
            if (isempty(next_open))
                tmp = neighbors;
                neighbors = [tmp; extenders];
                next_open = nbrlen+1; 
            end
            ix=next_open; 
        end
        neighbors.n(ix)       = stnix;
        neighbors.stnID(ix)   = stnID;
        neighbors.nbrID(ix)   = tbl.stnID(i);
        neighbors.stnName(ix) = tbl.stnName(i);
        neighbors.lat(ix)     = tbl.lat(i);
        neighbors.lon(ix)     = tbl.lon(i);
        neighbors.dist(ix)    = tbl.dist(i);
        neighbors.az(ix)      = tbl.az(i);
        
        stix = find(strcmpi(tbl.Properties.VariableNames, "State"),1);
        if (~isempty(stix)), neighbors.state(ix) = tbl{i,stix}; end
        elix = find(strcmpi(tbl.Properties.VariableNames, "elev"),1);
        if (~isempty(elix)), neighbors.elev(ix) = tbl{i,elix}; end

        
        stdix = 10+(itbl-1)*5+1;
        endix = stdix+1;
        pctix = stdix+2;
        yrix  = stdix+3;
        ixix  = stdix+4;
        
        
            % if the table is one of our station tables, it will have the following columns.
            
        colix=find(strcmp(tbl.Properties.VariableNames, "startDate"),1);
        if (~isempty(colix))
            if (isnumeric(tbl.startDate(i)))
                neighbors.(vns{stdix})(ix) = datestr(tbl.startDate(i),'yyyy-mm-dd');
            else
                neighbors.(vns{stdix})(ix) = tbl.startDate(i);
            end
        else
            neighbors.(vns{stdix})(ix) = "unknown";
        end
        colix=find(strcmp(tbl.Properties.VariableNames, "endDate"),1);
        if (~isempty(colix))
            if (isnumeric(tbl.endDate(i)))
                neighbors.(vns{endix})(ix) = datestr(tbl.endDate(i),'yyyy-mm-dd');
            else
                neighbors.(vns{endix})(ix) = tbl.endDate(i);
            end
        else
            neighbors.(vns{endix})(ix) = "unknown";
        end
        colix=find(strcmp(tbl.Properties.VariableNames, "pctValid"),1);
        if (~isempty(colix)), neighbors.(vns{pctix})(ix) = tbl.pctValid(i); end
        colix=find(strcmp(tbl.Properties.VariableNames, "years"),1);
        if (~isempty(colix)), neighbors.(vns{yrix})(ix) = tbl.years(i); end
        colix=find(strcmp(tbl.Properties.VariableNames, "index"),1);
        if (~isempty(colix)), neighbors.(vns{ixix})(ix)  = tbl.index(i); end
    end
end

function keepers = onlyAlls(nbrs)
% returns boolean array flagging rows where all siteTbls have valid data.

    n=size(nbrs,1);
    keepers = true(n,1);
    
    ntbls = (length(nbrs.Properties.VariableNames)-10)/5;
    
    for ix=1:n
        for i=1:ntbls
            colix = 10 + 4 + 5*(i-1);
            if (isnan(nbrs{ix,colix}) || nbrs{ix,colix} == 0)
                keepers(ix)=false;
                continue;
            end
        end
    end
end

function fullname = locate (siteFile)

    [~,hostname]=system('hostname');
    hostname = strtrim(hostname);
    on_dev_system  = contains(hostname,"icsf");
    on_hpcc_system = contains(hostname,"compute");
    
    if (on_dev_system)
        mypath = './';          % these need to be single quotes.  fullfile doesn't like strings yet.
    elseif (on_hpcc_system)
        mypath = '/lustre/scratch/iscottfl/obs/stations_netcdf';
    else
        mypath = '/data/obs/stations_netcdf';
    end
    
    siteFiles=["Tmin","tmax","prec","global"] ;
    ix = find(strcmpi(siteFiles, siteFile));
    
    if (isempty(ix)), error("error:  invalid siteFile:  %s\n", siteFile); end
       
    if (ix == 4)
        fullname = fullfile(mypath, 'all_locations.csv');
    else
        fullname = fullfile(mypath, sprintf('stations.%s.ncamerica.1850.2017.nc',siteFiles(ix)));
    end
end


function [stnID, siteTbl, nnbrs, maxdist, searchType, minpct, dateRange, fullYears, minYears, onlyAll, keepSameLoc, matchID] = initParams(stnID, siteTbl, varargin)

    search_types = ["stnID","stnName","any", "like_stnID", "like_stnName", "like_any"];
                    % parse input for DA_title
    p = inputParser;
    p.StructExpand = true;
    p.KeepUnmatched = true;        

    addRequired(p,"stnID",                  @(s) ischars(s));
    addRequired(p,"siteTbl",                @(s) ischars(s) || (istable(s) && isQCstntbl(tbl)));
    addOptional(p,"nnbrs",        25,       @(s) isnumeric(s) && s>0);
    addOptional(p,"maxdist",      0,        @(s) isnumeric(s) && s>=0);
    addParameter(p,"searchType", "stnID",   @(s) ischar_s(s) && any(strcmpi(s,search_types)));
    addParameter(p,"minpct",      0,        @(s) isnumeric(s) && s>=0 && s<=100);
    addParameter(p,"dateRange",   [],       @(s) isempty(s) || isnumeric(s) && any(numel(s)==[2,6]));
    addParameter(p,"fullYears",   false,    @(s) islogical(s) || (isnumeric(s) && any(s==[0,1])));
    addParameter(p,"minYears",    0,        @(s) numeric(s) && s>=0);
    addParameter(p,"onlyAll",     0,        @(s) isnumeric(s) && any(s==[0,1,2]));
    addParameter(p,"keepSameLoc", false,    @(s) islogical(s) || (isnumeric(s) && any(s==[0,1])));
    addParameter(p,"matchID",     false,    @(s) islogical(s) || (isnumeric(s) && any(s==[0,1])));
    addParameter(p,"showDates",  false,     @(s) islogical(s) || (isnumeric(s) && any(s==[0,1])));
    addParameter(p,"loadTblVars", true,     @(s) islogical(s) || (isnumeric(s) && any(s==[0,1])));

    parse(p, stnID, siteTbl, varargin{:});
    Parms = p.Results;
    Unmatched = p.Unmatched;        % save rest of input params for later
    
    if (~isempty(fields(Unmatched)))
        msg = strcat("error:  unexpected input keywords:  \n\t", sprintf(2,"%s ", fields(Unmatched)));
        error(msg);
    end
    
    stnID       = Parms.stnID;
    siteTbl     = Parms.siteTbl;
    nnbrs       = Parms.nnbrs;
    maxdist     = Parms.maxdist;
    searchType  = Parms.searchType;
    minpct      = Parms.minpct;
    dateRange   = Parms.dateRange;
    fullYears   = Parms.fullYears;
    minYears    = Parms.minYears;
    onlyAll     = Parms.onlyAll;   
    keepSameLoc = Parms.keepSameLoc;   
    matchID     = Parms.matchID;   
    showDates   = Parms.showDates;
    loadTblVars = Parms.loadTblVars;
            
    if (ischar_s(siteTbl))
        siteTbl = QC_get_site_table(siteTbl,start_date, end_date, 'fullYears',fullYears, "minYears", minYears, "loadTblVars",loadTblVars, "showDates", showDates);
    else
        if (~isQCstntbl(siteTbl)), error("error:  input table is not a QC_stntbl"); end
        calendar = siteTbl.Properties.UserData.calendar;
        if (~isempty(start_date) || ~isempty(end_date)) 
%             if (isempty(start_date)), my_start_date = siteTbl.Properties.UserData.day1; end
%             if (isempty(  end_date)),   my_end_date = siteTbl.Properties.UserData.day1 + siteTbl.Properties.UserData.npts-1; end
            if (isempty(start_date)), my_start_date = siteTbl.Properties.UserData.dates(1); end
            if (isempty(  end_date)),   my_end_date = siteTbl.Properties.UserDatadates(end); end
            if (~fullYears)
                keepers = (siteTbl.startDate <= my_end_date & siteTbl.endDate >= my_start_date);
            else
                keepers = (siteTbl.startDate <= my_start_date & siteTbl.endDate >= my_end_date);
            end
            if (sum(keepers) == 0)
                error('no valid data in date range %s to %s', datestr_cal(my_start_date,calendar, 'yyyy-mm-dd'),datestr_cal(my_end_date,calendar, 'yyyy-mm-dd')); 
            end
            siteTbl = siteTbl(keepers,:);
            siteTbl.Properties.UserData.allSites = false;
        end
        if (minYears > 0)
            bounded_startDate = max(siteTbl.startDate, start_date);
            bounded_endDate   = min(siteTbl.endDate,   end_date);

                % years rounded so leap-year issues don't bite us.  If divide by 365.25, but we're short 1 leap-day,
                % then test of >= nyears might fail incorrectly otherwise.

            stn_years = round((bounded_endDate - bounded_startDate)/calendar_length(calendar),2);
            keepers = stn_years >= minYears;       % exclude sites with less than minYears data
            siteTbl = siteTbl(keepers,:);
            siteTbl.Properties.UserData.allSites = false;
        end

    end
    
    if (~isempty(dateRange))
        if (length(dateRange)==2)
            dateRange = [dateRange(1),1,1; dateRange(2),12,31];
        end
    end
        
    if (ischar_s(stnID) && isfile(stnID)) 
        fname = stnID; 
        stnID = readtable(stnID);
    end
    
    if (istable(stnID))
        if (isfield(stnID.Properties.VariableNames,"stnID"))
            stnID = stnID.stnID;
            searchType = "stnID";
        elseif (isfield(stnID.Properties.VariableNames,"stnName"))
            stnID = stnID.stnName;
            searchType = "stnName";
        else
            stnID = stnID(:,1);     % No stnID or stnName in table.  use 1st column.  
            if (~ischars(stnID)), error("error:  No stnID or stnName column, and 1st column is not strings");end
            fprintf("using 1st column of %s for stnIDs", fname);
        end
    end
end