function [outTbl] = QC_to_netcdf_hourly_daily(varNames, stnInfo, yrRange, varargin)
% Program to convert csv or txt data to netcdf files.
%   written originally for Anne's hourly data for her 60 stations to netcdf files, but supports daily as well as hourly
%   data.
%
%   varNames        Variable name or array of variable names for csv column(s) with data to process. 
%                       Can be cell or string array of variable names,if the file contains more than one variable
%   stnInfo         either name of file with list of sites ('station_list.csv') or table created from it.
%                       file must contain the following columns
%                           stnID       station ID (often GHCN station ID
%                           stnName     station name 
%                           lat         lat,lon of station
%                           lon 
%                           filename    name of file with data
%
%                Possible keyword/value pairs of additional parameters
%
%               Many of these are given for specific datasets in function getParams.  If so, no need to specify them.
%               "Presets" for several known input sets are given in the function getParams(...), and will override
%               entries given here (sorry...change the presents as needed!)
%               See function getParams for more details.  
%
%   yrRange         start, end year of data to keep                 (must be specified;  cannot be -inf or inf
%                       given as [startyr,endyr]
%   min_yrs         mininum # of years of data required to include station in file. [default: 25]
%   basedir         base folder name where data files are located     default: "."
%   outdir          base folder of where to put output file.          default: "."
%   stnIDColumn     name of column  in stnInfo with stnIDs.           default: "stnID"
%   stnNameColumn   name of column  in stnInfo with stnNames.         default: "stnName"
%   latColumn       name of columns in stnInfo with lats & lons       default: "lat" and "lon"
%   lonColumn
%   filenameColumn  name of column  in stnInfo with filenames.        default: "filename"
%   elevColumn      name of column  in stnInfo with elevations.       default: "elev"
%   TZColumn        name of column  in stnInfo with Time Zone info.   default: "TZ"  
%   timeColumn      cell array of time column names in data file     default: ["year","month","day"]
%   ncComments      struct with keyword/string pairs of info to put in netcdf global comments, or
%                       name of a file containing keyword/value pairs, or
%                       a single string or char array containing 1 comment to add to netcdf as "comments"
%   units           units for input data (string array, 1 for each varName)
%   outUnits        units for output netcdf file (string array, one for each varName)
%                       note:  will use units and outUnits to determine toDegF, toDegC, toInches and toMM if not given.
%   long_names       long name for output units (string array, one for each varName)
%   outVarNames     output variable names (names in netcdf file).   (string array, one for each varName)
%   ncName          output netcdf filename
%   isHourly        boolean flag.  If true, data is hourly.
%   NA_in           string or value in input file representing NA (nan)  default:  "NA"
%                       data written out will have the standard 1e20 in output netcdf file.
%   local_ghcn      boolean.  if true, use local ghcn stations (only the ones in this file) to find nearest GHCN station
%   ghcn_stninfo    stninfo or stninfo filename to use for nearest GHCN stations
%   
%                       6/2019:  added the following to work with Anne's merra data
%   isUTC          true/false.  Use TZ info from station_info file to adjust times to UTC for netcdf of local-time data.
%   toUTC          true/false.  add TZ to time to adjust to UTC if isUTC is true.
%   source          "merra" or blank, for now.  Used to specify columns, units, etc.
%
%   max_stations    for debugging, stop processing after this # of stations.
%   istep           alternately, for debugging, # of stations to skip.
%
%   toInches, toMM, toDegF, toDegC  logical.  set to true to override preassigned output units where needed.
%                       
%                       
%
%   returns augmented QC station table
%   writes netcdf file to ncName
%
%   This program will generally need a little editing to add support for new csv/txt input files.
%
%   For now, supports Anne Stoner's station files for her transportation studies and lake (ntmwd).
%
%
%
%--------------TO DO----------------
%   This needs an outCalendar to specify whether to create 360-day, 365-day or standard calendar output.
%   For now, this always write a standard calendar netcdf.
%
%   This needs an inCalendar so it can handle csv's of non-standard (leap-day) calendars.
%
%   needs automatic calculation of yrRange
%
%   This needs a way to put multiple variables for the same stations into a single netcdf file
%   This needs a way to add or update a variable to an existing netcdf file.
%   This Needs automatic lookup of timezone for input that doesn't supply timezone info 
%       Google API now charges for large numbers of timezone calls;  timezoneDB API is slow, so this will slow it down.
%       Probably want to add "timezone_lookup" flag.
%
%---------------------QC Table output----------------
%   QC table (output) columns:
%           stnID           station ID
%           lat, lon        location  (decimal degrees, usually to 4 decimals;  some only to 2 or 3)
%           stnName         station Name
%           elev            elevation, in meters
%           startDate       start timestamp of first valid record for site, as Matlab datenum (for data's calendar type)
%           endDate         end   timestamp of last  valid record for site, as Matlab datenum (for data's calendar type)
%           pctValid        percent of valid points between startDate and endDAte
%           index           index of record in original input file
%           time_zone       (currently set to NA)
%           is_ghcn         boolean.  true if site is a GHCN station
%           nearest_ghcn    stnID of nearest GHCN station, if is_ghcn is false
%           ghcn_elev       elevation of nearest GHCN stations
%           ghcn_dist       distance to nearest GHCN station
%           ghcn_az         direction (degrees AZ from true North) to nearest GHCN station.
%
%   QC table Properties.UserData fields & example values:
%       Properties.UserData contains useful table-wide data extracted from netcdf file's global attributes.
%
%         ncName: "/Volumes/lacie_1/jcsf_data/data/obs/stations_netcdf/stations.Tmax.ncamerica.1850.2018.20_yrs_min.nc"
%                                                       original input file, full pathname
%        varName: 'TMAX'                                variable name
%           npts: 24106                                 # of data points
%      daterange: "1950-01-01 2015-12-31"               date range of extracted data
%      fullRange: 0                                     boolean.  flags whether all dates of original file included
%                                                           if false, data is extracted date range from original file
%       minYears: 0                                     if>0 only keeping stations with at least this # of years of data
%     boundRange: 0                                     boolean. If true, 
%         NAFlag: 1.0000e+20                            file's original NA value ('missing or invalid data' flag)
%       longName: 'daily max temp'                      variable long name
%          units: 'degC'                                variable units
%         nsites: 70                                    # of sites being returned
%      timeunits: 'days since 1850-01-01 00:00:00 UTC'  netcdf time units
%       calendar: 'standard'                            calendar type.  365-day, 360-day, standard, etc.
%           day1: 675700                                matlab datenum of 1st day 
%          isUTC: 1                                     boolean. flags whether timestamps are UTC or local time.
%     isQCstntbl: 1                                     boolean.  true if table is a QC station table.
%          dates: [24106×1 double]                      timesteps, as datenum_cal's (datenum based on calendar)
%
%-----------------------------------

        % RP is RunParams struct.
    Parms = getParams(varNames, yrRange, varargin{:});
    nvars = length(Parms.varNames);
    
    %-------------------variable and column identification.  Update this for your specific files.
    
        % these are from Anne's Transportation Study files
        
    if (isempty(stnInfo)), stnInfo = "station_list.csv"; end
    if (~istable(stnInfo))
        stnInfo = readtable(stnInfo);
    end
    
    workTbl = add_QC_columns(stnInfo, Parms);
    NAFlag = 1e20;

    %------------------------Input, Output filenaming-----------------Update this for your specific files----------
    if (Parms.isHourly)
         timestr="Hourly";
    else
        timestr="Daily";
    end
%   fprintf('creating netcdf\n');
    [nc, tstamps] = create_ncdf(Parms.ncName, workTbl, Parms.yrRange, Parms.min_yrs, Parms.outVarNames, Parms.long_names, Parms.outUnits, NAFlag, Parms.ncComments, timestr, Parms.local_ghcn);  
%   fprintf('done creating netcdf\n');

    %---------------------Output filename
%     fprintf('creating netcdf file %s\n', outname);
%     fprintf('done creating nc file\n');
    nstns = length(workTbl.stnID);
    empty = false(nstns,1);


    %---------------------Process and add each station.

    iloop = 0;
    nloops = nstns * nvars;
    nstns_out = 0;
    keepers = false(nstns,1);
    if (~isempty(Parms.max_stations))
        istep = floor(nstns/(Parms.max_stations-1));
    else
        Parms.max_stations = ceil(nstns/Parms.istep);
        istep = Parms.istep;
    end
    
    fprintf('.........QC_to_netcdf %s,   %d to %d,  %s\n\n',Parms.ncName, yrRange(1), yrRange(2), datestr(now,'yyyy-mm-dd HH:MM:SS'));
    for i=1:nvars
        if (istep > 1)
            fprintf('Variable:  \t%s  -> %s (%d stations only, stepping by %d stations\n', Parms.varNames{i}, Parms.outVarNames{i}, Parms.max_stations, istep); 
        else
            fprintf('Variable:  \t%s  -> %s\n', Parms.varNames{i}, Parms.outVarNames{i}); 
        end
    end        

    for istn=1:istep:nstns            % istep is for for debugging
        for iv=1:nvars
            iloop = iloop + 1;
            ffnm = workTbl.filename{istn};
            fnm = fullfile(Parms.basedir,ffnm);

            if (~isfile(fnm))
                fprintf('%s:  no data for %s\n', ffnm, Parms.varNames{iv}); 
                continue; 
            end
            [vals, min_date, max_date, nvalid, mpts, file_startdate,file_enddate] = read_csv(fnm, Parms.varColumn{iv}, Parms.timeColumn, tstamps,timestr, Parms.NA_in, Parms.isUTC, Parms.toUTC, workTbl.time_zone(istn));
            nyrs = (max_date-min_date)/365;
            if (isempty(file_startdate) || isempty(file_enddate))
                empty(istn)=true;
                fprintf('%6d of %6d (%6.2f):           no valid data for %s\n', iloop, nloops, 100.0*iloop/nloops, ffnm); 
                continue;
            elseif (nvalid == 0)
                empty(istn)=true;
                fprintf('%6d of %6d (%6.2f):           no data in range for %s.  File dates: %s to %s\n', iloop, nloops, 100.0*iloop/nloops, ffnm, datestr(file_startdate,'yyyy-mm-dd HH:MM:SS'), datestr(file_enddate,'yyyy-mm-dd HH:MM:SS')); 
                 continue;
            elseif (nyrs < Parms.min_yrs)
                fprintf('%6d of %6d (%6.2f):           insufficient data for %s.  File dates: %s to %s\n', iloop, nloops, 100.0*iloop/nloops, ffnm, datestr(file_startdate,'yyyy-mm-dd HH:MM:SS'), datestr(file_enddate,'yyyy-mm-dd HH:MM:SS')); 
                continue;
            end
            workTbl.startDate(istn,iv)= min_date;
            workTbl.endDate(istn,iv)  = max_date;
            workTbl.pctValid(istn,iv) = nvalid / mpts * 100.0;       
            nstns_out = nstns_out+1;
            keepers(istn,1) = true;
            fprintf('%6d of %6d (%6.2f):  %6d: %-10s %20s %-32s %s to %s  %7.3f %% valid  from file %s\n', iloop, nloops, 100.0*iloop/nloops, nstns_out, varNames(iv), workTbl.stnID{istn}, workTbl.stnName{istn}, datestr(min_date,'yyyy-mm-dd'), datestr(max_date,'yyyy-mm-dd'), workTbl.pctValid(istn, iv), ffnm);

            vals=jc_units_conversion(vals, Parms.units(iv), Parms.outUnits(iv));

            write_stn_data(Parms.outVarNames(iv), vals, istn, nstns_out, nc, workTbl, Parms.yrRange(1), iv, Parms.local_ghcn);
        end
    end
    
    workTbl = workTbl(keepers,:);
    
    if (Parms.local_ghcn)
        workTbl = find_closest_ghcns(workTbl);
        rewrite_nearest_ghcn_data(nc, workTbl)
    else
        workTbl = find_closest_ghcns(workTbl, Parms.ghcn_stninfo);
        rewrite_nearest_ghcn_data(nc, workTbl)        
    end

    fprintf('%d stations written to netcdf file %s \n', nstns_out, Parms.ncName);
    if (nargout > 0)
        outTbl = workTbl(keepers,:);
    else
        outTbl = [];
    end
end

function [nc, tstamps] = create_ncdf(outName, workTbl, yrRange, min_yrs, outVarNames, long_names, outUnits, NAFlag, ncComments,timestr, local_ghcn)

% timestr is "Daily" or "Hourly"

    nc = ncdf('', 'Filename',outName, 'Format','netcdf4','create_ok',true);      % create an empty ncdf
    nc.putatt('title',sprintf('TTU %s Station Data, created %s', timestr, datestr(now,'yyyy-mm-dd')));
    nc.putatt('institution','Texas Tech Climate Science Center');
    nc.putatt('description',sprintf('TTU %s station data for %s', timestr, join(outVarNames,', ')));
    nc.putatt('date_range',sprintf('%s to %s', datestr(datenum(yrRange(1),1,1),'yyyy-mm-dd'), datestr(datenum(yrRange(2),12,31),'yyyy-mm-dd')));
    if (min_yrs > 1)
        nc.putatt('min_years',sprintf('including only stations with minimum of %d years in date range', min_yrs));
    end
    [~,uname] = getusername();
    nc.putatt('creation_date',datestr(now,'yyyy-mm-dd HH:MM:SS'));
    nc.putatt('creator',uname);
    nc.putatt('source_code',mfilename);
    if (strcmpi(timestr,"hourly"))
        nc.putatt('hourly_station_data','yes');
        nc.putatt('daily_station_data','no');
    else
        nc.putatt('hourly_station_data','no');
        nc.putatt('daily_station_data','yes');
    end
    if (local_ghcn)
        nc.putatt('nearest_ghcn_info','nearest_ghcn is local to the file (is present in this file)');
    else
        nc.putatt('nearest_ghcn_info','nearest_ghcn is from full GHCN database, and station may not be present in file');
    end
    nc.putatt('data_variables',join(outVarNames,'|'));

    if (isstruct(ncComments))
        fn=fieldnames(ncComments);
        for i=1:length(fn)
            nc.putatt(fn{i}, ncComments.(fn{i}));
        end
    else
        if (length(ncComments)==1)
            nc.putatt('comment', ncComments{1});
        else
            for i=1:length(ncComments)
                comment_lbl=sprintf("comment%d",i);
                nc.putatt(comment_lbl,ncComments{i});
            end
        end
    end    
    
%     nstations = size(workTbl,1);

    ndays = daysdif(datenum([yrRange(1),1,1]),datenum([yrRange(2)+1,1,1]));
    day1 = datenum([yrRange(1), 1, 1]);
    if (strcmpi(timestr,"hourly"))
        timevals    = (0:(24*ndays-1))/24;
    else
        timevals    = (0:(ndays - 1));
    end
    ntimes = length(timevals);
    tstamps = day1 + timevals;
    
            % char(...) conveniently space-pads the strings to the same length.
%   stnIDs = (char(workTbl.stnID))';      % transpose needed here to order the chars properly when writing out!
    stnIDLen = max(strlength(workTbl.stnID));

%   stnNames = (char(workTbl.stnName))';
    stnNameLen = max(strlength(workTbl.stnName));
    
%   nearest_ghcn = (char(workTbl.nearest_ghcn))';
    nearest_ghcnLen = max(11, max(strlength(workTbl.nearest_ghcn))); % size(nearest_ghcn,1);

    nvars=length(outVarNames);
    coutVarNames = (char(outVarNames))';
    varNameLen = size(coutVarNames,1);
    
    
        % define dimensions
    stn_dim = Dimension('stn_num',Inf,true);      % make stn_dim an unlimited dimension.
    time_dim = Dimension('time',ntimes);
    idLen_dim = Dimension('stnIDLength',stnIDLen);
    nmLen_dim = Dimension('stnNameLength',stnNameLen);
    ghcnLen_dim = Dimension('ghcnIDLength',nearest_ghcnLen);
    varNameLen_dim = Dimension('varNameLength',varNameLen);
    nvar_dim = Dimension('nvars',nvars,true);     % make nvar_dim unlimited, so we can add other variables later.
    
    nc.putdim(stn_dim);
    nc.putdim(time_dim);
    nc.putdim(idLen_dim);
    nc.putdim(nmLen_dim);
    nc.putdim(ghcnLen_dim);
    nc.putdim(varNameLen_dim);
    nc.putdim(nvar_dim);
    
        % Create Variables
        
        % time
    timeUnits = sprintf('days since %s UTC', datestr(day1,'yyyy-mm-dd HH:MM:SS'));
    calendar='standard';
    
    timeVar = Variable('time', timevals,'Dimensions',time_dim);
    timeVar.putatt('units',timeUnits);
    timeVar.putatt('calendar',calendar);
    timeVar.putatt('UTC_or_local','UTC');
    
        % latitude, longitude, elevation, timezone
        
    latVar = Variable('lat', 'Dimensions',stn_dim,'Datatype','double','FillValue',NAFlag);
    latVar.putatt('long_name','latitude');
    latVar.putatt('standard_name','latitude');
    latVar.putatt('units','degrees north');
    lonVar = Variable('lon','Dimensions',stn_dim,'Datatype','double','FillValue',NAFlag);
    lonVar.putatt('long_name','longitude');
    lonVar.putatt('standard_name','longitude');
    lonVar.putatt('units','degrees east');
    elevVar = Variable('elevation', 'Dimensions',stn_dim,'Datatype','single','FillValue',single(NAFlag));
    elevVar.putatt('units','meters');
    elevVar.putatt('long_name','elevation, meters above sea level');
    tzVar = Variable('time_zone','Dimensions',stn_dim,'Datatype','single','FillValue',single(NAFlag));
    tzVar.putatt('units','hours');
    tzVar.putatt('long_name','TZ hours east of GMT');
   

    st_dateVar = Variable('start_date','Dimensions',{stn_dim, nvar_dim},'Datatype','double','FillValue',NAFlag);
    st_dateVar.putatt('units',timeUnits);
    end_dateVar = Variable('end_date','Dimensions',{stn_dim, nvar_dim},'Datatype','double','FillValue',NAFlag);
    end_dateVar.putatt('units',timeUnits);
    pct_validVar = Variable('pct_valid','Dimensions',{stn_dim, nvar_dim},'Datatype','single','FillValue',single(NAFlag));
    pct_validVar.putatt('units','percent');

    
        % stations
        
    stnIDVar = Variable('stnID','Dimensions',[idLen_dim, stn_dim],'Datatype','char','FillValue',' '); 
    stnIDVar.putatt('long_name','Station ID');
    stnIDVar.putatt('description',sprintf('Station ID, stored as fixed length char strings %d chars long, padded with trailing spaces',stnIDLen));
    stnNameVar = Variable('stnName','Dimensions',[nmLen_dim, stn_dim],'Datatype','char','FillValue',' ');
    stnNameVar.putatt('long_name','Station Name');
    stnNameVar.putatt('description',sprintf('Station Name, stored as fixed length char strings %d chars long, padded with trailing spaces',stnNameLen));
    
    varNameVar = Variable('varName',coutVarNames,'Dimensions',[varNameLen_dim, nvar_dim],'FillValue',' '); 
    varNameVar.putatt('long_name','Data Variable Names');
    varNameVar.putatt('description',sprintf('Variable Names for data variables, stored as fixed length char strings %d chars long, padded with trailing spaces',varNameLen));
    
    is_ghcnVar = Variable('is_ghcn','Dimensions',stn_dim,'Datatype','uint8');
    is_ghcnVar.putatt('description','boolean, 1 if station is a GHCN station, 0 if not');
    
    nearest_ghcnVar = Variable('nearest_ghcn','Dimensions', [ghcnLen_dim, stn_dim],'Datatype','char','FillValue',' ');
    nearest_ghcnVar.putatt('long_name','nearest_ghcn_station');
    if (local_ghcn)
        nearest_ghcnVar.putatt('description','Station ID of nearest GHCN station in this file');
    else
        nearest_ghcnVar.putatt('description','Station ID of nearest TTU QC''d observation data from 2018 GHCN files');
    end
       
    
    ghcnDistVar = Variable('ghcn_distance', 'Dimensions',stn_dim,'Datatype','single','FillValue',single(NAFlag));
    ghcnDistVar.putatt('long_name','nearest_ghcn_distance');
    if (local_ghcn)
        ghcnDistVar.putatt('description','distance, in km, to nearest GHCN station in this file');
    else
        ghcnDistVar.putatt('description','distance, in km, to nearest TTU QC''d observation data from 2018 GHCN files');
    end
    ghcnDistVar.putatt('units','km');
    
    ghcnAzVar = Variable('ghcn_az', 'Dimensions',stn_dim,'Datatype','single','FillValue',single(NAFlag));
    ghcnAzVar.putatt('long_name','ghcn_az');
    ghcnAzVar.putatt('description','direction (from true North) from station to nearest GHCN station');
    ghcnAzVar.putatt('units','degrees from true north');
       
    ghcnElevVar = Variable('ghcn_elev', 'Dimensions',stn_dim,'Datatype','single','FillValue',single(NAFlag));
    ghcnElevVar.putatt('long_name','ghcn_elev');
    ghcnElevVar.putatt('description','elevation, in m, of nearest GHCN station');
    ghcnElevVar.putatt('units','m');
       
        % main variable(s)
            
    nc.putvar(timeVar);
    nc.putvar(latVar);
    nc.putvar(lonVar);
    for iv=1:nvars
%       fprintf('adding variable %s\n', outVarNames{iv});
        varVar = Variable(outVarNames{iv}, 'Dimensions',{time_dim, stn_dim},'Datatype','single','FillValue',NAFlag);
        varVar.putatt('units',outUnits{iv});
        varVar.putatt('long_name',long_names{iv});
        nc.putvar(varVar);
    end
    nc.putvar(stnIDVar);
    nc.putvar(elevVar);
    nc.putvar(tzVar);
    nc.putvar(stnNameVar);
    nc.putvar(st_dateVar);
    nc.putvar(end_dateVar);
    nc.putvar(pct_validVar);
    nc.putvar(is_ghcnVar);
    nc.putvar(nearest_ghcnVar);
    nc.putvar(ghcnDistVar);
    nc.putvar(ghcnAzVar);
    nc.putvar(ghcnElevVar);
    
        % write to the file and output all variables except the data variable.
    nc.writeschema(true,'Dimensions',true,'nofillmode',false);
    nc.writevars({'time'});
%     outvars = {'stnID','stnName','stnIDC','lat','lon','elevation','time_zone','nearest_ghcn', 'ghcn_distance','ghcn_az','ghcn_elev'};
%     nc.writevars(outvars);   
end
     
     
function write_stn_data(outVarName, vals, stn_ix, out_pos, nc, stn_tbl, st_yr, varnum, local_ghcn)
    stride = [1,1];
    try
        nc.writevar(outVarName, vals, [1,out_pos], stride);
    catch
        oops();
    end
    nc.writevar('stnID',char(stn_tbl.stnID(stn_ix))', [1,out_pos], stride);
    nc.writevar('stnName',char(stn_tbl.stnName(stn_ix))', [1,out_pos], stride);
    nc.writevar('lat',stn_tbl.lat(stn_ix), out_pos,1);
    nc.writevar('lon',stn_tbl.lon(stn_ix), out_pos,1);
    nc.writevar('elevation',stn_tbl.elev(stn_ix), out_pos,1);
    nc.writevar('time_zone',stn_tbl.time_zone(stn_ix), out_pos,1);
    nc.writevar('is_ghcn',stn_tbl.is_ghcn(stn_ix), out_pos, 1);
    if (~local_ghcn && ~isempty_s(stn_tbl.nearest_ghcn(stn_ix)))
        nc.writevar('nearest_ghcn', char(stn_tbl.nearest_ghcn(stn_ix))',[1,out_pos],stride);
        nc.writevar('ghcn_distance', stn_tbl.ghcn_dist(stn_ix), out_pos,1);
        nc.writevar('ghcn_az', stn_tbl.ghcn_az(stn_ix), out_pos,1);
        nc.writevar('ghcn_elev', stn_tbl.ghcn_elev(stn_ix), out_pos,1);
    end
    
    
    start_date = stn_tbl.startDate(stn_ix,varnum) - datenum(st_yr,1,1);
    end_date   = stn_tbl.endDate(stn_ix,varnum)   - datenum(st_yr,1,1);

    nc.writevar('start_date', start_date, [out_pos, varnum], [1,1]);
    nc.writevar('end_date',   end_date,   [out_pos, varnum], [1,1]);
    nc.writevar('pct_valid',  stn_tbl.pctValid(stn_ix,varnum), [out_pos, varnum], [1,1]);
        
end

function workTbl = find_closest_ghcns(workTbl, ghcnTbl)

    nstns = size(workTbl,1);
    if (~exist('ghcnTbl','var') || isempty(ghcnTbl))
        ghcnTbl = workTbl;
    end
    if (~any(strcmp('is_ghcn',fieldnames(ghcnTbl)))), return; end 
    ghcnTbl = ghcnTbl(logical(ghcnTbl.is_ghcn),:);
    if (size(ghcnTbl,1)==0), return; end 
    
    ngstns = size(ghcnTbl,1);
    radius = 6371;  % radius of earth in km
    km_per_degree = 2*pi*radius/360;
    
    for i=1:nstns
        if (~workTbl.is_ghcn(i))
                    % find the closest station in the inventory table
            stnlat = repmat(workTbl.lat(i), ngstns,1);
            stnlon = repmat(workTbl.lon(i), ngstns, 1);
            [dists,az] = distance(stnlat, stnlon, ghcnTbl.lat, ghcnTbl.lon);
            jx = find(dists == min(dists),1);
                % get distance and direction from inventory table
            workTbl.nearest_ghcn(i) = ghcnTbl.stnID(jx);
            workTbl.ghcn_dist(i) = dists(jx) * km_per_degree;
            workTbl.ghcn_az(i) = az(jx);
            workTbl.ghcn_elev(i) = ghcnTbl.elev(jx);
        end
    end
end

function rewrite_nearest_ghcn_data(nc, stn_tbl)
    nc.writevar('nearest_ghcn', char(stn_tbl.nearest_ghcn)');
    nc.writevar('ghcn_distance', stn_tbl.ghcn_dist);
    nc.writevar('ghcn_az', stn_tbl.ghcn_az);
    nc.writevar('ghcn_elev', stn_tbl.ghcn_elev);        
end

%function  workTbl = add_QC_columns(stnInfo, yrRange, RP)
function   workTbl = add_QC_columns(stnInfo,          RP)

    nstns = size(stnInfo,1);
    nv    = length(RP.varNames);
    
    workTbl = table( ...
                        strings(nstns,1), nan(nstns,1), nan(nstns,1), strings(nstns,1), nan(nstns,1), nan(nstns,1), nan(nstns,nv), nan(nstns,nv), nan(nstns,nv), false(nstns,1), repmat("unknown",nstns,1), nan(nstns,1), nan(nstns,1), nan(nstns,1), strings(nstns,1), ...
       'VariableNames',{'stnID',          'lat',        'lon',        'stnName',        'elev',       'time_zone',  'startDate',   'endDate',     'pctValid',    'is_ghcn',      'nearest_ghcn',            'ghcn_elev', 'ghcn_dist',       'ghcn_az',       'filename'});

    workTbl.stnID     = stnInfo.(RP.stnIDColumn);
    workTbl.lat       = stnInfo.(RP.latColumn);
    workTbl.lon       = stnInfo.(RP.lonColumn);
    workTbl.stnName   = stnInfo.(RP.stnNameColumn);
    workTbl.filename  = stnInfo.(RP.filenameColumn);
    
    if(ismember('is_ghcn',fieldnames(stnInfo))), workTbl.is_ghcn = logical(stnInfo.is_ghcn); end    
    if(ismember('ghcnID',fieldnames(stnInfo))), workTbl.nearest_ghcn = stnInfo.ghcnID; end  
    if(ismember('ghcn_dist',fieldnames(stnInfo))), workTbl.ghcn_dist = stnInfo.ghcn_dist; end
    if(ismember('ghcn_az',fieldnames(stnInfo))), workTbl.ghcn_az = stnInfo.ghcn_az; end
    if(ismember('ghcn_elev',fieldnames(stnInfo))), workTbl.ghcn_elev = stnInfo.ghcn_elev; end
    
    if (~isempty(RP.elevColumn) && ismember(RP.elevColumn, stnInfo.Properties.VariableNames))
        workTbl.elev  = stnInfo.(RP.elevColumn);
    else
%       workTbl.elev  = gdem2_elevations(workTbl.lat,workTbl.lon, true);
        workTbl.elev  = nan(nstns,1);        % We'll add elevations to the file later...
    end
    
    if (~isempty(RP.TZColumn) && ismember(RP.TZColumn, stnInfo.Properties.VariableNames))
        workTbl.time_zone = stnInfo.(RP.TZColumn);
    else
        
%       workTbl.time_zone = getTZ(workTbl.lat, workTbl.lon);
        workTbl.time_zone = nan(nstns,1);   % we'll add time zone info to the file later...
    end
    
            % hmmm. not quite sure what this was for.  I think it was for a dataset where I needed to get the data from
            % the nearest GHCN station to a given lat/lon.
%     minYrs = RP.min_yrs;
%     
%     if (RP.isghcn)
%         workTbl.nearest_ghcn = workTbl.stnID;
%         workTbl.ghcn_elev   = workTbl.elev;
%         workTbl.dist        = zeros(nstns,1);
%         workTbl.az          = zeros(nstns,1);
%     else
%             % this needs fixing, Ian!
%         nbrs=QC_all_stations_nearest_neighbors([workTbl.lat, workTbl.lon], "stations.Tmax.ncamerica.1850.2017.nc", 1, 100,[],yrRange, minYrs, false, false, false, false);
%     
%         if ((size(nbrs,1)) < size(workTbl,1))
%             fprintf("warning: can't find enough neighbors with at least %d years of data in %s", minYrs, "stations.Tmax.ncamerica.1850.2017.nc");
%             workTbl.nearest_ghcn = "none";
%             workTbl.ghcn_elev   = nan;
%             workTbl.dist        = nan;
%             workTbl.az          = nan;
% 
%         else
% 
%             missing_elevs = find(isnan(nbrs.elev));
%             if (~isempty(missing_elevs))
%                 nbrelevs=get_googlemaps_elevation(nbrs.lat(missing_elevs),nbrs.lon(missing_elevs));
%                 for i=1:length(missing_elevs)
%                     ie = missing_elevs(i);
%                     nbrs.elev(ie) = nbrelevs(i);
%                 end
%             end    
% 
%             workTbl.nearest_ghcn = nbrs.nbrID;
%             workTbl.ghcn_elev   = nbrs.elev;
%             workTbl.dist        = nbrs.dist;
%             workTbl.az          = nbrs.az;
%         end    
%     end
end

function [vals, min_date, max_date, nvalid, mpts, file_startdate, file_enddate] = read_csv(fnm, varColumn, timeColumn, tstamps, timestr, treatasempty, isUTC, toUTC, TZ)

% tstamps is matlab datenums.
% timeColumn is either name of single column with time column name, assumed here as "yyyymmddhh"
% Or can be names of year, month, day, hour, (and optionally minute and second columns)
    warning('off','MATLAB:table:ModifiedAndSavedVarnames');
    if (~isnumeric(treatasempty))
        tbl = readtable(fnm,'TreatAsEmpty',treatasempty);
        do_nans = false;
    else
        tbl = readtable(fnm);   % readtable(...) can't take a numeric value for the treatasempty, so we'll have to do this ourselves later.
        do_nans = true;
    end
        
    warning('on','MATLAB:table:ModifiedAndSavedVarnames');
    
    if (strcmpi(timestr,'hourly'))
        isHourly = true;
    else
        isHourly = false;
    end
        
    if (length(timeColumn)==1)     % for files w/ date/time column as yyyymmddhh
        v=tbl.(timeColumn{1});
        npts=length(v);
        yr=floor(v/1e6);
        mo=floor(mod(v,1000000)/10000);
        da=floor(mod(v,10000)/100);
        hh=floor(mod(v,100));
    else
        try
            yr=tbl.(timeColumn{1});
        catch
            vals = nan(length(tstamps),1);
            min_date = 0;
            max_date = 0;
            nvalid = 0;
            mpts = 0;
            return;
        end
        mo=tbl.(timeColumn{2});
        da=tbl.(timeColumn{3});
        if (length(timeColumn)>=4)
            hh=tbl.(timeColumn{4});
            if (length(timeColumn)>=5)
                mm=tbl.(timeColumn{5});
                if (length(timeColumn)==6)
                    ss=tbl.(timeColumn{6});
                end
            end
        end
        npts=length(yr);
    end
    
    bad=isnan(yr) | isnan(mo) | isnan(da);
    tbl=tbl(~bad,:);
    if (isempty(tbl))
        vals = nan(length(tstamps),1);
        min_date = 0;
        max_date = 0;
        nvalid = 0;
        mpts = 0;
        return;
    end
    if (~exist('hh','var')), hh=zeros(npts,1); end
    if (~exist('mm','var')), mm=zeros(npts,1); end
    if (~exist('ss','var')), ss=zeros(npts,1); end
        
    dnums=datenum([yr,mo,da,hh,mm,ss]);
    if (isHourly && ~isUTC && toUTC)
        dnums = dnums -  TZ/24;
    end
    dnums = dnums(~bad);
        
    file_startdate = dnums(1);
    file_enddate   = dnums(end);

    nout = length(tstamps);
    vals = nan(nout,1);
    
    if (isHourly)
        ix = 1+round((dnums-tstamps(1))*24);
    else
        ix = 1+round((dnums-tstamps(1)));
    end
    
%     try
%         if (any(ix < 1) || any(ix > nout))
%             fprintf('warning:  dates out of range for %s, %s to %s\n', fnm, datestr(dnums(1),'yyyy-mm-dd HH:MM:SS'), datestr(dnums(end),'yyyy-mm-dd HH:MM:SS'));
%         end
%     catch
%         oops();
%     end
    keepers = (ix > 0 & ix <= nout);        % flag the ones where the destination index is in range
    dest = ix(keepers);                     % dest is where they are going.
    vals(dest) = tbl.(varColumn)(keepers);    % move them to the right place.
    
    infs = isinf(vals);
    if (sum(infs) > 0)
        fprintf("%d infinite values found.  Setting to NAs\n", sum(infs));
        vals(infs) = nan;
    end
    mu = nanmedian(vals);
    delta = abs(1e-6*mu);
    nines = abs(vals-999.9) < delta;
    if (sum(nines) > 0)
        fprintf("%d 999.9's found.  Setting to NAs\n", sum(nines));
        vals(nines) = nan;        
    end
    if (do_nans)
        mu = nanmedian(vals);
        delta = abs(1e-6*mu);
        for i=1:length(treatasempty)
            mynans = abs(vals - treatasempty(i)) <= delta;
            fprintf("%d values matched %f\n", sum(mynans),treatasempty(i));
            vals(mynans) = nan;
        end
    end
    
    nvalid = sum(~isnan(vals));
    min_ix = find(~isnan(vals),1);
    max_ix = find(~isnan(vals),1,'last');
    min_date = tstamps(min_ix);
    max_date = tstamps(max_ix);
    mpts = max_ix - min_ix + 1;

end

function Parms = getParams(varNames, yrRange, varargin)

% yrRange, ncComments, toDegF, toMM, basedir, fnameCol

    p = inputParser;
    p.KeepUnmatched = true;        

    addRequired (p,'yrRange', @(x) isnumeric(x) && x(1)>=1800 && x(2)>=x(1) && x(2)<=2100);
    
    addParameter(p,'basedir','.');
    addParameter(p,'outdir','.');
    addParameter(p,'stnIDColumn','stnID');
    addParameter(p,'stnNameColumn','stnName');
    addParameter(p,'latColumn','lat');
    addParameter(p,'lonColumn','lon');
    addParameter(p,'elevColumn','elev');
    addParameter(p,'TZColumn','TZ');
    addParameter(p,'timeColumn',["year","month","day"]);
    addParameter(p,'ncComments', []);     % either a struct or name of a file to read keyword/pairs from for netCDF global comments.
    addParameter(p,'units',"");
    addParameter(p,'outUnits',"");
    addParameter(p,'long_names',"");
    addParameter(p,'outVarNames',"");
    addParameter(p,'ncName','out.nc');
    addParameter(p,'isHourly',false);
    addParameter(p,'isUTC',true);
    addParameter(p,'toUTC',true);
    addParameter(p,'source','');        % set to "merra" for Anne's mepdg stuff.
    addParameter(p,'min_yrs',25);
    
    addParameter(p,'filenameColumn',"");
    addParameter(p,'NA_in','NA');
    addParameter(p,'local_ghcn',false);
    addParameter(p,'ghcn_stninfo',[]);
    addParameter(p,'max_stations',[]);
    addParameter(p,'istep',1);
    
    addParameter(p,'toDegF',false);
    addParameter(p,'toDegC',false);
    addParameter(p,'toInches',false);
    addParameter(p,'toMM',false);

    parse(p, yrRange, varargin{:})
    
    Parms = p.Results;
    Parms.ncName = fullfile(Parms.outdir, Parms.ncName);
    Parms.varNames = string(varNames);
    if (strlength(Parms.outUnits)==0),   Parms.outUnits    = Parms.units;    end
    if (strlength(Parms.outVarNames)==0),Parms.outVarNames = Parms.varNames; end
    nv = length(Parms.varNames);
    nvo = length(Parms.outVarNames);
    if (nv ~= nvo),  error("error:  number outVarNames (%d) doesn't match number of varNames (%d)", nvo, nv); end
    
    unmatched=fieldnames(p.Unmatched);
    if (~isempty(unmatched))
        fprintf(2,"error:  unexpected input parameter(s):\n\t");
        fprintf(2,"%s ", unmatched);
        fprintf(2,'\n');
        error("please correct and rerun");
    end
    
        % read GHCN stninfo file if filename given.
    if (isempty(Parms.ghcn_stninfo))
        try
            load("/Volumes/lacie_1/data/obs/stations_netcdf/stations.Tmax.ncamerica.1850.2018.20_yrs_min.mat",'tbl');  
            Parms.ghcn_stninfo=tbl;
        catch
            fprintf("NOTE:  no ghcn info available\n");
        end
    elseif (ischar_s(Parms.ghcn_stninfo))
        [~,~,ext] = fileparts(Parms.ghcn_stninfo);
        if (strcmpi(ext,'.csv'))
            Parms.ghcn_stninfo = readtable(Parms.ghcn_stninfo);
        elseif (strcmpi(ext,'.mat'))
            load(Parms.ghcn_stninfo,'tbl');
            Parms.ghcn_stninfo = tbl;
        elseif (strcmpi(ext,'.nc'))
            Parms.ghcn_stninfo = QC_get_site_table(Parms.ghcn_stninfo);
        else
            error("error:  wrong file type for ghcn_stninfo: %s;  must be .csv or .nc", Parms.ghcn_stninfo);
        end
    end
    
            % some presets
    
    for iv = 1:nv
            % the following are for Anne's transportation study         NOTE:  this need updating for revised code!
            % I used temp_F, rh_F and prcp.in for Anne's merra (reanalysis) data as well, Only change needed is time
            % column for merra prcp.in.
        if (strcmpi(varNames(iv),"temp_F"))
            Parms.long_names(iv) = "hourly temperature";
            Parms.units(iv) = "DegF";
            Parms.outUnits(iv) = "DegF";
            Parms.varColumn(iv) = "temp_F";
            Parms.outVarNames(iv)="temp_F";
            Parms.timeColumn = "yyyymmddhh";
        elseif (strcmpi(varNames(iv),"rh_F"))
            Parms.long_names(iv) = "hourly relative humidity";
            Parms.units(iv) = "percent";
            Parms.outUnits(iv) = "percent";
            Parms.varColumn(iv) = "rh_F";
%           Parms.outVarNames(iv)="relativehumidity";       % why did I use this, Ian?
            Parms.outVarNames(iv)="rh_F";
            Parms.timeColumn = "yyyymmddhh";
        elseif (strcmpi(varNames(iv),"tasmax"))
            Parms.long_names(iv) = "daily max temp";
            Parms.units(iv) = "DegC";
            Parms.outUnits(iv) = "DegC";
            Parms.varColumn(iv) = "downscaled";
            Parms.timeColumn = ["year","month","day"];
        elseif (strcmpi(varNames(iv),"tasmin"))
            Parms.long_names(iv) = "daily min temp";
            Parms.units(iv) = "DegC";
            Parms.outUnits(iv) = "DegC";
            Parms.varColumn(iv) = "downscaled";
            Parms.timeColumn = ["year","month","day"];
        elseif (strcmpi(varNames(iv),"rhsmax"))
            Parms.long_names(iv) = "daily max relative humidity";
            Parms.units(iv) = "percent";
            Parms.outUnits(iv) = "percent";
            Parms.varColumn(iv) = "downscaled";
            Parms.timeColumn = ["year","month","day"];
        elseif (strcmpi(varNames(iv),"rhsmin"))
            Parms.long_names(iv) = "daily min relative humidity";
            Parms.units(iv) = "percent";
            Parms.outUnits(iv) = "percent";
            Parms.varColumn(iv) = "downscaled";
            Parms.timeColumn = ["year","month","day"];
        elseif (strcmpi(varNames(iv),"pr"))
            Parms.long_names(iv) = "daily precipitation";
            Parms.units(iv) = "mm";
            Parms.outUnits(iv) = "mm";
            Parms.varColumn(iv) = "downscaled";
            Parms.timeColumn = ["year","month","day"];
            
                % these are the standard QC output files from Ranjini's QC code
                % 
        elseif (strcmpi(varNames(iv),"Tmax"))
            Parms.long_names(iv) = "daily max temp";
            Parms.units(iv) = "DegC";
            Parms.outUnits(iv) = "DegC";
            Parms.varColumn(iv) = "Tmax";
            Parms.timeColumn = ["Year","Month","Day"];
            if (strlength(Parms.filenameColumn)==0)
                Parms.filenameColumn = "stnID";
            end
        elseif (strcmpi(varNames(iv),"Tmin"))
            Parms.long_names(iv) = "daily min temp";
            Parms.units(iv) = "DegC";
            Parms.outUnits(iv) = "DegC";
            Parms.varColumn(iv) = "Tmin";
            Parms.timeColumn = ["Year","Month","Day"];
            if (strlength(Parms.filenameColumn)==0)
                Parms.filenameColumn = "stnID";
            end
        elseif (strcmpi(varNames(iv),"Prec"))
            Parms.long_names(iv) = "daily precipitation";
            Parms.units(iv) = ".1mm";
            Parms.outUnits(iv) = "mm";
            Parms.outUnits(iv) = "mm";
            Parms.timeColumn = ["Year","Month","Day"];
            if (strlength(Parms.filenameColumn)==0)
                Parms.filenameColumn = "stnID";
            end

                % these are for the TX Lake Study files             NOTE:  these should be OK for revised code.

        elseif (strcmp(varNames(iv),'prcp.in'))   % prcp.in for merra data as well.
            Parms.units(iv)        = "inches";
            Parms.outUnits(iv)     = "mm";
            Parms.outVarNames(iv)  = "pr";
            Parms.long_names(iv)   = "total daily precipitation, mm";
            Parms.varColumn(iv)    = "prcp_in";     % changed '.' to '_' to match readtable's column-name conversion
            if (Parms.source == "merra")
                Parms.timeColumn  = "yyyymmddhh";
            else
                Parms.timeColumn  = ["year","month","day"];
            end
        elseif (strcmp(varNames(iv),'evap.in'))
            Parms.units(iv)        = "inches";
            Parms.outUnits(iv)     = "mm";
            Parms.varColumn(iv)    = "evap_in";     % changed '.' to '_' to match readtable's column-name conversion
            Parms.outVarNames(iv)  = "evap";
            Parms.long_names(iv)   = "evapotranspiration, mm";
            Parms.timeColumn  = ["year","month","day"];       
        elseif (strcmp(varNames(iv),'tmax.F'))
            Parms.units(iv)       = "DegF";
            Parms.outUnits(iv)    = "DegC";
            Parms.varColumn(iv)      = "tmax_F";     % changed '.' to '_' to match readtable's column-name conversion
            Parms.outVarNames(iv)  = "tasmax";
            Parms.long_names(iv)   = "daily max air temperature, deg C";
            Parms.timeColumn  = ["year","month","day"];       
        elseif (strcmp(varNames(iv),'tmin.F'))
            Parms.units(iv)       = "DegF";
            Parms.outUnits(iv)    = "DegC";
            Parms.varColumn(iv)      = "tmin_F";     % changed '.' to '_' to match readtable's column-name conversion
            Parms.outVarNames(iv)  = "tasmin";
            Parms.long_names(iv)   = "daily Min Temperature, deg C";
            Parms.timeColumn  = ["year","month","day"];       
        else    
            error("can't handle variable %s yet",varNames(iv));
        end  
        
        if (Parms.toDegF   && Parms.outUnits(iv)=="DegC"),   Parms.outUnits(iv) = "DegF";   end
        if (Parms.toDegC   && Parms.outUnits(iv)=="DegF"),   Parms.outUnits(iv) = "DegC";   end
        if (Parms.toInches && Parms.outUnits(iv)=="mm"),     Parms.outUnits(iv) = "inches"; end
        if (Parms.toMM     && Parms.outUnits(iv)=="inches"), Parms.outUnits(iv) = "mm";     end
        
            % For Anne's merra (mepdg) data:
        if (Parms.source == "merra")
            Parms.isUTC = false;
            Parms.toUTC = true;
        end
            
        
        if (strlength(Parms.filenameColumn)==0)
            Parms.filenameColumn = "filename";
        end
        
    end
    
    if (isrow(Parms.varNames)), Parms.varNames = Parms.varNames'; end
    if (isrow(Parms.outVarNames)), Parms.outVarNames = Parms.outVarNames'; end
        
    
    if (~isempty(Parms.ncComments))
        if (ischar_s(Parms.ncComments))
            if (exist(Parms.ncComments,'file')==2)
                Parms.ncComments = keyvalue_read(Parms.ncComments);
            else
                Parms.ncComments = struct('ncComments',Parms.ncComments);
            end
        elseif (~isstruct(Parms.ncComments))
            error('ncComments must be either a struct or a name of a file with keyword/value pairs for ncdf comments');
        end
    end
end
