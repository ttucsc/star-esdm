function ncOutName = QC_daily_to_hourly(ncHourly, ncMinvals, ncMaxvals, ncOutName, ncComments, ncAtts, toDegF,  useLocalTime, stns, yrRange_input, yrRange_output, outMinutes, inMinutes)
% Program to expand daily min/max into hourly, using an hourly temperature or humidity profile built from hourly observations.
%
%   Creates a 365-day set of daily smoothed profiles from the observations, giving an estimate of the mean temperature
%   or humidity for each time of day, each day of the year.  Then the max & min values are used to interpoloate that
%   profile to create a simulated hourly profile for the model data.  Temperature data puts the max in the aftermoon and
%   min near sunrise.  Humidity puts the max in early morning and min in late afternoon.
%
%   The program run_daily_to_hourly.m is a wrapper which will run this routine over multiple sets of data.
%
%   Inputs:
%       ncHourly            name of input file with hourly observation data
%       ncMinvals           names of input files containing daily min and max values to be expanded to hourly data
%       ncMaxvals
%       ncOutName           output netcdf file name
%       ncComments          Info to go into global attribute "comments"
%       ncAtts              Name Value array of string pairs of global attributes to add to file. Nx2 matrix of strings
%       toDegF              if true, convert temperatures from degC to degF
%       stns                (optional lists of stations to use.  Can be csv file or QC station table.  If file or table,
%                               must have column labeled "stnID".  Program will search for matching stnID, or if not
%                               found, matching station name in ncHourly, ncMinvals, ncMaxvals.
%       yrRange_input       (optional) years to use from hourly input  [hourly data's start date, hourly data's end date] 
%       yrRange_output      (optional) years to use from min, max input [1950,2100]
%       outMinutes          (optional) # of minutes between output samples [60 -- i.e., hourly]
%       inMinutes           (optional) # of minutes between input samples, [60]

    
    %-------------------variable and column identification.  Update this for your specific files.
    
    siteTblHourly  = QC_get_site_table(ncHourly);
    siteTblMinvals = QC_get_site_table(ncMinvals);
    siteTblMaxvals = QC_get_site_table(ncMaxvals);
    
    sources = [string(ncHourly), string(ncMinvals), string(ncMaxvals)];
    
    if (~exist('toDegF',      'var') || isempty(toDegF)),       toDegF       = false; end
    if (~exist('useLocalTime','var') || isempty(useLocalTime)), useLocalTime = false; end
    
    if (exist('stns','var') && ~isempty(stns))
        if (~istable(stns))
            if (size(stns,1) > 1)
                stns =  table(stns,'VariableNames',{'stnName'});
            else
                [~,~,ext] = fileparts(stns);
                if (strcmpi(ext,'.nc'))
                    stns = QC_get_site_table(stns);
                elseif (strcmpi(ext,'.txt') || strcmpi(ext,'.csv'))
                    stns = readtable(stns);
                    if (~any(strcmp(stns.Properties.VariableNames,'stnID')))
                        error('bad station list.  file must have column labeled "stnID"');
                    end
                end
            end
        end
        workTbl = QC_get_site_table(siteTblMinvals, "stnID", stns);
    else
        stns = siteTblHourly.stnID(:);
        workTbl = siteTblMinvals;
    end
    
%   nstns = size(stns,1);
    nstns = size(workTbl,1);
    
    if (~exist('yrRange_input','var') || isempty(yrRange_input))
        dv = datevec(siteTblHourly.startDate);
        if (length(unique(datenum(dv(:,1:3)))) ~= 1), error("error:  multiple start dates for hourly data;  please specify a start and end date"); end
        yr1 = dv(1,1);
        dv = datevec(siteTblHourly.endDate);
        if (length(unique(datenum(dv(:,1)))) ~= 1), error("error:  multiple end dates for hourly data;  please specify a start and end date"); end
        yr2 = dv(1,1);
        yrRange_input = [yr1,yr2];
    end
    if (~exist('yrRange_output','var') || isempty(yrRange_output))
        yrRange_output = [1950,2100];
    end
    
    if (~exist('outMinutes','var') || isempty(outMinutes))
        outMinutes = 60;
    end
    
    if (~exist('inMinutes','var') || isempty(inMinutes))
        inMinutes = 60;
    end
    
%     nterms_day = 6;
%     sigterm_day = 1;
%     nterms_yr = 6;
%     sigterm_yr = 3;
%     
    varName  = siteTblHourly.Properties.UserData.varName;
    units    = siteTblHourly.Properties.UserData.units;
    longName = siteTblHourly.Properties.UserData.longName;
    if (toDegF && strcmpi(units,'degC'))
        units='degF'; 
        longName = 'degrees Fahrenheit';
    end
            
    [nc, workTbl] = create_ncdf(ncOutName, workTbl, varName, units, longName, yrRange_output, ncComments, ncAtts, outMinutes, sources, useLocalTime);
    
    
    for istn=1:nstns
 %      stn = stns(istn,:);
        stnName = stns.stnName(istn);
        
                % kludge to work around Bristol, stn # 62, which only has valid data 1991-mid 2010. THISIS FOR THE OLD NTMWD DATA!  MEPDG just uses 24 stations
%         if (istn==62)
%             yrRange_input=[1991,2009];
%         end

                % look for stnName in table
        obs24  = QC_get_site_table(siteTblHourly, yrRange_input(1), yrRange_input(2),  "stnID", stnName, "removeLeaps",true, "searchType","any");
        if (isempty(obs24))    % try matching 1st 8 chars only
            obs24  = QC_get_site_table(siteTblHourly, yrRange_input(1), yrRange_input(2), "stnID", extractBefore(stnName,9),  "removeLeaps", true, "searchType","like_any");
            if isempty(obs24)
                fprintf(2, "warning:  cannot locate station %s in obs24 table", stnName);
                continue;
            else
                fprintf(2, "warning:  only partial match for station %s.  matched to %s\n", stnName, obs24.stnName(1));
            end
        end
                
        minday = QC_get_site_dadta(siteTblMinvals, yrRange_output(1), yrRange_output(2),  "stnID", stnName, "removeLeaps",true, "searchType","any");
        if (isempty(minday))
            minday  = QC_get_site_table(siteTblHourly,  yrRange_input(1), yrRange_input(2), "stnID", extractBefore(stnName,9), "removeLeaps", true, "searchType","like_any");
            if isempty(minday)      % try matching 1st 8 chars only
                fprintf(2, "warning:  cannot locate station %s in minday table", stnName);
                continue;
            else
                fprintf(2, "warning:  only partial match for station %s.  matched to %s\n", stnName, minday.stnName(1));
            end
        end
    
        maxday = QC_get_site_table(siteTblMaxvals, yrRange_output(1), yrRange_output(2), "stnID", stnName, "removeLeaps",true, "searchType","any");
        if (isempty(maxday))        % try matching 1st 8 chars only
            maxday  = QC_get_site_table(siteTblHourly,  yrRange_input(1), yrRange_input(2),  "stnID", extractBefore(stnName,9), "removeLeaps", true, "searchType","like_any");
            if isempty(maxday)
                fprintf(2, "warning:  cannot locate station %s in maxday table", stnName);
                continue;
            else
                fprintf(2, "warning:  only partial match for station %s.  matched to %s/n", stnName, maxday.stnName(1));
            end
        end

        obs_cal  =   obs24.Properties.UserData.calendar;
        work_cal = workTbl.Properties.UserData.calendar;

        isRH = strncmpi(obs24.Properties.UserData.varName, "rh",2) || strncmpi(obs24.Properties.UserData.varName, "relh",4);
        isUTC = obs24.Properties.UserData.isUTC;

        if (isnan(obs24.time_zone(1)))
            hrsEast = round(obs24.lon(1)/15);
        else
            hrsEast = obs24.time_zone(1);
        end
        
%       obsHourly = fixnans(obs24.data(1,:), nterms_day, sigterm_day);       % so few nans in Anne's data that we'll just interpolate over the missing data.
        obsHourly = obs24.data(1,:);
        if (isUTC)
            obsHourly = circshift(obsHourly, hrsEast);     % last few points are actually in previous day.
        end                                 % Circular shift adjusts the clock to local time.

%         minvals = fixnans(minday.data(1,:), nterms_yr, sigterm_yr);
%         maxvals = fixnans(maxday.data(1,:), nterms_yr, sigterm_yr);
        minvals = minday.data(1,:);
        maxvals = maxday.data(1,:);
        
                    % convert to degF from degC
        if(toDegF && strcmpi(obs24.Properties.UserData.units,'degC'))
            obsHourly = 32.0 + obsHourly * 9.0/5.0;
        end
        if(toDegF && strcmpi(minday.Properties.UserData.units,'degC'))
            minvals = 32.0 + minvals * 9.0/5.0;
        end
        if(toDegF && strcmpi(maxday.Properties.UserData.units,'degC'))
            maxvals = 32.0 + maxvals * 9.0/5.0;
        end
        
        hourly = daily_to_hourly(obsHourly, minvals, maxvals, isRH, outMinutes, inMinutes);
        
%         (need to update some fields in workTbl).
        
        if (~useLocalTime)
            hourly=shift_to_UTC(hourly, hrsEast, outMinutes);
        end
        
        workTbl = update_workTbl(workTbl, hourly, istn, outMinutes);      % update percent valid, and start, end dates.
        
        write_stn_data(varName, hourly, istn,   nc, workTbl);  

        try
            fprintf('%6d of %6d (%7.3f):  %-80s %-32s   hourly:  %s to %s  %7.3f %% valid, output %s to %s  %7.3f %% valid\n', istn, nstns, 100.0*istn/nstns, ncOutName, workTbl.stnName{istn}, ...
                    datestr_cal(  obs24.startDate(   1),  obs_cal, 'yyyy-mm-dd'), datestr_cal(  obs24.endDate(   1),  obs_cal, 'yyyy-mm-dd'),    obs24.pctValid(   1), ...
                    datestr_cal(workTbl.startDate(istn), work_cal, 'yyyy-mm-dd'), datestr_cal(workTbl.endDate(istn), work_cal, 'yyyy-mm-dd'), workTbl.pctValid(istn));
        catch me
            fprintf('oops!');
            msgtext=getReport(me);
            fprintf('%s\n', msgtext);
            fprintf('------\n');
        end            
    end
    
end

function workTbl = update_workTbl(workTbl, vals, istn, outMinutes)

    day1 = workTbl.Properties.UserData.day1;
    ix = find(~isnan(vals),1);
    if (isempty(ix))
        workTbl.startDate(istn) = nan;
        workTbl.pctValid(istn) = 0.0;
    else
        startdix = (ix-1)/24*60/outMinutes;
        workTbl.startDate(istn) = day1 + startdix;
        ix = find(~isnan(vals),1,'last');
        enddix = (ix-1)/24*60/outMinutes;
        workTbl.endDate(istn) = day1 + enddix;
        
        workTbl.pctValid(istn) = 100.0 * sum(~isnan(vals))/length(vals);
    end 
end

function [nc, workTbl] = create_ncdf(ncOutName, workTbl, varName, units,  longName, yrRange,        ncComments, ncAtts, outMinutes, sources, useLocalTime)
% function [nc, tstamps] = create_ncdf(outName, workTbl, yrRange, varName, longName, varUnits, NAFlag, ncComments)

%     varName  = workTbl.Properties.UserData.varName;
%     units    = workTbl.Properties.UserData.units;
%     longName = workTbl.Properties.UserData.longName;
%     NAFlag   = workTbl.Properties.UserData.NAFlag;    % this will be a float in most cases, but we want to use double(1e20)
    NAFlag   = 1e20;    % set to a double.  

    nc = ncdf('', 'Filename',ncOutName, 'Format','netcdf4');      % create an empty ncdf
    nc.putatt('title',sprintf('Generated Hourly %s Station Data', varName));
    nc.putatt('institution','Texas Tech Climate Science Center');
    nc.putatt('description',sprintf('TTU Generated Hourly station data for %s', varName));
    nc.putatt('date_range',sprintf('%s to %s', datestr(datenum(yrRange(1),1,1),'yyyy-mm-dd'), datestr(datenum(yrRange(2),12,31),'yyyy-mm-dd')));
    [~,uname] = getusername();
    nc.putatt('creation_date',datestr(now,'yyyy-mm-dd HH:MM:SS'));
    nc.putatt('creator',uname);
    nc.putatt('source_code',mfilename);
    nc.putatt('hourly_station_data',int32(1));
    nc.putatt('daily_station_data',int32(0));
    nc.putatt('data_variables',varName);
    nc.putatt('data_source',sprintf('%s ', sources{:}));


    for i=1:length(ncComments)
        if (i==1)
            comment_lbl="comment";
        else
            comment_lbl=sprintf("comment%d",i);
        end            
        nc.putatt(comment_lbl,ncComments{i});
    end
    for i=1:size(ncAtts,1)
        nc.putatt(ncAtts(i,1), ncAtts(i,2));
    end
    
    nstations = size(workTbl,1);
    
        % add time info
    st_yr = yrRange(1);
    end_yr = yrRange(2);
    nyrs = end_yr - st_yr + 1;
    ndays = 365*nyrs;
    day1 = datenum365([st_yr, 1, 1]);
    dayend = datenum365([end_yr, 12, 31]);
    nptsPerDay = 24*60/outMinutes;
    npts = ndays * nptsPerDay;
    timevals  = (0:(npts-1))/nptsPerDay;
    ntimes = length(timevals);
%     tstamps = day1 + timevals;
    
    
            % char(...) conveniently space-pads the strings to the same length.
    stnIDs = (char(workTbl.stnID))';      % transpose needed here to order the chars properly when writing out!
    stnIDLen = size(stnIDs,1);

    stnNames = (char(workTbl.stnName))';
    stnNameLen = size(stnNames,1);
    
    if (any(strcmp('nearestGHCN',fieldnames(workTbl))))
        is_ghcn=workTbl.is_GHCN;
        nearest_ghcn = (char(workTbl.nearestGHCN))';
        nearest_ghcnLen = size(nearest_ghcn,1);
        ghcn_distance = workTbl.GHCN_distance;
        ghcn_az = workTbl.GHCN_az;
        ghcn_elev = workTbl.GHCN_elev;
    elseif (any(strcmp('nearest_ghcn',fieldnames(workTbl))))
        is_ghcn=workTbl.is_ghcn;
        nearest_ghcn = (char(workTbl.nearest_ghcn))';
        nearest_ghcnLen = size(nearest_ghcn,1);
        ghcn_distance = workTbl.ghcn_distance;
        ghcn_az = workTbl.ghcn_az;
        if (any(strcmp('ghcn_elev',fieldnames(workTbl))))
            ghcn_elev = workTbl.ghcn_elev;
        else
            ghcn_elev = nan(nstations,1);
        end
    else
        nearest_ghcn=repmat("unknown",nstations,1);
        nearest_ghcnLen = 11;
        is_ghcn=false(nstations,1);
        ghcn_distance=nan(nstations,1);
        ghcn_az=nan(nstations,1);
        ghcn_elev=nan(nstations,1);
    end
        
    
        % define dimensions
    stn_dim = Dimension('stn_num',nstations);
    time_dim = Dimension('time',ntimes);
%     source_dim = Dimension('sources',3);
    idLen_dim = Dimension('stnIDLength',stnIDLen);
    nmLen_dim = Dimension('stnNameLength',stnNameLen);
    ghcnLen_dim = Dimension('ghcnIDLength',nearest_ghcnLen);
    
    nc.putdim(stn_dim);
    nc.putdim(time_dim);
%     nc.putdim(source_dim);
    nc.putdim(idLen_dim);
    nc.putdim(nmLen_dim);
    nc.putdim(ghcnLen_dim);
    
        % Create Variables
        
        % time
    if (useLocalTime), UTCflag = 'local';
    else,              UTCflag = 'UTC';  end
    timeUnits = sprintf('days since %s %s', datestr365(day1,'yyyy-mm-dd HH:MM:SS'), UTCflag);
    calendar='365-day';
    
    timeVar = Variable('time', timevals,'Dimensions',time_dim);
    timeVar.putatt('units',timeUnits);
    timeVar.putatt('calendar',calendar);
    timeVar.putatt('UTC_or_local',UTCflag);
    
        % latitude, longitude, elevation
        
    latVar = Variable('lat', workTbl.lat, 'Dimensions',stn_dim);
    latVar.putatt('long_name','latitude');
    latVar.putatt('standard_name','latitude');
    latVar.putatt('units','degrees north');
    lonVar = Variable('lon',workTbl.lon,'Dimensions',stn_dim);
    lonVar.putatt('long_name','longitude');
    lonVar.putatt('standard_name','longitude');
    lonVar.putatt('units','degrees east');
    elevVar = Variable('elevation', workTbl.elev,'Dimensions',stn_dim,'FillValue',NAFlag);
    elevVar.putatt('units','meters');
    elevVar.putatt('long_name','elevation, meters above sea level');
    tzVar = Variable('time_zone',workTbl.time_zone,'Dimensions',stn_dim);
    tzVar.putatt('units','hours');
    tzVar.putatt('longName','TZ hours east of GMT');
    
    st_dateVar = Variable('start_date','Dimensions',stn_dim,'Datatype','double','FillValue',NAFlag);
    st_dateVar.putatt('units',timeUnits);
    end_dateVar = Variable('end_date','Dimensions',stn_dim,'Datatype','double','FillValue',NAFlag);
    end_dateVar.putatt('units',timeUnits);
    pct_validVar = Variable('pct_valid','Dimensions',stn_dim,'Datatype','single','FillValue',NAFlag);
    pct_validVar.putatt('units','percent');
    
        % stations
        
    stnIDVar = Variable('stnID',stnIDs,'Dimensions',[idLen_dim, stn_dim]); 
    stnIDVar.putatt('long_name','Station ID');
    stnIDVar.putatt('description',sprintf('Station ID, stored as fixed length char strings %d chars long, padded with trailing spaces',stnIDLen));
    stnNameVar = Variable('stnName',stnNames,'Dimensions',[nmLen_dim, stn_dim]);
    stnNameVar.putatt('long_name','Station Name');
    stnNameVar.putatt('description',sprintf('Station Name, stored as fixed length char strings %d chars long, padded with trailing spaces',stnNameLen));    

    is_ghcnVar = Variable('is_ghcn',is_ghcn, 'Dimensions',stn_dim,'Datatype','uint8');
    is_ghcnVar.putatt('description','boolean, 1 if station is a GHCN station, 0 if not');
    
    nearest_ghcnVar = Variable('nearest_ghcn',nearest_ghcn, 'Dimensions', [ghcnLen_dim, stn_dim]);
    nearest_ghcnVar.putatt('long_name','nearest_GHCN_station');
    nearest_ghcnVar.putatt('description','Station ID of nearest TTU QC''d observation data from 2017 GHCN files');

    ghcnDistVar = Variable('ghcn_distance', ghcn_distance,'Dimensions',stn_dim);
    ghcnDistVar.putatt('long_name','nearest_GHCN_distance');
    ghcnDistVar.putatt('description','distance, in km, to nearest TTU QC''d observation data from 2017 GHCN files');
    ghcnDistVar.putatt('units','km');

    ghcnAzVar = Variable('ghcn_az', ghcn_az,'Dimensions',stn_dim);
    ghcnAzVar.putatt('long_name','ghcn_az');
    ghcnAzVar.putatt('description','direction (from true North) from station to nearest GHCN station');
    ghcnAzVar.putatt('units','degrees from true north');

    ghcnElevVar = Variable('ghcn_elev', ghcn_elev,'Dimensions',stn_dim);
    ghcnElevVar.putatt('long_name','ghcn_elev');
    ghcnElevVar.putatt('description','elevation, in m, of nearest GHCN station');
    ghcnElevVar.putatt('units','m');

        % main variable
        
    varVar = Variable(varName, 'Dimensions',[time_dim, stn_dim],'Datatype','single','FillValue',NAFlag);
    varVar.putatt('units',units);
    varVar.putatt('long_name',longName);
    
    nc.putvar(timeVar);
    nc.putvar(latVar);
    nc.putvar(lonVar);
    nc.putvar(varVar);
    nc.putvar(elevVar);
    nc.putvar(tzVar);
    nc.putvar(stnIDVar);
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
    
        % make sure all your variables are listed here, except the ones we'll write later:  varvar, start_date,
        % end_date, pct_valid, and time (already written)
    outvars = {'time','lat','lon','stnID','elevation','time_zone','stnName','start_date','end_date','pct_valid','is_ghcn', 'nearest_ghcn', 'ghcn_distance','ghcn_az','ghcn_elev'};
    nc.writevars(outvars);   
    
    
        % update the workTbl to reflect the new times, calendar, etc.
    workTbl.startDate = max(workTbl.startDate, day1);
    workTbl.endDate   = min(workTbl.endDate, dayend);
    
    workTbl.Properties.UserData.day1 = day1;
    workTbl.Properties.UserData.npts = npts;
    workTbl.Properties.UserData.dates = timevals + day1;
    workTbl.Properties.UserData.NAFlag = NAFlag;
    workTbl.Properties.UserData.varName = varName;
    workTbl.Properties.UserData.longName = longName;
    workTbl.Properties.UserData.units = units;
    workTbl.Properties.UserData.timeunits = timeUnits;
    workTbl.Properties.UserData.calendar = calendar;
        
end
     
%        write_stn_data(varName, hourly, istn,   nc, workTbl);     
function write_stn_data(varName, vals,   stn_ix, nc, workTbl)

    start=[1, stn_ix];
    stride = [1,1];
    nc.writevar(varName, vals, start, stride);
    
    start_date = workTbl.startDate(stn_ix) - workTbl.Properties.UserData.day1;
    end_date   = workTbl.endDate(stn_ix)   - workTbl.Properties.UserData.day1;
    nc.writevar('start_date', start_date, stn_ix, 1);
    nc.writevar('end_date',   end_date,   stn_ix, 1);
    nc.writevar('pct_valid',  workTbl.pctValid(stn_ix),  stn_ix, 1);    
    
end

function vals=shift_to_UTC(vals, hrsEast, outMinutes, useNAs)
%   Shifts the data series from local time to UTC time.
%   if useNAs is true, then inserts NAs at beginning or end.  Otherwise duplicates part of first or last day.
%   If east of Long 0, then fill in first few values and chop off last few.
%   If west of Long 0, then chop off first few values and fill in last few.
    
    if (~exist('useNAs','var')), useNAs=false; end
    if (isrow(vals))
        wasrow = true;
        vals=vals';
    else 
        wasrow = false;
    end
    
    nptsPerDay = round(24*outMinutes/60);
    shift = abs(round(hrsEast*outMinutes/60));
    len = length(vals);
    
    if (hrsEast < 0)
        if (useNAs)
            vals = [nan(shift,1); vals(1:end-shift)];     % chop off the last few points and append NAs.
        else
            ix=(nptsPerDay-shift+1):nptsPerDay;           % duplicate the end of the first day to start the series
            vals = [vals(ix); vals(1:(end-shift))];       % and chop off the last few points
        end
    else
        if (useNAs)
            vals = [vals((shift+1):end); nan(shift,1)];     % and chop off the first few points.
        else
            ix=len-nptsPerDay + (1:shift);                  % duplicate the beginning of the last day to end the series
            vals = [vals((1+shift):end); vals(ix)];         % and chop off the first few points.
        end
    end
    
    if (wasrow), vals = vals'; end

end
