function QCtbl_to_nc(tbl, ncName, do_matfile, ncComments)

    tblCols = tbl.Properties.VariableNames;
    calendar = tbl.Properties.UserData.calendar;
        % translate startDate, endDate columns back to dnums if they were text.
    if (any(strcmp(tblCols,"startDate")) && ~isnumeric(tbl.startDate))
        tbl.startDate = datenum_cal(tbl.startDate,calendar);
        tbl.endDate   = datenum_cal(tbl.endDate,  calendar);
    end



    tblCols = tbl.Properties.VariableNames;
    varName = tbl.Properties.UserData.varName;
    fprintf('creating netcdf %s\n', ncName);
    [nc, day1] = create_ncdf(ncName, tbl, ncComments);  
    fprintf('done creating netcdf\n');
    
    if (any(strcmp(tblCols, varName)))
        vname = varName;
    elseif (any(strcmp(tblCols, "data")))
        vname = "data";
    else
        tbl = QC_get_site_table(tbl,"loadData",true);
        vname = "data";
    end

    %---------------------Output filename

    nstns = length(tbl.stnID);

    for istn=1:nstns            % istep is for for debugging
        fprintf("stn %4d of %d (%6.1f%%): %s %s\n", istn, nstns, 100*istn/nstns, tbl.stnID{istn}, tbl.stnName{istn});
        write_stn_data(varName, vname, istn, nc, tbl, day1);
    end

    if (do_matfile)
        [d,fbase,fext] = fileparts(ncName);
        matname = fullfile(d,fbase,".mat");
        
    end

end

function [nc, day1] = create_ncdf(outName, tbl, ncComments)% , yrRange, min_yrs, outVarNames, long_names, outUnits, NAFlag, ncComments,timestr, local_ghcn)

    calendar  = tbl.Properties.UserData.calendar;   
    timeunits = tbl.Properties.UserData.timeunits;
    startvec = nc_parse_date_str(timeunits); 
    day1 = datenum_cal(startvec, calendar);
    tstamps = tbl.Properties.UserData.dates;
    timevals = tstamps - day1;
    ntimes = length(timevals);
    
    NAFlag    = tbl.Properties.UserData.NAFlag;
    units     = tbl.Properties.UserData.units;
    long_name = tbl.Properties.UserData.longName;
    
    if (tbl.Properties.UserData.isUTC == 1)
        UTCstr = 'UTC';
    else
        UTCstr = 'local';
    end
    
    dvecs = datevec_cal(tstamps, calendar);%    yrRange = [min(dvecs(:,1)), max(dvecs(:,1))];
    varName = tbl.Properties.UserData.varName;
    [~,uname] = getusername();
    srcfilename = tbl.Properties.UserData.ncName;
    
    globalAttributes = tbl.Properties.UserData.globalAttributes;
    
    globalAttributes.date_range = sprintf("%s to %s", datestr(dvecs(1,:)), datestr(dvecs(end)));
    globalAttributes.creation_date = datestr(now,'yyyy-mm-dd HH:MM:SS');
    globalAttributes.creator = uname;
    globalAttributes.source_code = mfilename;
    globalAttributes.data_variables = varName;
    globalAttributes.original_source = tbl.Properties.UserData.ncName;
    
    nc = ncdf('', 'Filename',outName, 'Format','netcdf4','create_ok',true);      % create an empty ncdf
    fields = fieldnames(globalAttributes);
    for i=1:length(fields)
        nc.putatt(fields{i}, globalAttributes.(fields{i}));
    end

    if (isstruct(ncComments))
        fn=fieldnames(ncComments);
        for i=1:length(fn)
            nc.putatt(fn{i}, ncComments.(fn{i}));
        end
    else
        ncComments = string(ncComments);
        if (length(ncComments)==1)
            nc.putatt('comment', ncComments{1});
        else
            for i=1:length(ncComments)
                comment_lbl=sprintf("comment%d",i);
                nc.putatt(comment_lbl,ncComments{i});
            end
        end
    end    
    
%     nstations = size(tbl,1);

%     ndays = daysdif(datenum([yrRange(1),1,1]),datenum([yrRange(2)+1,1,1]));
%     day1 = datenum([yrRange(1), 1, 1]);
%     if (strcmpi(timestr,"hourly"))
%         timevals    = (0:(24*ndays-1))/24;
%     else
%         timevals    = (0:(ndays - 1));
%     end
%     ntimes = length(timevals);
%     tstamps = day1 + timevals;
    
    
            % char(...) conveniently space-pads the strings to the same length.
%   stnIDs = (char(tbl.stnID))';      % transpose needed here to order the chars properly when writing out!
    stnIDLen = max(strlength(tbl.stnID));

%   stnNames = (char(tbl.stnName))';
    stnNameLen = max(strlength(tbl.stnName));
    
%   nearest_ghcn = (char(tbl.nearest_ghcn))';
    nearest_ghcnLen = max(11, max(strlength(tbl.nearest_ghcn))); % size(nearest_ghcn,1);

    nvars=1; % length(outVarNames);
    coutVarNames = (char(varName))'; % (char(outVarNames))';
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
    timeVar = Variable('time', timevals,'Dimensions',time_dim);
    timeVar.putatt('units',timeunits);
    timeVar.putatt('calendar',calendar);
    timeVar.putatt('UTC_or_local',UTCstr);
    
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
    st_dateVar.putatt('units',timeunits);
    st_dateVar.putatt('calendar',calendar);
    end_dateVar = Variable('end_date','Dimensions',{stn_dim, nvar_dim},'Datatype','double','FillValue',NAFlag);
    end_dateVar.putatt('units',timeunits);
    end_dateVar.putatt('calendar',calendar);
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
    nearest_ghcnVar.putatt('description',sprintf('Station ID of nearest TTU QC''d observation data from %s', srcfilename));

       
    
    ghcnDistVar = Variable('ghcn_distance', 'Dimensions',stn_dim,'Datatype','single','FillValue',single(NAFlag));
    ghcnDistVar.putatt('long_name','nearest_ghcn_distance');

    ghcnDistVar.putatt('description',sprintf('distance, in km, to nearest TTU QC''d observation data from %s', srcfilename));

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

%       fprintf('adding variable %s\n', outVarNames{iv});
    varVar = Variable(varName, 'Dimensions',{time_dim, stn_dim},'Datatype','single','FillValue',NAFlag);
    varVar.putatt('units',units);
    varVar.putatt('long_name',long_name);
    nc.putvar(varVar);

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
          
function write_stn_data(outVarName, vname, stn_ix, nc, stn_tbl, day1)
    stride = [1,1];
    vals = stn_tbl.(vname)(stn_ix,:)';
    try
        nc.writevar(outVarName, vals, [1,stn_ix], stride);
    catch
        oops();
    end
    nc.writevar('stnID',char(stn_tbl.stnID(stn_ix))', [1,stn_ix], stride);
    nc.writevar('stnName',char(stn_tbl.stnName(stn_ix))', [1,stn_ix], stride);
    nc.writevar('lat',stn_tbl.lat(stn_ix), stn_ix,1);
    nc.writevar('lon',stn_tbl.lon(stn_ix), stn_ix,1);
    nc.writevar('elevation',stn_tbl.elev(stn_ix), stn_ix,1);
    nc.writevar('time_zone',stn_tbl.time_zone(stn_ix), stn_ix,1);
    nc.writevar('is_ghcn',stn_tbl.is_ghcn(stn_ix), stn_ix, 1);
    nc.writevar('nearest_ghcn', char(stn_tbl.nearest_ghcn(stn_ix))',[1,stn_ix],stride);
    nc.writevar('ghcn_distance', stn_tbl.ghcn_distance(stn_ix), stn_ix,1);
    nc.writevar('ghcn_az', stn_tbl.ghcn_az(stn_ix), stn_ix,1);
    nc.writevar('ghcn_elev', stn_tbl.ghcn_elev(stn_ix), stn_ix,1);
    
    
    start_date = stn_tbl.startDate(stn_ix,1) - day1;
    end_date   = stn_tbl.endDate(stn_ix,1)   - day1;

    nc.writevar('start_date', start_date, [stn_ix, 1], [1,1]);
    nc.writevar('end_date',   end_date,   [stn_ix, 1], [1,1]);
    nc.writevar('pct_valid',  stn_tbl.pctValid(stn_ix,1), [stn_ix, 1], [1,1]);
        
end
