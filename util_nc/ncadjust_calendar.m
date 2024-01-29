function [ out_start_vec, out_days_since, out_data, keepers, inserts] = ncadjust_calendar( start_vec, days_since, data, calendar, out_calendar, do_random, out_start_vec, dimid  )
    % Recodes netcdf time info and data for a new calendar type.
    %
    %   Possible calendars are standard 365.25-day (various names), 365-day (no leap), and 360-day (12 30-day months)
    %
    %   Inputs:
    %       start_vec       datevec (or just start year) from time units, "days_since..."
    %       days_since      vector of days_since from netcdf file's time variable
    %       data            1-d, 3-d data, (lon, lat, time) corresponding to days_since
    %       calendar        calendar string
    %       out_calendar    desired output calendar
    %       do_random       optional boolean.  If true, then inserts or deletes occur randomly
    %       out_start_vec   optional datevec (or just start year) to adjust all days_since to
    %       dimid           optional dimension to operate on.  If not present, operates on the longest dimension.
    %   Outputs
    %       out_start_vec   output start vec (or just start year is so specified)
    %       out_days_since  output days_since array
    %       out_data        output data w/ NAs inserted or with leap & possibly 5 extra days deleted
    %       keepers         logical array, true where input date retained, false where removed.
    
    %  note:  do_random probably not right yet, ian!
    
    % fprintf("ncadjust_calendar: calendar %s out_calendar: %s do_random: %d\n", calendar, out_calendar, do_random);
    
            % save the current state of the random number generator, and switch to multFibonacci RNG.    
    if (isnumeric(do_random))
        orig_rng = rng(do_random, 'multFibonacci');
    else
        orig_rng = rng(0,'multFibonacci');
    end
    
    
    incal_len  = calendar_length(calendar);
    outcal_len = calendar_length(out_calendar);
    npts = length(days_since);

        % svec & ovec are datevec's from start_vec and out_start_vec.
    if (length(start_vec) == 1)
        svec = [start_vec,1,1];
    else
        svec = start_vec;
        if (incal_len == 360 && (svec(2)>1 || svec(3)>1))
            sv1 = svec(1) + ((svec(2)-1)*30 + svec(3))/360;
            svec = datevec(round(datenum(sv1,0,0)));            % svec as a standard calendar date
        end
    end
    if (~exist('out_start_vec','var') || isempty_s(out_start_vec))
        out_start_vec = start_vec;
    end
    if (length(out_start_vec) == 1)
        ovec = [out_start_vec,1,1];
    else
        ovec = out_start_vec;
        if (outcal_len == 360 && (ovec(2)>1 || ovec(3)>1))
            ov1 = ovec(1) + ((ovec(2)-1)*30 + ovec(3))/360;
            ovec = datevec(round(datenum(ov1,0,0)));            % ovec as a standard calendar date
        end
    end
    if (~exist('dimid','var')), dimid=[]; end
    
    svec=svec(1:3); 
    ovec=ovec(1:3); 
    

    if (incal_len == outcal_len)        % same calendar.  output is input, but we might adjust days_since later.
        out_start_vec = start_vec;
        out_days_since = days_since;
        out_data = data;
        keepers = true(size(days_since));
        inserts = false(size(keepers));
    else

                    % get list of years in the data and a boolean array, true if day is leapday.
                    % also calculates # of leap days in the period between start_vec and 1st date and # of days to delete
                    % if going from 365 to 360-day calendar. 
      
        [npreadds_360, leapdays, npreleaps, npredels_360, ndels_360] = find_yrs(svec, days_since, calendar);        

        if (incal_len == 365.25)        % std. calendar to 365- or 360-day 
                    % first:  remove leap day
            if (do_random)
                ndels = sum(leapdays);
                keepers = find_keepers(true(npts,1), ndels, do_random);
            else
                keepers = ~leapdays;
            end
                    % now remove 5 days from each year if going to 360-day calendar.
            if (outcal_len == 360)
                keepers = find_keepers(keepers,ndels_360,do_random);
                day1 = days_since(1) - npreleaps - npredels_360;
            else
                day1 = days_since(1) - npreleaps;            
            end
            [out_data,nout] = keep_data(data, keepers, dimid);
            out_days_since = day1 + (0:(nout-1));

        elseif (incal_len == 365)

            if (outcal_len == 360)          % reducing to 360-day.  figure out which days to keep
                keepers = find_keepers(true(npts,1),ndels_360,do_random);
                [out_data, nout] = keep_data(data, keepers, dimid); 
                day1 = days_since(1) - npredels_360;
                out_days_since = day1 + (0:(nout-1));
                inserts = false(size(out_days_since));

            else                            % going to 365.25 day.  insert leaps.
                days_since = to_standard(svec, calendar, days_since, npreleaps, npreadds_360, do_random);
                [out_data, nout, inserts] = insert_data(data, days_since);
                day1 = days_since(1);
                out_days_since = day1 + (0:(nout-1));
            end

        else                                % 360-day to 365- or 365.25-day.
                                % go to 365-day .
            if (outcal_len == 365)
                days_since = to_365(days_since, calendar, npreadds_360, do_random);
                [out_data, nout, inserts] = insert_data(data, days_since, dimid);
                day1 = days_since(1) + npreadds_360;
                out_days_since = day1 + (0:(nout-1));
            else                % go to standard calendar
                days_since = to_standard(svec, calendar, days_since, npreleaps, npreadds_360, do_random);
                [out_data, nout, inserts] = insert_data(data, days_since, dimid);
                day1 = days_since(1);
                out_days_since = day1 + (0:(nout-1));
            end 
        end
    end


            % adjust days_since to out_start_vec if not the same.
             
    if (~isequal(svec, ovec))            
        if (outcal_len == 360)
            offset = days360(datenum(svec), datenum(ovec));
        elseif (outcal_len == 365)
            offset = days365(datenum(svec), datenum(ovec));
        else
            offset = datenum(ovec) - datenum(svec);
        end
        out_days_since = out_days_since - offset;           
    end
    
    rng(orig_rng);
end

function [data, npts] = keep_data(data, keepers, dimid)
    
    sz = size(data);
    nd=length(sz);
    if (isempty(dimid))
        dimid=find(sz==max(sz));
    end
    
    if (nd==2)
        if (dimid==1)
            data = data(keepers,:);
        else
            data = data(:,keepers);
        end
    elseif (nd==3)
        if (dimid==1)
            data=data(keepers,:,:);
        elseif (dimid == 2)
            data=data(:,keepers,:);
        else
            data=data(:,:,keepers);
        end
    else
        throw(MException('ICSF:BAD_DIMS',sprintf('error: %s:  too many dimensions.  dim(data)=%d',mfilename,nd)));
    end
    
    npts=sum(keepers);
end

function [outdata, npts, inserts] = insert_data(data, dix, dimid)
% creates a 2D or 3D matrix of NAs, and inserts the data into the matrix at locations identified by vector dix (day-index)
% locations in dix are relative;  min(dix) is subtracted from dix before inserting, 
% so that dix can be be either datenums or days_since.
%
%   Inputs:
%       data        input data to have nan's inserted
%       dix         (relative) day-indexes of where data should be copied
%                       length(dix) must equal length(data), along data's the time dimension (usually longest dim)
%       dimid       dimension where data should be inserted (usually the time axis dimension)
%   Returns:
%       outdata     data, with nan's inserted at locations not included in dix
%       inserts     logigal array, size (1,nout), true where data inserted, false where original data copied in.
    
    sz = size(data);
    nd=length(sz);
    if (isempty(dimid))
        dimid=find(sz==max(sz));
    end
    
    npts=max(dix)-min(dix)+1;
    sz(dimid)=npts;               % biggest dimension.
    outdata=nan(sz);
    
    ix=dix-min(dix) + 1;
    
    if (isempty(data))
        outdata = [];
    else
        if (nd==2)
            if (dimid==1)
                outdata(ix,:) = data;
            else
                outdata(:,ix) = data;
            end
        elseif (nd==3)
            if (dimid==1)
                outdata(ix,:,:) = data;
            elseif (dimid == 2)
                outdata(:,ix,:) = data;
            else
                outdata(:,:,ix) = data;
            end
        else
            throw(MException('ICSF:BAD_DIMS',sprintf('error: %s:  too many dimensions.  dim(data)=%d',mfilename,nd)));
        end
    end
    inserts = true(1,npts);
    inserts(ix) = false;
end

function keepers = find_keepers(keepers, ndels, do_random)        
    
    npredels = sum(~keepers);
    npts = length(keepers);
    step = npts/ndels;
    istep = floor(step);
    hstep = ceil(istep/2);
%    fprintf("-----find_keepers:  first 50 deletes are:");
    for i=0:(ndels-1)
        i1=floor(i*step)+1;
        if (do_random)
            i2=min(i1+istep-1, npts);
            try
                ix = i1  + find_random(keepers(i1:i2));
            catch me
                fprintf(2, 'oops! %s\n', me.message);
            end
        else
            ix = i1 + hstep;
        end
        keepers(ix)=false;
%        if (i < 50), fprintf("%5d ", ix);end
    end
%    fprintf(" -----\n");
    odels = sum(~keepers);
    if (odels ~= ndels + npredels)
        throw(MException('ICSF:OOPS',sprintf('Bug here...deleting %d instead of %d',odels, ndels+npredels)));
    end

end

function ix = find_random(flags)

    ixflags = find(flags);      % list of days we can select from  (flag is true if selectable, false if not.)
    if (isempty(ixflags))
        throw(MException('ICSF:NODATA','error:  ncadjust_calendar_ic:  nothing left to delete'));
    end
    nix = length(ixflags);      
    ixx = randi(nix);           % choose one of the days
    ix = ixflags(ixx);          % return index of selected day
end

function days_since = to_standard(svec, calendar, days_since, npreleaps, npreadds_360, do_random)
        % adjusts days_since to a standard calendar.
        
    cal_len = calendar_length(calendar);
    if (cal_len == 365.25)
        return;
    end
    
            % if 360-day, then make it 365 day 1st.
    if (cal_len == 360)
        days_since = to_365(days_since, calendar, npreadds_360, do_random);
    end
    
        % daynoleap2datenum needs days since 1/1/pivot_year..
        
    day_off=datenum(svec)-datenum(svec(1),1,1);
    if (isleap(svec(1) && svec(2)>2))
        day_off = day_off-1;
    end
    mydays_since = days_since + day_off;             % now relative to start of year of svec
    pivot_yr = svec(1);
    
    npts = length(mydays_since);
    dn1 = daynoleap2datenum(min(mydays_since),pivot_yr);
    dn2 = daynoleap2datenum(max(mydays_since),pivot_yr);

    if (~do_random)
        leapers = find(isleapday(dn1:dn2));
        nleaps = length(leapers);
    else
%        fprintf("-----to_standard(%d)  first 50 leaps are ",calendar_length(calendar));
        nleaps = sum(isleapday(dn1:dn2));        
        leapers = zeros(nleaps,1);
        step = npts/nleaps;
        istep = floor(step);
        for i=1:nleaps
            jx = randi(istep);
%            if (i<50), fprintf("%5d ", jx); end 
            leapers(i) = i*istep - jx;
        end
%        fprintf(" -----\n");
    end
    for i=nleaps:-1:1   % start at the end and work backwards inserting days.
        mydays_since(leapers(i):end) = mydays_since(leapers(i):end)+1;
    end
    days_since = mydays_since - day_off + npreleaps;
end

function outdays_since = to_365(days_since, calendar, npreadds_360, do_random)
%   returns days_since svec for a 365-day calendar for 360-day days_since.
    cal_len = calendar_length(calendar);
    if (cal_len ~= 360)
        throw(MException('ICSF:BAD_CALENDAR', sprintf('error:  to_365():  calendar %s not 360-day', calendar)));
    end
    
    dmin = min(days_since);
    dmax = max(days_since);
    drange = dmax - dmin + 1;
    nadds = round(5*(dmax-dmin+1)/360);
    step = drange/nadds;
    istep = floor(step)-1;
    hstep = floor(step/2);
    outdays_since = days_since + npreadds_360;
%    fprintf("-----to_365:  first 50 inserts are: ");
    for i=1:nadds
        if (~do_random)
            iday = dmin + floor(i*step - hstep);
        else
            iday = dmin + floor(i*step  - randi(istep));
        end
        ix = find(days_since>=iday,1);
%        if (i <= 50), fprintf("%5d ", ix); end
        outdays_since(ix:end)=outdays_since(ix:end)+1;
    end
%    fprintf(" -----\n");
    
end

function [npreadds_360, leapdays, npreleaps, npredels_360, ndels_360] = find_yrs(svec, days_since, calendar)

    cal_len = calendar_length(calendar);
    
    if (cal_len == 365.25)      % standard calendar
        npreadds_360 = [];
        
        dnums = datenum(svec) + days_since;
        leapdays = isleapday(dnums);                            % boolean vector, true if day is a leapday
        
        prenums = datenum(svec) + 0:(days_since(1)-1);
        npreleaps = sum(isleapday(prenums));                     % # of leap days years between svec & day1
                
        npre_yrs = (days_since(1)-npreleaps)/365;
        day1 = floor(datenum(svec+[npre_yrs,0,0]));
        npreleaps = sum(isleapday(datenum(svec):(day1-1)));      % # of leap days needed to add at start to go from 365 to standard calendar
        
        npredels_360 = round(5*npre_yrs);                        % # of days to delete when going to 360-day year from a 365-day year.
                
        nyrs = (length(days_since) - sum(leapdays))/365;
        ndels_360 = round(nyrs*5);
        
    elseif (cal_len == 365)     % 365-day calendar
        npreadds_360 = [];
        
        leapdays = false(size(days_since));
        npreleaps = 0;
        
        npre_yrs = days_since(1)/365;
        npredels_360 = round(5*npre_yrs);                        % # of days to delete when going to 360-day year from a 365-day year.
        
        ndels_360 = round(5*length(days_since)/365);
        
    else                        % 360-day calendar
        npre_yrs = days_since(1)/360;
        day1 = floor(datenum(svec+[npre_yrs,0,0]));
        npreleaps = sum(isleapday(datenum(svec):(day1-1)));     % # of leap days needed to add at start to go from 365 to standard calendar
        npreadds_360 = round(5*npre_yrs);                       % # of days to add to days_since to go from 360 to 365-day calendar.
        
        leapdays = false(size(days_since));
        
        npredels_360 = 0;
        
        ndels_360 = 0;
    end
end
