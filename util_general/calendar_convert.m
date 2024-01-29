function [out_dnums, inserted, deleted, to_ix, keepers] = calendar_convert(indates, incal, outcal, varargin)
%
%   Inputs:             
%       indates         input dates, in incal space.  Can be datevecs, or datenums or datestr's (such as "2020/04/30")
%                           if (only 2 dates, then dates are treated as start date and end date. 
%       incal           input  calendar.  string e.g. "360-day", "365-day", "standard", "julian", etc.
%       outcal          output calendar  (like incal)
%       varagin         optional keyword/value pairs:  do_random (true/false or rng seed) & start_date (if dates are ncdf-style days_since w/ start_date) 
%   Outputs:
%       out_dnums       datenums (in outcal space)
%       inserted        indexes of inserted dates (if outcal > incal)
%       deleted         indexes of deleted dates (if outcal < incal).   can use:  my_outdata = my_indata; my_outdata(deleted)=[];
%       to_ix           index of where to put data, if calendar_length(outcal) >= calendar_length(incal).    (empty if not)                  
%                           if calendar_length(outcal) >= calendar_length(incal), then you can use:  my_outdata(to_ix) = my_indata;   outdata(inserted) = nan; 
%       keepers         logical array of values to keep, if calendar_length(outcal) <= calendar_length(incal).    (empty if not)       
%                           if calendar_length(outcal) <= calendar_length(incal), then you can use:  my_outdata = my_indata(keepers);
%                           ( keepers(deleted) are false;  keepers(~deleted) are true )
%
% for 360-day<->365-day:  insert or delete 2/4, 4/17, 6/28, 9/8 & 11/19
%                                   doy:    35,  107,  179, 251 & 323.

% for standard, add or delete 2/29 in leap years.

% do_random:  if true, days insertion or deletion dates are perturbed by adding gaussian noise with mean of 0 &stdev of 7  (limited to +/- 35) 

%   Ian - this needs optional start_date, to handle indates as days_since start_date.
%           make do_random and start_date keyword/value pairs.
%         ALSO:  need to test do_random.

%   [datenums_out, inserted, deleted] = datenum_cal(datenums_in, incal, outcal)

    % convert indates to datenums.
    
    if (nargin < 2), error("error:  usage:  [out_dnums, inserted, deleted, from_ix, keepers] = calendar_convert(indates, incal, outcal [+kwd/val for  do_random, start_date])"); end
    if (~exist("outcal",    "var") || isempty(outcal)), outcal = incal; end    % this lets us easily generate a set of dates from start and end dates (as strings, datevecs or datenums)

 
    p = inputParser;
    addParameter(p,"do_random",false, @(s) islogical(s) || (isnumeric(s) && mod(s,1)==0));     
    addParameter(p,"start_date",[],   @(s) isempty(s) || isstring(s) || (isnumeric(s) && any(numel(s)==[1,3,6])));

    parse(p, varargin{:});
    do_random = p.Results.do_random;
    start_date = p.Results.start_date;     
    
        % handle case where dates are  from netcdf, with "days since (start_date)".  start_date can be datevec or string
        % date, with or without "days since " at start.
    if (~isempty(start_date))
        if (isstring(start_date) || ischar(start_date))
            if (strncmpi(start_date,"days since ",11))
                start_date = extractAfter(start_date,11);
            end
        end
        start_dnum = datenum_cal(start_date, incal);
        indates = start_dnum + indates;
    end
    
   if (isrow(indates) && ~any(size(indates,2)==[3,6]))
        indates = indates';
    end
    if (any(size(indates,2)==[3,6]) || isstring(indates) || ischar(indates) || (iscell(indates) && ischar(indates{1})))
        in_dnums = datenum_cal(indates, incal);
    else
        in_dnums = indates;
    end
    if (length(in_dnums)==2)
        in_dnums = (in_dnums(1):in_dnums(2))';
    end
    
    incal_len  = calendar_length(incal);
    outcal_len = calendar_length(outcal);

    if (incal_len == outcal_len)
        out_dnums = in_dnums;
        inserted = [];
        deleted  = [];
        to_ix = (1:length(in_dnums))';
        if (nargout > 4), keepers = true(size(to_ix)); end
        return;
    end
    
    npts_in = size(in_dnums,1);
    
        % remove and store any fractional days.
        % calc average fractional day (for inserted days)        
    
    dayfracs = mod(in_dnums,1.0);
    if (any(dayfracs ~= 0))
        has_dayfracs = true;
        in_dnums = floor(in_dnums);
        if (outcal_len > incal_len), mean_dayfrac = mean(dayfracs); end
    else
        has_dayfracs = false;
    end
    
        % get starting and ending datevecs, and create output datenums vector.
    dvecs_in = datevec_cal(in_dnums, incal);    
    dv1  = dvecs_in(1,:);
    dv2i = dvecs_in(end,:);
    dv2o = dv2i;
            % if ending on last day of a year, fix day-of-month for last day (dv2[i|o]) if necessary
    if (incal_len==360)
        if (dv2i(2)==12 && dv2i(3)==30)
            dv2o(3)=31;                     % end output on Dec 31st if going from 360 to longer calendar.
        end
    elseif (outcal_len==360)
        if (dv2i(2)==12 && dv2i(3)==31)
            dv2o(3)=30;                     % end output on Dec 30th if going to a 360-day calendar. 
        end
    end
    
    % find missing dates in input.

    if (length(in_dnums) < (in_dnums(end)-in_dnums(1))+1 )
        in_dnums_all = in_dnums(1):in_dnums(end);
        missing_dnums_in = setdiff(in_dnums_all,in_dnums);
    else
        missing_dnums_in = [];
    end
        
        % complete list of output datenums.  (will delete any missing input dates later)
    out_dnums_all = (datenum_cal(dv1, outcal):datenum_cal(dv2o, outcal))';
    out_dnums = out_dnums_all;
    npts_out  = length(out_dnums_all);
    dvecs_out = datevec_cal(out_dnums_all, outcal);
        
        % convert to out_dnum values.
    if (outcal_len < incal_len)
        if (outcal_len == 360)
                    % remove any dates matching:  2/4/????, 2/29/????, 4/17/????, 6/28/????, 9/8/???? & 11/19/????, 
            deleted = find((dvecs_in(:,2)==2 & dvecs_in(:,3)==4)  | (dvecs_in(:,2)==2 & dvecs_in(:,3)==29) | (dvecs_in(:,2)==4  & dvecs_in(:,3)==17) | ...
                           (dvecs_in(:,2)==6 & dvecs_in(:,3)==28) | (dvecs_in(:,2)==9 & dvecs_in(:,3)== 8) | (dvecs_in(:,2)==11 & dvecs_in(:,3)==19));
        else
            % remove any leap days 
            deleted = find((dvecs_in(:,2)==2 & dvecs_in(:,3)==29));
        end            
            
        if (do_random)
            if (islogical(do_random))
                oldseed=rng(mod(round(86400*now()), 2^32));       % truly random...based on current clock.
            else
                oldseed=rng(do_random);         % use do_random as the seed.  This is so users can duplicate results later.
            end
            offsets=round(randn(size(deleted))*7);      % Normally distributed offsets w/ sigma of 7.  prob(abs(offsets)) > 34 is 3e-7.
            deleted = max(2,min(length(out_dnums)-1, deleted + offsets)); % And make sure we don't delete the first or last points.
            rng(oldseed);  % put randon number generator back where it was.
        end
        inserted = [];
        to_ix    = [];
        keepers  = true(npts_in,1);
        keepers(deleted) = false;
        
    else    

        if (incal_len == 360)
            inserted = find((dvecs_out(:,2)==2 & dvecs_out(:,3)==4)  | (dvecs_out(:,2)==2 & dvecs_out(:,3)==29) | (dvecs_out(:,2)==4  & dvecs_out(:,3)==17) | ...
                            (dvecs_out(:,2)==6 & dvecs_out(:,3)==28) | (dvecs_out(:,2)==9 & dvecs_out(:,3)== 8) | (dvecs_out(:,2)==11 & dvecs_out(:,3)==19));
        else
            % inserted are any leap days 
            inserted = find((dvecs_out(:,2)==2 & dvecs_out(:,3)==29));
        end

        if (do_random)
            if (islogical(do_random))
                oldseed=rng(mod(round(86400*now()), 2^32));       % truly random...based on current clock.
            else
                oldseed=rng(do_random);         % use do_random as the seed.  This is so users can duplicate results later.
            end
            offsets=round(randn(size(inserted))*7);      % Normally distributed offsets w/ sigma of 7.  prob(abs(offsets)) > 34 is 3e-7.
            inserted = max(2,min(length(out_dnums)-1, inserted + offsets)); % And make sure we don't delete the first or last points.
            rng(oldseed);  % put randon number generator back where it was.
        end
        to_ix   = (1:npts_out);
        to_ix(inserted) = [];
        keepers = logical([]);
        deleted = [];
    end
        
    
        % finally, remove any dates that were missing in the input
    
    if (~isempty(missing_dnums_in))
        if (incal_len == 360)  % can't use datevecs, so have a few hoops to jump through to get there.
            scaling = length(out_dnums_all)/length(in_dnums_all);
            [~,ix] = intersect(in_dnums_all,missing_dnums_in);
            ix = round(ix*scaling); 
        elseif (outcal_len == 360)
            scaling = length(out_dnums_all)/length(in_dnums_all);
            [~,ix] = intersect(in_dnums_all,missing_dnums_in);
            ix = round(ix*scaling);             
        else
                    % this is easy.  just delete the equivalent datevecs.
            missing_dvecs = datevec_cal(missing_dnums_in, incal);
                    % don't delete leap days if going from std to 365-day...
            if (outcal_len == 365)
                ix = find(missing_dvecs(2)==2 & missing_dvecs(3) == 29);
                if (~isempty(ix))
                    missing_dvecs(ix,:) = [];
                end
            end
            missing_dnums_out = datenum_cal(missing_dvecs, outcal);
            [~,ix] = intersect(out_dnums, missing_dnums_out);
        end
        if (~isempty(to_ix))
            to_ix(ix) = [];
            out_dnums(ix) = [];
        end
        if (~isempty(keepers))
            out_dnums(ix) = [];
            keepers(ix) = [];
        end
    end
    
        % finally finally, add fractional days back in.
    if (has_dayfracs)
        if (incal_len > outcal_len)
            out_dnums = out_dnums + dayfracs(keepers);
        else
            out_dnums(to_ix)  = out_dnums(to_ix)  + dayfracs;
            out_dnums(inserted) = out_dnums(inserted) + mean_dayfrac;
        end
    end
end

