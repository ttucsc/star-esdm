function [vals, min_date, max_date, nvalid, mpts, file_startdate, file_enddate] = read_csv(fnm, varColumn, timeColumn, tstamps, timestr, treatasempty, isUTC, toUTC, TZ)

% tstamps is matlab datenums.
% timeColumn is either name of single column with time column name, assumed here as "yyyymmddhh"
% Or can be names of year, month, day, hour, (and optionally minute and second columns)
    warning('off','MATLAB:table:ModifiedAndSavedVarnames');
    tbl = readtable(fnm,'TreatAsEmpty',treatasempty);
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
    
    nvalid = sum(~isnan(vals));
    min_ix = find(~isnan(vals),1);
    max_ix = find(~isnan(vals),1,'last');
    min_date = tstamps(min_ix);
    max_date = tstamps(max_ix);
    mpts = max_ix - min_ix + 1;

end
