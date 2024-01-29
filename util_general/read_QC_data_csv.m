function [vals, dates, min_date, max_date, nvals] = read_QC_data_csv(stnID, data_dir, start_date, end_date)
% reads QC data for station stnID from QC csv file
% returns dates, values, min & max dates, & # of valid datapoints.
%   all missing values are filled with NAs, so you are guaranteed to get entire
%   date range requested with no gaps.
%
%   Inputs:
%       stnIDs      station ID  of station to read
%       data_dir    folder where data is found.  leave blank for current folder
%       start_date  1st date to extract.  Use just year, 1950, or full date [1950,1,1].  If empty, returns all available
%                       data (with gaps filled with NAs).  Do not leave empty if reading multiple files.
%       end_date    last date or year to extract.
%       varName     variable name.  (not needed fo TTU CSC QC'd csv files, since variable is specified in data_dir.
%       ncfile      if data from QC station netcdf file, name of file.  If blank, will look for csv files with same
%                       name as stnID.
%
%   Outputs as labeled, with following notes:
%       dates       are matlab datenums (days since Jan 0, year 0).
%
%   if nargout is 1 (only 1 output argument), returns a table with datenums and values as variable "dates".
%       column labels are "dnum" and "val".
%   otherwise returns data in separate variables.
%       mindate, maxdate are returned as strings, yyyy-mm-dd.

    
    if (isstring(stnID)),    stnID=char(stnID);       end
    if (isstring(data_dir)), data_dir=char(data_dir); end
    
    fn=fullfile(data_dir, stnID);
    
    if (~exist('start_date','var') || isempty_s(start_date)), start_date=[]; end
    if (~exist(  'end_date','var') || isempty_s(  end_date)),   end_date=[]; end
    
    if (length(start_date) == 1), start_date = [start_date, 1, 1]; end
    if (length(  end_date) == 1),   end_date = [  end_date,12,31]; end
        % get rid of NAs in file.  This is faster than having readtable do it.
        % then read data in to table t.
    tmpname = tempfilename('/tmp','QC_csv','txt');
    cmd = sprintf('grep -v NA %s > %s', fn, tmpname);
    system(cmd);
    t = readtable(tmpname);
    system(sprintf('rm %s', tmpname));

    file_dnums = datenum(t.(1),t.(2),t.(3));

        % get start, end date if not provided
    if (isempty_s(start_date))
        start_dnum = nanmin(file_dnums);
        
        if (isempty(start_dnum))        % bail out if no valid data.
            dates=[];
            vals=[];
            min_date=[];
            max_date=[];
            nvals=0;
            return;
        end

    else
        start_dnum = datenum(start_date);        
    end
    
    if (isempty_s(end_date))
        end_dnum=nanmax(file_dnums);
    else
        end_dnum = datenum(end_date);
    end
    
    day0 = start_dnum - 1;
    dates=(start_dnum:end_dnum)';
    npts = end_dnum - day0;
    vals = nan(length(dates),1);

    v = t.(4);
    if (~isnumeric(v))
        fprintf(2, 'error reading file %s\n', stnID); 
        return;
    end
    ix = file_dnums - day0;
    keepers = ix>0 & ix <= npts;
    ix = ix(keepers);
    v = v(keepers);
    vals(ix) = v;
    
        % set dates for all days w/ no data to nan's
    kept = ~isnan(vals);
    
    if (nargout == 1)
        return;
    elseif (nargout == 2)
        dates=table(dates,vals,'VariableNames',{'datenum','val'});
        return;
    else
        
        nvals=sum(~isnan(vals));            

        if (nvals == 0)
%             fprintf(2, '%s: no valid data in %d lines of input\n', fn, size(t,1));
            min_date = [];
            max_date = [];
        else
            imin=find(kept,1);
            imax=find(kept,1,'last');

            min_date = string(datestr(dates(imin),'yyyy-mm-dd'));
            max_date = string(datestr(dates(imax),'yyyy-mm-dd'));
        end
    end    
end
