function [date_range, date_ranges] = ncdf_get_date_range(fnames)
%returns date range of valid data in one or more filenames
% outputs:
%   date_range      overall date range, min & max, [yyyy,mm,dd; yyyy, mm, dd]
%   date_ranges     date range for each file, size:  2 x 3 x nfiles  
    
    nfiles = numel(fnames);
    date_ranges = zeros(2,3,nfiles);
    dnums = nan(nfiles,2);
    for i=1:nfiles
        [~,calendar,startvec, days_since] = ncdf_get_time_info(fnames{i});
        dnum1 = datenum_cal(startvec, calendar)+days_since(1);
        dnum2 = datenum_cal(startvec, calendar)+days_since(end);
        dvecs = datevec_cal([dnum1;dnum2],calendar);
        date_ranges(:,:,i) = dvecs(:,1:3);
        dnums(i,:) = datenum(dvecs);
    end
    if (nfiles == 1)
        date_range = squeeze(date_ranges);
    else
        date_range = [datevec(min(dnums(:,1))); datevec(max(dnums(:,2)))];
    end
end

