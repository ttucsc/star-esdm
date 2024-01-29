function dnum360 = datenum360(dvec)
% %   returns days since [0,0,0,0,0,0] in 360-day years.
% %   NOTE:  results are wrong for dvec of [0,0,0] or earlier, due to bug in days360.
%     ndays = size(dvec,1);
%     
%     dnum = datenum(dvec);
%     frac = mod(dnum,1.0);
%     dnum = floor(dnum);
% 
%             % bug in days360 when using date of 0 !!!! days is off by 30!
% %     dvec0 = zeros(ndays,1);  % days since 0, 0, 0000
% %     dnum360 = days360(datenum(dvec0),datenum(dvec));
%     
%     dnum1 = ones(ndays,1);  % using days since Jan 1, 0000 and then adding 1
%     dnum360 = days360(dnum1,dnum) + frac + 1;
%                                 % make datenum360 act the same way datenum does for day 0 of any month.
%     z = dvec(:,3)==0;           % days360 doesn't work as expected for day 0's of any months
%     dnum360 = dnum360 - z;      % subtracting a boolean vector which is 1 wherever the day of the month is zero.
% end

    % this should be able to return a datevec if the input is a string.
    % solution to that will take a little work, since datevec(...) returns the wrong datevec for feb 29 (non-leap year) or feb 30
    %
    if (size(dvec,2)==3)
        dnum360=360*dvec(:,1) + 30*(dvec(:,2)-1) + dvec(:,3);
    elseif (size(dvec,2)==4)
        dnum360=360*dvec(:,1) + 30*(dvec(:,2)-1) + dvec(:,3) + dvec(:,4)/24;
    elseif (size(dvec,2)==5)
        dnum360=360*dvec(:,1) + 30*(dvec(:,2)-1) + dvec(:,3) + dvec(:,4)/24 + dvec(:,5)/1440;
    elseif (size(dvec,2)==6)
        dnum360=360*dvec(:,1) + 30*(dvec(:,2)-1) + dvec(:,3) + dvec(:,4)/24 + dvec(:,5)/1440 + dvec(:,6)/86400;
    else
        error("datenum360: error:  invalid number of columns for input.  must be 3-6 columns representing [yyyy, mm, dd, HH, MM, SS].  Size is %d x %d", size(dvec));
    end
end