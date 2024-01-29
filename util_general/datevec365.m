function dvec = datevec365(dnum365)
% returns datevec [yy,mm,dd,hh,mm,ss] from a dnum365 -- a datenum in 365-day years (inverse of datenum365(dvec)
%     month_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

    if (isstring(dnum365))
        dvec = datevec(dnum365);
        return;
    end
     mdoy365 = [31,59,90,120,151,181,212,243,273,304,334,365];
    
    dvec = datevec(mod(dnum365,1));    % regular datevec.  gets us the hh, mm, ss fields
%       now get yr, month, day
    dnum365_1 = floor(dnum365-1);
    yr = floor(dnum365_1/365);
    dnum365_1 = mod(dnum365_1,365);
    
    doy365_1 = mod(dnum365_1,365);

    for i=1:length(dnum365_1(:))
        mo = find(mdoy365>doy365_1(i),1);
        if (mo == 1)
            da = doy365_1(i) + 1;
        else
            da = doy365_1(i) - mdoy365(mo-1)+1;
        end     
        dvec(i,1:3) = [yr(i),mo,da];
    end
end