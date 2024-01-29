function dvec = datevec360(dnum360)
% returns datevec [yyyy,mm,dd,hh,mm,ss] from a dnum360 -- a datenum in 360-day years (inverse of datenum360(dvec)

    % this should be able to return a datevec if the input is a string.
    % solution to that will take a little work, since datevec(...) returns the wrong datevec for feb 29 (non-leap year) or feb 30
    %

    dvec = datevec(mod(dnum360,1));    % regular datevec.  gets us the hh, mm, ss fields
    
%       now get yr, month, day
    dnum360_1 = floor(dnum360-1);
    yr = floor(dnum360_1/360);
    dnum360_1 = mod(dnum360_1,360);
    mo = floor(dnum360_1/30) + 1;
    da = mod(dnum360_1,30) + 1;
    
    dvec(:,1:3) = [yr,mo,da];
end