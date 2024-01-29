function dnum365 = datenum365(dvec)
%   returns days since [0,0,0,0,0,0] in 365-day years.
%       NOTE:  bug in matlab's days365 for start_date of 0/0/0

    if (isstring(dvec))
        dvec = datevec(dvec);
    end
    ndays = size(dvec,1);
    
    dnum = datenum(dvec);
    frac = mod(dnum,1.0);       % to get hh,mm,ss
    dnum = floor(dnum);

            % bug in days365 when using date of 0 !!!! days365 crashes
%     dvec0 = zeros(ndays,1);  % days since 0, 0, 0000
%     dnum365 = days365(dvec0,dnum);
    
    dnum1 = ones(ndays,1);  % using days since Jan 1, 0000 and then adding 1
    dnum365 = days365(dnum1,dnum) + frac + 1;
end