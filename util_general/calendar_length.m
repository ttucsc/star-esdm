function cal_length = calendar_length( cal )
    % yr_length = calendar_length( cal )
    %  
    %   Returns length of year, based on string calendar
    %   Either:  360, 365, or 365.25
    cal365={'day-365','day365','day_365','365-day','365_day','365day','noleap','365days','365_days','days365','days_365'};
    calstd={'std','standard','julian','gregorian','proleptic','proleptic_gregorian'};
    cal360={'day-360','day_360','day360','360-day','360_day','360day'};

    if (isnumeric(cal))
        cal_length = cal;
    elseif (any(strcmpi(cal365,cal)))
        cal_length = 365;
    elseif (any(strcmpi(calstd,cal)))
        cal_length = 365.25;
    elseif (any(strcmpi(cal360,cal)))
        cal_length = 360;
    else
        throw(MException('ICSF:BAD_CALENDAR',sprintf('error:  unknown calendar type: %s', cal)));
    end
    
    
end

