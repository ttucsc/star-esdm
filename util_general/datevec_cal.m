function dvec = datevec_cal(dnum, calendar)
% dvec = datevec_cal(dnum, calendar) : converts datenums or datestr's into datevecs.  Datestr's should be of the form yyyy-mm-dd
% 
    if (ischars(dnum))
        dvec = datevec(dnum);
        if (calendar_length(calendar)==360)
            dvec(dvec(:,3)>30,3)=30;
        end
    elseif (calendar_length(calendar) == 360)
        dvec = datevec360(dnum);
    elseif (calendar_length(calendar) == 365)
        dvec = datevec365(dnum);
    else
        dvec = datevec(dnum);
    end
end