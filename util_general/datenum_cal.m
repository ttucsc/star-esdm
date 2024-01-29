function dnum_cal = datenum_cal(dvec, calndar)

    if (ischar(dvec) || isstring(dvec) || (iscell(dvec) && (ischar(dvec{1}) || isstring(dvec{1}))))
        dvec = datevec_cal(dvec,calndar);
    end
    if (calendar_length(calndar) == 360)
        dnum_cal = datenum360(dvec);
    elseif (calendar_length(calndar) == 365)
        dnum_cal = datenum365(dvec);
    else
        dnum_cal = datenum(dvec);
    end
end