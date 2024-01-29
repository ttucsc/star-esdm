function yesno = isleapday( d )
    %isleapday( d )  returns true if d is a leapday, false if not
    %   returns boolean vector or matrix if d is vector or matrix.
    if (size(d,2)==3 || size(d,2)==6)
        yesno = d(:,2) == 2 && d(:,3)==29;
    else
        dv = datevec(d);
        yesno = (dv(:,2)==2 & dv(:,3)==29);
    end

end

