function cal_type = calendar_type( cal_length )
    % yr_length = calendar_length( calendar )
    %  
    %   Returns "360-day", "365-day","standard", "366-day", or "unknown" based on cal_length
    %   input:  Either:  array of calendar lengths.  usually 360, 365, or 365.25 or 366.  
    %   output:  string (see below).  If not valid length, returns "unknown"
    
    cal_lengths = [360, 365, 365.25, 366];
    cal_types=["360-day","365-day","standard","366-day"];

    ncals = numel(cal_length);
    cal_type = strings(size(cal_length));
    for i=1:ncals
        ix = find(cal_length(i)==cal_lengths,1);
        if (isempty(ix))
            if (cal_length(i) >=365.2 && cal_length(i) <= 365.3)
                cal_type(i) = cal_types(3);     % standard.
            else
                cal_type(i) = "unknown";
            end
        else                
            cal_type(i) = cal_types(ix);
        end
    end
    
    
end
