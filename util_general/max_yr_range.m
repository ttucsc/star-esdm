function outrange = max_yr_range(varargin)
% returns outer range (largest range) of years defined by input ranges, where input ranges are 2-element pairs,
% Ranges can be [start_year, end_year], [start_year; end_year], or 2 datevecs [ yy,mm,dd[,hh,mm,ss]; yy,mm,dd[,hh,mm,ss] ]

    outrange = [];
    for i=1:length(varargin)
        if (~isempty(varargin{i}))
            if (isempty(outrange)), outrange = nan(1,2); end
            range2 = make_range(varargin{i});
            outrange = [nanmin(outrange(1),range2(1)),nanmax(outrange(2),range2(2))];
        end
    end
end

function myrange = make_range(v)
    if (isempty(v))
        myrange = nan(1,2);
    elseif (numel(v) <= 2)
        myrange = [v(1);v(end)];
    elseif (isrow(v))
        myrange = [v(1);v(end)];
    else
        myrange = [v(1,1);v(end,1)];
    end
end


