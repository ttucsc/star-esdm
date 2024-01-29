function outrange = min_yr_range(varargin)
% returns inner range (smallest range) of years defined by input ranges, where input ranges are 2-element pairs,
% Ranges can be [start_year, end_year], [start_year; end_year], or 2 datevecs [ yy,mm,dd[,hh,mm,ss]; yy,mm,dd[,hh,mm,ss] ]
% return value is [largest_start_year, smallest_end_year];

    outrange = make_range(varargin{1});
    for i=2:length(varargin)
        range2 = make_range(varargin{i});
        outrange = [nanmax(outrange(1),range2(1)),nanmin(outrange(2),range2(2))];
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


