function [slope,intcp] = quick_lsq(y,x, sumx, sumxsq)
%   calculates lsq slope, intercept for (x,y)
%   if x is 0, then x <- 0:len(y)-1
%   if x is 1 or missing, then x <- 1:len(y)
%   otherwise, x is vector of length(y);
%

    n = length(y);
    y=squeeze(y);
    if (~exist('x','var') || isempty(x) || (length(x)==1 && x==1))
        x=1:n;
        if (~exist('sumx','var') || isempty(sumx))
            sumx = round(.5*(n+.5)^2);
        end
    elseif (length(x)==1 && x==0)
        x=0:(n-1);
        if (~exist('sumx','var') || isempty(sumx))
            sumx = round(.5*(n-.5)^2);
        end
    elseif (~exist('sumx','var') || isempty(sumx))
        sumx = sum(x);
    end
    if (~exist('sumxsq','var') || isempty(sumxsq))
        sumxsq = sum(x.*x);
    end
    
    if (iscolumn(y) ~= iscolumn(x)), x=x'; end
    
    sumy = sum(y);
%   sumysq = sum(y.*y);
    
    sumxy = sum(x.*y);
    
    slope = (n*sumxy - sumx*sumy)/(n*sumxsq - sumx*sumx);
    
    intcp = sumy/n - slope*sumx/n;

end

