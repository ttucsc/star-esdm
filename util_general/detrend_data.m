function [det_y, avg, trend_params, trend] = detrend_data(y, order)
% [det_y, avg, trend] = detrend_data(y, order)
% calculates trend from original data, and returns detrended signal, avg value, and trend.
%   Inputs:
%       y           data to detrend
%       order       polynomial order to use (3 is good for 150-year data.  2nd order, parabolic, doesn't fit data well enough.
%                                           (1 is good for 40- or 50-year historical period)
%                                            0 for straight average
%                                           -1 :  data already detrended.  
%
%   Outputs:
%       det_y         original signal, with trend subtracted
%       avg           mean value of y
%       trend_params  polyfit params to recreate the long term trend.  trend = polyval(p,(0:(len-1))/365,[],mu);
%       trend         daily trend over entire period 
%                       Note:  trend & trend_params are zero-mean...i.e. with mean removed..

        % do polynomial fit to data of specified order,
        % then subtract the trend to separate out the long term trend and detrended signal
    if (order >= 0)
        x=(0:(length(y)-1))'/365;       % x is in (fractions of) years
        avg=mean(y);
        y = y - avg;
%         [p,~,mu] = polyfit(x,y,order);
%         trend = polyval(p,x,[],mu);
        p = polyfit(x,y,order);
        trend = polyval(p);
%         p2=polyfit(x,y,order);
%         trend2 = polyval(p2,x);
%         figure;
%         plot(x,trend,x, trend2+.01, x, trend-trend2); legend('trend','trend2','diff'); title(sprintf('%.20f',max(abs(trend-trend2))));
        det_y = y - trend;
    else                        
        det_y = y;
        p = 0;
        avg = 0;
        trend = zeros(size(y));
        mu=0;        
    end
    
%     trend_params = struct('p',p,'len',length(y),'mu',mu);     % to reconstruct trend:  trend = polyval(p,0:(len-1)/365,[],mu);
    trend_params = struct('p',p,'len',length(y),'x','x=(0:(len-1))/365; % x is in (fractions of) years);  to reconstruct trend:  trend = polyval(p,0:(len-1)/365);');     % to reconstruct trend:  trend = polyval(p,0:(len-1)/365);
end    

