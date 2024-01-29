function [y,sg_stats, gn_cdf] = generalized_normal(mu, alpha, beta, x)
% returns generalized_normal(alpha) type-1, as described at https://en.wikipedia.org/wiki/Generalized_normal_distribution
%
%   It returns a symmetric distribution with kurtosis.
%
%       If beta = 2, it produces a gaussian with standard deviation of sqrt(alpha/2).
%       range of x should be at least +/- 5*alpha about the mean, mu
%
%   Inputs:
%       mu          location, real
%       alpha       scale (positive, real)  controls standard deviation
%       beta        shape (positive, real)  controls kurtosis
%       x           range over which to evaluate generalized_normal.       
%
%       y           generalized gaussian
%       sg_stats    generalized-gaussian stats struct, with fields mu, sigma, skew and kurt.
%                               note:  xkurt is excess kurtosis (kurtosis-3)  xkurt==0 for pure gaussian (alpha=0)
%       gn_cdf      cdf of y (cumsum(y)), normalized so y(end) == 1.
%
    dx = x(2)-x(1);
    y = beta/(2*alpha*gamma(1/beta)) * exp(-((abs(x-mu))/alpha).^beta);
    y = y / (sum(y)*dx);
    
    if (nargout > 1)
        sg_stats = struct;
        [sg_stats.mu, sg_stats.sigma, sg_stats.skew, sg_stats.xkurt] = pdf_stats(y,x,[],true);
        
        if (nargout > 2)
            gn_cdf = cumsum(y);
            if (gn_cdf(end) ~= 1.0)
                gn_cdf = gn_cdf / gn_cdf(end);
            end
        end
    end

end

