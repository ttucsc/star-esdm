function [y,sgn_stats, sgn_cdf] = skew_generalized_normal(mu, alpha, beta, x)
% Returns a skewed, generalized normal distribution
%   Based on Azzalini, A. (1985). "A class of distributions which includes the normal ones". Scandinavian Journal of Statistics. 12: 171?178.
%   as described simply in https://en.wikipedia.org/wiki/Skew_normal_distribution
%   which notes that the skew_normal can be generalized to use any symmetric pdf & cdf to produce a skewed distribution.
%   Here we use the generalized normal distribution  (see https://en.wikipedia.org/wiki/Generalized_normal_distribution)
%   which allows us to add a specific amount of kurtosis.
%
%   Uses:  generalized_normal(mu, alpha, beta, x)
%       
%   Inputs:
%       mu          location, real
%       alpha       scale (positive, real), controls standard deviation
%       beta        shape (positive, real), controls kurtosis
%       x           range over which to evaluate generalized_normal. 
%       gamma       skew  (real)  controls skew
 %
%       y           skewed gaussian
%       sg_stats    skewed-gaussian stats struct, with fields mu, sigma, skew and kurt.
%                               note:  xkurt is excess kurtosis (kurtosis-3)  xkurt==0 for pure gaussian (alpha=0)

%     aa = abs(alpha);
%     [y1, ~, y2] = generalized_normal(mu, alpha, beta, x); 
%     if (alpha < 0), y2 = 1-y2; end
    y1 = generalized_normal(mu, 1, beta, (x-mu)); 
    y2 = 1/2*(1+erf(alpha*(x-mu)/sqrt(2)));
    
    y = 2 * y1 .* y2;
    
    if (nargout > 1)
        sgn_stats = struct;
        [sgn_stats.mu, sgn_stats.sigma, sgn_stats.skew, sgn_stats.xkurt] = pdf_stats(y,x,[],true);
        
        if (nargout > 2)
            sgn_cdf = cumsum(y);
            if (sgn_cdf(end) ~= 1.0)
                sgn_cdf = sgn_cdf / sgn_cdf(end);
            end
        end
    end
end

