function [ab_pdf, ab_stats, ab_cdf] = calc_model_pdfs(alpha, beta, mdl_stats, x)
% returns calculated pdf, using skew_generalized_normal, rescaled to match mdl_stats' mean and standard deviation.
%   Inputs:
%       alpha, beta     parameters for skew_generalized_normal
%       mdl_stats       struct with mu, sigma, skew, excess_kurtosis for model
%       x               x-locations (bins) to evaluate pdf and cdf at.
            
    [ab_pdft,ab_stats, ~] = skew_generalized_normal(0,alpha,beta,x);  %  generate the pdf, and get its stats.
    
        % scaling between skew_generalized_normal and model pdfs.
    k = ab_stats.sigma / mdl_stats.sigma;
            
%     ab_bins = x * k;     % rescale the mdl_bins
%     [bb_pdft,ab_stats, ~] = skew_generalized_normal(0,alpha,beta,ab_bins);  %  generate the skew_generalized_normal pdf using scaled bins.
    ab_pdf = interp1(x/k-ab_stats.mu/k, ab_pdft, x-mdl_stats.mu,'makima',0);    % resample at desired bins and adjust mean to match model's mean. 
    ab_pdf = ab_pdf / nansum(ab_pdf);  % rescale so pdf sums to 1.
    ab_cdf = cumsum(ab_pdf);           % create cdf
    ab_cdf = ab_cdf / ab_cdf(end);     % make sure cdf ends exactly at 1.

    [ab_stats.mu, ab_stats.sigma, ab_stats.skew, ab_stats.xkurt] = pdf_stats(ab_pdf,x, [], true);

end
