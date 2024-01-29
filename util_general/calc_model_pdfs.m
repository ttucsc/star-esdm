function [ab_pdf, ab_cdf] = calc_model_pdfs(alpha, beta, mdl_sigs, mdl_mus, x)
% returns calculated pdf, using skew_generalized_normal, rescaled to match mdl_stats' mean and standard deviation.
%   Inputs:
%       alpha, beta     parameters for skew_generalized_normal
%       mdl_stats       struct with mu, sigma, skew, excess_kurtosis for model
%       x               x-locations (bins) to evaluate pdf and cdf at.

    nc = length(mdl_sigs);
    nx = length(x);
    ab_pdf = zeros(nx,nc);

    [ab_pdft,ab_stats, ~] = skew_generalized_normal(0,alpha,beta,x);  %  generate the pdf, and get its stats.
    
    for i=1:nc
    
        % scaling between skew_generalized_normal and model pdfs.
        k = ab_stats.sigma / mdl_sigs(i);

        ab_pdf(:,i) = interp1(x/k-ab_stats.mu/k, ab_pdft, x-mdl_mus(i),'makima',0);    % resample at desired bins and adjust mean to match model's mean. 
        ab_pdf(:,i) = ab_pdf(:,i) / nansum(ab_pdf(:,i));  % rescale so pdf sums to 1.
    end
    ab_cdf = cumsum(ab_pdf);           % create cdf
    ab_cdf = ab_cdf ./ ab_cdf(end,:);     % make sure cdf ends exactly at 1.

end
