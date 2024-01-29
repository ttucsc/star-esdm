function find_model_pdf_coeffs_verify(figbase, ab_dx)

    if (~exist('figbase','var') || isempty(figbase)), figbase=101;  end
    if (~exist('ab_dx','var')   || isempty(ab_dx)),   ab_dx= .1;    end
    
%     [models,vars] = get_model_info();    
%     nmodels = length(models);
%     nvars   = length(vars);
%     sk = load("skew_kurt.mat");     % generated skewed & kurtosed pdfs
    sk = load("model_pdf_coefs.mat");     % table with best sk_normal coefs for each model/var.
    
%     alpha = nan(nmodels, nvars) ;
%     beta  = nan(nmodels, nvars);
%     sigma = nan(nmodels, nvars);
    alpha = sk.alpha ;
    beta  = sk.beta;
%     sigma = sk.sigma;
    models = sk.models;
    vars = sk.vars;
    
%     model_stats = cell(nmodels, nvars);
    
    for m=1:length(models)
        model = models(m);

        for v=1:length(vars)
            varname = vars(v);
            cdfname = sprintf('basic_cdfs/cdf_basic.%s.%s.mat',model, varname);
            load(cdfname,'mdl_stats','mdl_cdf','mdl_bins','okflags')
            
            ab_bins = -50:ab_dx:50;

            fprintf('%-20s %-8s \tmean %8.4f \tsigma %8.4f \tskew %8.4f \tkurt %8.4f\n', model, varname, mdl_stats.mu, mdl_stats.sigma, mdl_stats.skewness, mdl_stats.excess_kursosis);
            mdl_pdf = diff([0,mdl_cdf]);    % get the model's cdf
            dx = mdl_bins(2)-mdl_bins(1);
            nbins = length(mdl_bins);

            %flags = ~isnan(mdl_pdf);
            flags = okflags;
            ixmin = max(1, find(flags,1)-10);
            ixmax = min(nbins, find(flags,1,'last')+10);
            xlims = mdl_bins([ixmin,ixmax]);

    %       npts = length(mdl_bins);
            nrml_pdf = normpdf(mdl_bins,mdl_stats.mu, mdl_stats.sigma);
            nrml_cdf = normcdf(mdl_bins,mdl_stats.mu, mdl_stats.sigma);

%             % skew & kurtosis are independent of scaling.
%                 % find closest match in skew & kurtosis.  
%             skewdiff = mdl_stats.skewness - sk.skews;
%             kurtdiff = mdl_stats.excess_kursosis - sk.kurts;
%             dist = sqrt(skewdiff.*skewdiff + kurtdiff.*kurtdiff);
% 
%             abix = find(dist()==min(dist(:)));
%             [ixalpha,ixbeta] = ind2sub([2001,2001],abix);
% 
%             fprintf('closest match: alpha(%d): %8.4f  beta(%d): %8.4f\n', ixalpha, sk.alphas(ixalpha), ixbeta, sk.betas(ixbeta));
%             fprintf('closest match: std: %8.4f skew %8.4f kurt:  %8.4f\n',  sk.sigmas(ixalpha,ixbeta), sk.skews(ixalpha,ixbeta), sk.kurts(ixalpha, ixbeta))
% 
%             alpha(m,v) = sk.alphas(ixalpha);
%             beta(m,v)  = sk.betas(ixbeta);
%             sigma(m,v) = sk.sigmas(ixalpha,ixbeta);
%             
%             model_stats{m,v} = mdl_stats;
            
%             [ab_pdf, ab_stats, ab_cdf] = calc_model_pdf(alpha(m,v), beta(m,v), model_stats{m,v}, mdl_bins);
            [ab_pdf, ab_stats, ab_cdf] = calc_model_pdf(alpha(m,v), beta(m,v), mdl_stats, ab_bins);
            
%             ab_bins = mdl_bins * k;
% 
%             [ab_pdft,ab_stats, ~] = skew_generalized_normal(0,sk.alphas(ixalpha),sk.betas(ixbeta),ab_bins);
%             %ab_pdf = interp1(mdl_bins-ab_stats.mu/k, ab_pdft, mdl_bins+mdl_stats.mu,'makima',0);
%             ab_pdf = interp1(mdl_bins-ab_stats.mu/k, ab_pdft, mdl_bins-mdl_stats.mu,'makima',0);
%             ab_pdf = ab_pdf / nansum(ab_pdf);
%             ab_cdf = cumsum(ab_pdf);
%             ab_cdf = ab_cdf / ab_cdf(end);

            fprintf('%-20s %-8s \tmean %8.4f \tsigma %8.4f \tskew %8.4f \tkurt %8.4f\n', 'calculated pdf',' ',ab_stats.mu, ab_stats.sigma, ab_stats.skew, ab_stats.xkurt)

            figure(figbase+m-1);
            subplot(2,3,1+3*(v-1));
            plot(mdl_bins, nrml_pdf,'g','linewidth',2);
            hold on;
            plot(mdl_bins, mdl_pdf/dx, 'b',ab_bins, ab_pdf/ab_dx,'r');
            hold off;
            xlim(xlims);
            title(sprintf('%s %s',model,varname),'interpreter','none');
            legend('normal','model','generated');
            grid on;

            subplot(2,3,2+3*(v-1));
            semilogy(mdl_bins, nrml_pdf,'g','linewidth',2);
            hold on;
            semilogy(mdl_bins, mdl_pdf/dx, 'b', ab_bins, ab_pdf/ab_dx,'r');
            hold off;
            xlim(xlims);
            title(sprintf('%s pdf',varname));
            grid on;

            subplot(2,3,3+3*(v-1));
            semilogy(mdl_bins, nrml_cdf,'g',mdl_bins, 1-nrml_cdf,'g','linewidth',2);
            hold on;
            semilogy(mdl_bins, mdl_cdf, 'b', mdl_bins,1-mdl_cdf,'b', ab_bins, ab_cdf,'r', ab_bins, 1-ab_cdf,'r');
            hold off;
            xlim(xlims);
            title(sprintf('%s cdf, 1-cdf',varname));
            grid on;
        end
    end
    
%     README = get_README(); %#ok<NASGU>
%     save('model_pdf_coefs','models','vars','alpha','beta','sigma','mdl_stats','README');
%     
end

% function [models,vars] = get_model_info()
% 
%     finfo=dir("basic_cdfs/*.mat");
%     fnames = string({finfo.name});
%     nfnames = length(fnames);
%     models = strings(nfnames,1);
%     vars   = strings(nfnames,1);
%     for i=1:length(fnames)
%         fn = split(fnames(i),'.');
%         models(i) = string(fn{2});
%         vars(i)   = string(fn{3});
%     end
%     
%     models = unique(models);
%     vars   = unique(vars);
%         % exclude Tmax and Tmin, which are for station data only, but are duplicated for stations with tasmax and tasmin
%     ix = find(strcmpi(vars,"Tmax"));
%     if (~isempty(ix)), vars(ix)=[]; end
%     ix = find(strcmpi(vars,"Tmin"));
%     if (~isempty(ix)), vars(ix)=[]; end
% end

% function README = get_README()
%     README = ['This file contains the coefficients for estimating the extended pdfs and cdfs for each model/var combination.  ' ...
%               'It was generated with the function find_model_pdf_coeffs( ).'  ...
%               'To generate an estimate of the overall pdf and cdf for a model & variable, call calc_model_pdf( ) '...
%               'with the alpha, beta, sigma and pdf_stats struct for the model & var.  '...
%               'The output are a pdf and cdf with the same standard deviation, skew and kurtosis of the model & var''s overall pdf ' ...
%               'from combining the normalized model output for the historical period 1950-2005 for multiple gridcells across '...
%               'N. America for each model & variable.  See the comments in find_model_pdf_coeffs( ) and calc_model_pdf( ) ' ...
%               'for more details'];
% end

