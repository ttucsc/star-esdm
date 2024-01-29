function [pdfs, CDFs, clim_sigs, slvrsigs] = kde_surf(obj, hist_tbl, bins)
%
%   This does (only) kernel-density smoothing on the histogram table,
%   using the parameters in obj.RP (the run parameters struct).
%
%   hist_tbl can be a single (nbins x ndays) histogram surface, or a set of rolling histograms (nbins x ndays x nsets)
%
%   It calculates the climatology of the daily histograms, then does Silverman's rule-of-thumb gaussian smoothing
%   on the daily histograms, but using the climatology-smoothed standard deviations.
%
%   Output is a set of daily pdf's (& CDFs) based on the original histogram table.
%
%   Resulting pdf's can either be gaussian-smoothed along the day-of-year axis, or can be smoothed by calculating a
%   dense set of probability iso-lines and taking their climatology.  This latter approach will avoid fitting to the
%   noise of limited sampling size.  Gaussian smoothing produces bumps in the tails of the pdfs when the random
%   sampling happens to result in several points in the tails near the same date.  It also smooths along fixed 
%   temperature or precipitation values, which distorts the probabilties during transition times.  
%   The alternative is to smooth along fixed probability values, which requires mapping to a probability space,
%   smoothing in that space, and mapping back to temperature/precip space.
%
%   Calculates CDFs from pdfs only if nargout > 1.
   

    if (~exist('bins','var') || isempty(bins))
        bins = obj.RP.bins;
    end
%   dx = obj.RP.bins(2) - obj.RP.bins(1);
    dx = bins(2)-bins(1);

    [nbins,yrlen,nsets] = size(hist_tbl);
    ss_days = max(1,60*yrlen/365);      % for silverman's sigma calculations.  We'll smooth via low-pass filtering, with 
                                        % approx. equivalent of 2 months (60 days).  so we'll multiply the number of 
                                        % datapoints by ss_days.  Otherwise we'll be oversmoothing.

        % normalize histograms to sum to 1

    n = sum(hist_tbl);
    pdf_tbl = hist_tbl ./ n;
    n = squeeze(n);         % get rid of singleton 1st dimension.  (we needed it 1 x yrlen x nsets for previous statement)
    if (isrow(n)), n=n'; end
    if (isrow(bins)), bins=bins'; end

    GS         = zeros(nbins,yrlen);    % for gaussian kernels in fourier domain.
    pdfs = zeros(size(hist_tbl));
    if (nargout > 1)
        CDFs = zeros(size(hist_tbl));
    end

    for j=1:nsets
        pdfs_t = squeeze(pdf_tbl(:,:,j));
        nn   = n(:,j)*ss_days;
             
%         mus = sum(bins .* pdfs_t);                       % mean daily value, s/b around midpoint of bins.
% 
%         sigs = sqrt(sum(bins.^2 .* pdfs_t) - mus(1,:).^2); % sigmas will be noisy, so we need to low-pass-filter them
%         
        sigs = pdf_sigmas(pdfs_t, bins, [], false);
        clim_sigs = obj.DA_climatology(sigs, obj.RP.anom_nterms(1), obj.RP.anom_sig_terms(1), yrlen)';  % we use these sigmas in calculating the silverman ROT sigmas.

        silversigs = 1.06 * clim_sigs .* (nn'.^(-.2))/dx;  % silverman's rule of thumb.  divide by dx so units is bins, not the data's units.
        silvermin = silversig_min(nbins);                   % returns a minimum value of silversigs for the given number of bins.
        if (min(silversigs) < silvermin)
            obj.DP.warn_log(sprintf('warning: kde_surf:  silverman sigmas too small;  adjusted to %.4f (min was %.4f)  location: %s %s\n', silvermin, min(silversigs), obj.DP.stnName, obj.DP.stnID));
        end
        silversigs = max(silversigs, silvermin);              % If silversigs are too small, the filter ends up distorted because the gaussian in the fft doesn't get small enough.
        slvrsigs = silversigs*dx;                             % elsewhere we need silversigs in the data's units, so multiply by dx.
            % make fourier-domain gaussians based on silversigs
        for i=1:yrlen, GS(:,i) = GAUSS_fft(nbins, silversigs(i), true)'; end
            % convolve with pdfs
        pdfs_t = abs(ifft(GS .* fft(pdfs_t)));     % s/b purely real, but may have tiny imag. part, so take abs of inverse fft here.

            % renormalize so we have pdfs still.
        ss=sum(pdfs_t);
        pdfs_t = pdfs_t./ss;

            % put the results into the output pdf & CDF matrices.
        if (nsets > 1)
            pdfs(:,:,j) = pdfs_t;
        else
            pdfs = pdfs_t;
        end
                % Now integrate to create CDFs
        if (nargout > 1)
            if (nsets > 1)
                CDFs(:,:,j) = cumsum(pdfs(:,:,j));
            else
                CDFs = cumsum(pdfs);
            end
        end
    end

end
