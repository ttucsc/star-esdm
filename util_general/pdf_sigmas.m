function [sigmas,mus] = pdf_sigmas(pdfs, bins, mus, do_normalize)
%
%   This returns std. deviations of pdfs.
%   Also calculates means, if empty or not provided
%
%   pdfs are along dimension 1 (down each column).

    if (~exist('mus','var')), mus=[]; end
    if (~exist('do_normalize','var') || isempty(do_normalize)), do_normalize=false; end

    if (do_normalize)
        pdfs = pdfs ./nansum(pdfs);
    end
    [sigmas,mus] = pdf_moment(pdfs, 2, bins, mus, false);
    sigmas = sqrt(sigmas);
    
end

