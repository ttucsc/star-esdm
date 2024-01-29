function [kurtosis, mus, sigmas] = pdf_kurtosis(pdfs, bins, mus, sigmas, do_normalize)
%
%   This returns Excess Kurtosis of pdfs.  
%   It also calculates means and sigmas, if empty or not provided.

    if (~exist('mus','var')), mus=[]; end
    if (~exist('sigmas','var')), sigmas=[]; end
    if (~exist('do_normalize','var') || isempty(do_normalize)), do_normalize=false; end
    
    if (do_normalize)
        pdfs = pdfs ./sum(pdfs);
    end
    
    if (isempty(mus)),    mus    = pdf_moment(pdfs, 1, bins, mus, false);   end
    if (isempty(sigmas)), sigmas = pdf_sigmas(pdfs, bins, mus, false);      end
        
    kurtosis = pdf_moment(pdfs, 4, bins, mus, false);
    kurtosis = kurtosis ./ (sigmas .^ 4) - 3;
        
end

