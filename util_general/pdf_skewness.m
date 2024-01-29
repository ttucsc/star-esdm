function [skewness,mus, sigmas] = pdf_skewness(pdfs, bins, mus, sigmas, do_normalize)

    if (~exist('mus','var')), mus=[]; end
    if (~exist('sigmas','var')), sigmas=[]; end
    if (~exist('do_normalize','var') || isempty(do_normalize)), do_normalize=false; end
    
    if (do_normalize)
        pdfs = pdfs ./sum(pdfs);
    end
    
    if (isempty(mus)),    mus    = pdf_moment(pdfs, 1, bins, mus, false);   end
    if (isempty(sigmas)), sigmas = pdf_sigmas(pdfs, bins, mus, false);      end
        
    skewness = pdf_moment(pdfs, 3, bins, mus, false);
    skewness = skewness ./ (sigmas .^ 3);
        
end

