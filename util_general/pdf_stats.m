function [mus, sigmas, skewness, kurtosis] = pdf_stats(pdfs, bins, mus, do_normalize)
%   [mus, sigmas, skewness, kurtosis] = pdf_stats(pdfs, bins, mus, do_normalize)
%
%   returns stats for pdfs, using bins.
%   Kurtosis is "excess kurtosis" -- i.e., kurtosis - 3.0;
%   If pdfs is matrix, then it expects 
%
%   Inputs:
%       pdfs        one histogram or pdf, or a 2D array w/ pdfs running down the columns.
%       bins        bins used in histogramming.  Bins should be centered, not the edges array!
%       mus         optional mean(s) to use.  If missing or empty, will calculate mean(s)
%       do_normalize    boolean.  If true, will normalize a histogram into a pdf.  Otherwise assumes input is already a 
%                       pdf that sums to 1.
%
%   Outputs:
%       mus         mean(s)
%       sigmas      std. deviation(s)
%       skewness    3rd moment (skewness)
%       kurtosis    excess kurtosis (i.e., 4th_moment/sigma^4 - 3)
%
%
%   For individual values, you can call one of these directly:
%
%   pdf_moment(...)
%   pdf_sigmas(...)
%   pdf_skewness(...)
%   pdf_kurtosis(...)
%       all of which call pdf_moment(...) to do their work.
%
%--------------------------------------

    
    if (~exist('do_normalize','var') || isempty(do_normalize)), do_normalize = false; end

    if (isrow(pdfs)), pdfs=pdfs'; end
    if (isrow(bins)), bins=bins'; end
    
    if (nargout == 1)
        mus = pdf_moment(pdfs, 1, bins,[],do_normalize);
        return;
    end
    
    if (~exist('mus','var')), mus=[]; end
    
    [sigmas,mus] = pdf_sigmas(pdfs,bins, mus, do_normalize);
    
    if (nargout > 2)
        skewness = pdf_skewness(pdfs, bins, mus, sigmas, do_normalize);
        
        if (nargout > 3)
            kurtosis = pdf_kurtosis(pdfs, bins, mus, sigmas, do_normalize);
        end
    end
    
        
end

