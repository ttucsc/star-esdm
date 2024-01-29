function mus = pdf_mus(pdfs, bins, do_normalize)
%
%   This returns mean pdfs.
%   Also calculates means, if empty or not provided
%
%   pdfs are along dimension 1 (down each column).

    if (~exist('do_normalize','var') || isempty(do_normalize)), do_normalize=false; end

    if (do_normalize)
        pdfs = pdfs ./nansum(pdfs);
    end
    mus = pdf_moment(pdfs, 1, bins, [], false);
    
end

