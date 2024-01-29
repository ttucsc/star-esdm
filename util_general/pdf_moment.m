function [m,mus] = pdf_moment(pdfs, nth, bins, mus, do_normalize)
    % returns nth moment about the mean for each column in pdf (works down dimension 1).
    % if pdfs is 2-D array, returns nth moment for each column.
    % 
    %   assumes:  pdf's sum to 1, unless do_normalize is true.
    %
    %   Will calculate means if mus is not provided.
    %
    %   pdfs must be a column if only 1 pdf.
    
    if (~exist('do_normalize','var') || isempty(do_normalize)), do_normalize = false; end
    
    [nr, nc, nsets]=size(pdfs);
    
    if (nr == length(bins))
        nbins = nr;
        yrlen = nc;
    elseif (nc == length(bins))
        nbins = nc;
        yrlen = nr;
        pdfs = pdfs';
    else
        error('error:  pdf_moment:  bin/pdf size mismatch'); 
    end
    
    if (isrow(bins)), bins=bins'; end

    if (do_normalize)
        pdfs = pdfs ./ nansum(pdfs);
    end
    
    if (~exist('mus','var') || isempty(mus))
        mus = nansum(bins .*pdfs);
    end
    
    if (nth == 1)
        mus=squeeze(mus);       % remove singleton dimension if 3D
        m = mus;
        return;
    end
    
    if (size(pdfs,3)>1 && ismatrix(mus)), mus = reshape(mus,1, size(mus,1), size(mus,2)); end
    
    bb = repmat(bins, 1, yrlen,nsets);
    mm = repmat(mus,nbins,1,1);
    x = bb-mm;
    m = nansum((x.^nth).*pdfs);
    mus=squeeze(mus);       % remove singleton dimension if 3D
    m=squeeze(m);       % remove singleton dimension if 3D
end
