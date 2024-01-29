function [pdf, cdf] = make_pdf(hCounts, dim)
%       returns pdf and cdf for histogram in hcounts, along dimension dim
%       note:  pdf(i,:), cdf(i,:) set to uniform if sum(pdf(i,:)) is zero -- no data for given day.

    if (dim == 1)
        hCounts = hCounts';
        do_transpose = true;
    else
        do_transpose = false;
    end
    
    nx = size(hCounts, 2);
                % make it a pdf
    hCounts(isnan(hCounts))=0;          
    s = sum(hCounts,2);
    
    zz = s==0;            % if any rows all zero, then set to nans.
    if (any(zz))
        o = nan(size(hCounts));
        hCounts(zz,:) = o(zz,:);
        s(zz) = nx;
    end
    pdf = hCounts ./ s;
    
    cdf = cumsum(pdf, 2);

    if (do_transpose)
        pdf = pdf';
        cdf = cdf';
    end
    
end

