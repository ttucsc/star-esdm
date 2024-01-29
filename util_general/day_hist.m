function [hh, oorange] = day_hist(d, yrlen, edges, do_normalize)
% returns histogram of counts for each day of the year.
%   oorange is count of out-of-range values [nlow, nhigh];
%
%   Inputs:
%       d               data vector
%       yrlen           length of year.  numel(d) must be multiple of yrlen
%       edges           histogram edges.
%       do_normalize    boolean.  If true, histogram is normalized to a PDF, with each row summing to 1.
%                           This adjusts for the presence of NAs on some days but not others.
%
%   Outputs
%       hh              histogram, of size yrlen rows X nbins columns, where nbins = length(edges)-1.
%                           each row is a day's histogram.
%       oorange         2-element counts of out-of-range values.  (1) = # low ( < edges(1),  (2) = # high ( > edges(end)
%
%-------------------------------

    if (~exist('do_normalize','var') || isempty(do_normalize)), do_normalize = false; end
    
    npts = numel(d);
    nyrs = npts/yrlen;
    
    nbins = length(edges)-1;
    hh=zeros(yrlen, nbins); 
    d=reshape(d,yrlen,nyrs);

    for i=1:yrlen                       % histogram is rotated, so nbins rows by 365 columns
        hh(i,:)=histcounts(d(i,:),edges);
    end
    oorange=[sum(d<edges(1)),sum(d>edges(end))];
    
            % this also works, but isn't any faster.
%         x=repmat(1:yrlen,1,nyrs)';
%         hh=histcounts2(x, d, 1:(yrlen+1), edges);  

    if (do_normalize)       % normalize counts so each day has same number of counts.
                            % this adjusts for missing days (NAs)
        sums = nansum(hh,2);
        s = sum(sums)./max(sums);
        hh = hh .* s;
    end
end
