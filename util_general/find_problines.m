function [problines, pdf_xlines, cdf_xlines, linear_flags, nan_flags] = find_problines(pdf, cdf, probs, edges, is_hcounts, pdfdim, cdf_thresh, report_bad, nterms, sig_terms)
%       function to find probability lines on a pdf surface.
%           finds location of probs for each row of a series of pdfs.
%
%       pdf,cdf         pdf, cdf, or histcounts to find probability lines for.  (see is_hcounts and is_cdf)
%                           if histcounts, pdf will be generated from histcounts.
%                           pdf does not need to be equally spaced in x & y...values for each column are given by bins
%               y       day-of-year dimension  (normally;  see pdfdim below)
%               x       value dimension (temp, precip, whatever)
%       probs(1,n)      n probabilities to find in CDF 
%       edges           edges of bins for x-dimension
%       is_hcounts      set this to true if pdf is histogram counts, rather than being a pdf.  
%                           If true, code will normalize each row to area of 1.
%       pdfdim          pdf dimension.  routines works across rows, so if dim is 1, data is rotated, then rotated back.
%       cdf_thresh      threshold for useful part of CDF.  Usually 1e-6 or smaller.
%       report_bad      boolean flag.  If true, prints some error messages to console (for debugging)
%       nterms,         LPF parameters.  If present and not empty, output lines will be low-pass filtered.
%       sig_terms
    % output
%       problines(i,j)  probability for probs(j) on day i
%       pdf_xlines(i,j) pdf value for probs(j) on day i, for drawing problines on pdf surface
%       linear_flags    boolean flags identifying where pchip (spline) interpolation produced NA, and linear interpolation was used.
%       nan_flags       boolean flags identifying where both pchip and linear interpolation produced NAs

    if (~exist('report_bad','var') || isempty(report_bad)), report_bad = false; end
    if (pdfdim == 1)
        pdf = pdf';     % we work with the pdf running along the rows, not down the columns.
        cdf = cdf';
    end

    yrlen = max(size(pdf,1),size(cdf,1));   % should be the same size, unless 1 is empty.

    if (exist('is_hcounts','var') && is_hcounts)
        [pdf, cdf] = make_pdf(pdf, 2);
    elseif (isempty(pdf))
        pdf = diff([zeros(yrlen,1),cdf],2);
    end
%   okflags = ~isnan(pdf) & ~isnan(cdf);
    pdf(isnan(pdf))=0;

    if (isempty(cdf))
        cdf = cumsum(pdf,2);
    end
    
    not_ok = abs(1.0 - cdf(:,end)) > min(cdf_thresh/10, 1e-15); % sanity check, in case we don't have true CDFs at this point.
    if (any(not_ok))                                            % should never happen unless caller said it was a CDF but it wasn't really.
        cdf = cdf ./ cdf(:,end);
    end
    if (length(edges) > size(pdf,2))
        edges = edges(2:end);
    end
    nprobs = length(probs);
    nbins  = length(edges);
    
    problines   =   nan(yrlen,nprobs);
    linear_flags= false(yrlen,nprobs);
    nan_flags   = false(yrlen,nprobs);
    pdf_xlines  =   nan(yrlen,nprobs);
    cdf_xlines  =   nan(yrlen,nprobs);

    flags = ~isnan(cdf) & cdf >= cdf_thresh & cdf <= 1-cdf_thresh;
    flags = flags & pdf > max(1e-15,cdf_thresh/10); % flag where CDF is not monotonically increasing (may be flat)
                                                    % we also check pdf here;  if pdf very near zero, then CDF may
                                                    % be flat after integrating.
                                                    
    cdf(:,1) = 0;                       % fill in any gaps so we can interpolate properly.
    cdf(:,end) = 1;
    cdf = fillmissing(cdf,'linear',2);

    pones = probs == 1.0;
    pzers = probs == 0.0;
    prbflags = ~pones & ~pzers;
    for i=1:yrlen
        cc = max(0,cdf(i,:));       % copy 1 line, and make sure they're all non-negative
        yy = edges;
        try
            xx = interp1(cc(flags(i,:)), yy(flags(i,:)), probs(prbflags), 'makima',nan);     % locate probabilities on this row of CDF  pchip does spline curve to interpolate
            marginal = isnan(xx);
            if (any(marginal))
                xx(marginal) = interp1(cc(flags(i,:)), yy(flags(i,:)), probs(marginal), 'makima',nan);
                still_bad = isnan(xx);
                marginal(still_bad) = false;       % reset those which are still NAs after trying linear interpolation
                linear_flags(i,prbflags) = marginal;
                nan_flags(i,prbflags) = still_bad;
                if (report_bad)
                    fprintf('row %d:  %d linears, %d NAs\n', i, sum(marginal), sum(still_bad));
                end
            end
            problines(i,prbflags) = xx;
            if (any(pones))
                last = find(flags(i,:),1,'last');
                problines(i,pones) = edges(min(nbins,last+1));
            end
            if (any(pzers))
                first = find(flags(i,:),1);
                problines(i,pzers) = edges(max(1,first-1));
            end
        catch
            if (report_bad)
                fprintf(2, 'oops find_problines.  problem with interpolating row %d\n',i);
            end
        end
    end    
            % low-pass filter if nterms & sig_terms provided
    if (exist('nterms','var') && ~isempty(nterms))
        for j=1:nprobs
            problines(:,j) = climatology(problines(:,j), nterms, sig_terms, yrlen);
        end
    end
    
            % check everything is monotonic
    mydifs = diff(problines,1,2);
    ncross = sum(mydifs<0,2,'omitnan');
    if (any(ncross>0))
        fprintf(2, 'warning:  find_problines:  %d probline crossings detected\n', sum(ncross));   % this better not happen!
        problines(mydifs<0) = nan;
    end
    
    
        % find xlines only if calling code is receiving them back.
    if (nargout > 1)
        okflags = ~isnan(pdf);
        for i=1:yrlen
            flags = okflags(i,:);
            try
               pdf_xlines(i, :) = interp1(edges(flags), pdf(i,flags), problines(i,:),'pchip',0);
            catch me
                report_me_error(me,mfilename());
                oops();
            end
        end
    end
    if (nargout > 2)
        okflags = ~isnan(cdf);
        for i=1:yrlen
            flags = okflags(i,:);
            cdf_xlines(i,:) = interp1(edges(flags),      cdf(i,flags), problines(i,:),'pchip',0);
        end
    end
end

