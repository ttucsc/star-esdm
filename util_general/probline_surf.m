function [psurf1, nanmap, psurf2, negatives_count] = probline_surf(psurf, bins, sdvals, pdfdim, fwd, is_cdf, nterms, sig_term, norm_sigmas)
%       function to find probability line surface of a pdf surface.
%       Maps the pdf surface to a dimension of fixed probability values, where the the z dimension is temperature,
%       rather than probability.
%       Forward transform:  
%           input surface:  probability surface, size: yrlen X nbins (probability value for day-of-year X temp/precip)
%               row     day-of-year
%               col     temp/precip bin
%               value   probability for given day (each column sums to 1)
%           output surface:  probline surface, size:  yrlen X nprobs (temp/precip value for day-of-year X #-of-std-dev's
%               row     day-of-year
%               col     probability bin, evenly spaced in standard-deviations-from-mean, (sdvals) 
%                           for example:  -7.5 to 7.5 sigmas, step .25
%                           where sigma is the climatology-smoothed std dev for the day
%               value   temp/precip, 
%                           psurf1:  raw (unsmoothed) probline surface
%                           psurf2:  smoothed surface, smoothed with (nearly-ideal) frequency-domain low-pass filter
%
%       Reverse transform
%           input surface:  probline surface (see above)    nbins x nprobs
%           output surface: std. probability surface        nbins x nbins
%               row     day-if-year
%               col     temp/precip bin
%               value   probability for given day and bin.
%
%   Inputs:
%       surf        pdf, cdf or probline surface, of size ndays X [bins or sdvals]
%                       for best forward results, pdf surface should be kde'd histograms, 
%                       which must be normalized to sum to 1 along each row.
%       bins        temp or precip bins
%       sdvals      standard deviations.  Usually [-7:.05:7].  -7 std dev's -> 1.3e-12 probability
%       pdfdim      1 or 2.  If 2, we have to transpose at start and retranspose back on exit
%       fwd         boolean.  true or false.
%                       If fwd,  then translate from pdf/cdf surface to problines surface
%                       if ~fwd, then translate from problines surface to pdf & cdf surface.
%       is_cdf      true if input surface is a CDF surface rather than a pdf surface.  Needed only for forward.
%       nterms      low-pass filter param for climatology calculation (used for fwd only)
%       sig_term    low-pass filter param for climatology calculation (used for fwd only)
%
%   Outputs:
%       psurf1
%       nanmap
%       psurf2
%               If fwd,  then psurf1 is raw probline surf and psurf2 is smoothed probline surf
%               if ~fwd, then psurf1 is pdf and psurf2 is CDF.

    if (~islogical(fwd))
        if (strcmpi(fwd,'fwd') || strcmpi(fwd,'forward'))
            fwd = true;
        elseif (strcmpi(fwd,'rev') || strcmpi(fwd,'reverse') || strcmpi(fwd,'back'))
            fwd = false;
        else
            error('probline_surf:  fwd must be logical or string ''fwd'' or ''back''');
        end
    end
    
    if (pdfdim == 1)
        psurf = psurf';
        do_transpose=true;
    else
        do_transpose=false;
    end
    if (fwd)
        [psurf1, nanmap, psurf2, negatives_count, minval] = prbline_surf(psurf, bins, sdvals, norm_sigmas, is_cdf, nterms, sig_term);
    else
        [psurf1, nanmap, psurf2, negatives_count, minval] = pdfline_surf(psurf, bins, sdvals);
    end
    
    if (do_transpose)
        psurf1 = psurf1';
        nanmap = nanmap';
        psurf2 = psurf2';
    end

%   if (any(negatives_count) > 0)
    if (minval <  -1e-17)
        fprintf(2, 'warning: probline_surf:  negatives:  %4d mival=%.5g **************\n', sum(negatives_count), minval);
    elseif (any(negatives_count > 0))
        fprintf(      'note: probline_surf:  negatives:  %4d mival=%.10g **************\n', sum(negatives_count), minval);
    end
end

% forward.  create psurf and psurf_smooth from a pdf surface.
function [psurf, nanmap, psurf_smooth, negatives_count, minval] = prbline_surf(pdf, bins, sdvals, norm_sigmas, is_cdf, nterms, sig_term)
    [yrlen, nbins] = size(pdf);
    nproblines = length(sdvals);
    negatives_count = 0;    % not used in fwd;  only in reverse
    minval = 0;             % not used in fwd;  only in reverse

    if (~is_cdf)
        pdf(isnan(pdf))=0;
        cdf = cumsum(pdf, 2);
    else
        cdf=pdf;
    end
        % now go back to pdf from cdf, so we can determine where CDF is flat and avoid that region.
        % NOTE:  prepend column of zeros, to avoid shifting the pdf by 1.
    pdf=diff([zeros(yrlen,1),cdf],1,2);
%     pdf = [pdf,zeros(yrlen,1)];
    
    if (nbins ~= length(bins)), error('error:  probline_surface:  bin size mismatch'), end
        
    pvals = normcdf(sdvals);
    psurf=zeros(yrlen,nproblines);
%   okflags = pdf>0;
    okflags = pdf>1e-17;
    for i=1:yrlen
        flags=okflags(i,:) & pdf(i,:) > 1e-17;
        dd = diff(cdf(i,flags));
        if (any(dd < 1e-17))       % 
            flags = okflags(i,:) & pdf(i,:) > 1e-16;
        end
        try
            mybins = bins .*norm_sigmas(i);
            psurf(i,:) = interp1(cdf(i,flags),mybins(flags),pvals,'pchip',nan);
%           psurf(i,:) = interp1(cdf(i,flags),bins(flags),pvals,'makima',nan);
        catch me
            t = getCurrentTask();
            if (~isempty(t))
                workerID = t.ID;            
                fprintf("prbline_surf:  worker %3d, i=%d\n", workerID);
            else
                fprintf("prbline_surf\n");
            end
            rethrow(me);
        end
    end
    
    nanmap = isnan(psurf);
    if (~exist('nterms','var'))
        psurf_smooth=[];
    else
        if (sum(isnan(psurf(:)))>0)
%           psurf=fillmissing(psurf,'pchip');
            psurf=fillmissing(psurf,'linear');
        end
        psurf_smooth = zeros(size(psurf));
        for i=1:nproblines
            if (any(isnan(psurf(:,i))))
                psurf_smooth(:,i) = nan;
            else
                try
                    psurf_smooth(:,i) = climatology(psurf(:,i), nterms, sig_term, yrlen);
                catch
                    oops();
                end
            end
        end
        psurf_smooth(nanmap) = nan;
    end
end

% reverse.  Create pdf surface from probline surface
function [pdf, nanmap, cdf, negatives_count, minval] = pdfline_surf(psurf, bins, sdvals)

    [yrlen, nsdvals] = size(psurf);
    nbins = length(bins);
    
%   psurf = fillmissing(psurf,'pchip');
    psurf = fillmissing(psurf,'linear');

    if (nsdvals ~= length(sdvals)), error('error:  probline_surface:  sdvals size mismatch'), end
    pvals = normcdf(sdvals);
    cdf=zeros(yrlen,nbins);
    okflags = ~isnan(psurf);
    for i=1:yrlen
        flags = okflags(i,:);
        try           
            cdf(i,:) = interp1(psurf(i,flags),pvals(flags),bins,'makima',nan);
        catch
            t = getCurrentTask();
            if (~isempty(t))
                workerID = t.ID;            
                oops(sprintf("pdfline_surf:  worker %3d, i=%d\n", workerID, i));
            else
                oops(sprintf("pdfline_surf\n"));
            end
            rethrow(me);
        end
    end
    
        % normalized to true cdf (0-1)
%   cdfmax = nanmax(cdf,[],2);
    cdfmax = max(cdf,[],2, "omitnan");
    cdf = cdf ./ cdfmax;
    
    if (nargout<2)
        return
    end
    nanmap = isnan(cdf);
        % sanity check.  Create pdf;  if any points are negative, truncate to zero, renormalize, and recreate CDF
    
            % first, fill in nan's with zeros (at bottom) or 1's (at top) of CDF
%     mid=floor(nbins/2);
    if (sum(isnan(cdf(:)))>0)
        cdf(:,1) = 0;
        cdf(:,end) = 1;
%       cdf = fillmissing(cdf,'pchip',2);
        cdf = fillmissing(cdf,'linear',2);
        
%         yrlen = size(cdf,1);
%         for i=1:yrlen
%             ix=find(okmap(i,:),1);
%             cdf(i,1:ix)=0;
%             jx=find(okmap(i,:),1,'last');
%             cdf(i,jx:end)=1;
%         end
%         cdf(:,    1:mid)=fillmissing(cdf(:,    1:mid),'constant',0);
%         cdf(:,mid+1:end)=fillmissing(cdf(:,mid+1:end),'constant',1);
    end
    pdf = diff([zeros(yrlen,1),cdf],1,2);       % pre-pend zeros to avoid shift when differentiating.
    flags = pdf<0;
            % check for inconsistency...negative probabilities.
            % can occur since ideal filter-smoothing can create oscillations, so resulting smoothed CDF may not be
            % monotonically increasing everywhere.
    negatives_count = sum(flags);
    if (negatives_count > 0)
        minval = min(pdf(:));
        pdf(flags) = 0;
    else
        minval = 0;
    end
            % normalize so they sum to 1
%   pdfsum = nansum(pdf,2);
    pdfsum = sum(pdf,2, "omitnan");
    pdf = pdf ./ pdfsum;
            % rework the CDF if we had any negative pdf values.
    if (negatives_count > 0)
        cdf = cumsum(pdf);
        cdf(nanmap) = nan;
    end
    pdf(nanmap) = nan;
end

