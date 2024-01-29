function problines = find_problines_from_histcounts(probs, histcounts, edges, is_cdf)
% returns counts # of datapoints < each probability line.
%       If pdf was calculated correctly, then there should be .1% of points < .001 probline, 50% of points < .50
%       probline, etc.
%
%       This is a barebones routine, for finding problines from a set of histograms.  
%       Use find_problines(...) for a smarter one which checks for various problems and is safer for working with 
%       pdfs and cdfs where the edges to near zero but not quite all the way.
%
%   Inputs:
%       probs(nprobs,1)          probability values to find lines for
%       histcounts               histogram of datapoints, of size yrlen x nbins
%                                   can also be a pdf or cdf.
%       edges                     bin edges
%       is_cdf                   if present, flags whether histcounts is a cdf.  if false, histcounts can be histogram
%                                   or pdf. 

%   Outputs
%       problines(yrlen,nprobs)  anomaly values (fractional, interpolated from edges) of for reach probability for each 
%                                   day of year.  Size yrlen x nprobs

    nprobs = length(probs);
    [yrlen,nbins] = size(histcounts);
    problines = zeros(yrlen,nprobs);
    
    if (~is_cdf)
        cumhist = cumsum(histcounts,2,'omitnan');
        cdf = cumhist ./ cumhist(:,end);
    end
    pdf = diff([zeros(yrlen,1),cdf]);   % create the pdf from the cdf so we can check for zeros. stuff near one might end up being
    allflags = pdf ~= 0;

    for j=1:yrlen
        flags = allflags(j,:);
        ix=max(2,find(flags,1));
        flags(ix-1)=true;
        problines(j,:) = interp1(cdf(j,flags),edges(flags),probs,'makima',nan);
        bad = isnan(problines(j,:));
        if (any(bad))
            problines(probs(bad)<cdf(j,1))=1;
            problines(probs(bad)>cdf(j,end))=nbins;
        end
    end
end

