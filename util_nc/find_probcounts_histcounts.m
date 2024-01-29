function [probCounts, totCount] = find_probcounts_histcounts( problines, probs, histcounts, bins, anoms, weights)
% function [probCounts, totCount] = find_probcounts_histcounts( problines, probs, histcounts, bins, anoms, weights)
%
%   This probably isn't tested yet, Ian!
%
% returns counts # of datapoints < each probability line.
%       If pdf was calculated correctly, then there should be .1% of points < .001 probline, 50% of points < .50
%       probline, etc.
%
%   Inputs:
%       probs(nprobs,1)          probability values to do counts for
%       problines(yrlen, nprobs) anomaly value for each prob for each day of the year.
%                                   if empty, then problines are calculated from histcounts
%       histcounts               histogram of datapoints, of size yrlen x nbins
%       bins                     bin centers
%       anoms(ndays,nsets)       is column vector (1 input set) or matrix (multiple input sets, with weights),
%                                   where ndays = nyrs * yrlen.
%       weights                  weighting for each set (if missing or empty, weights set to all ones)
%   Outputs
%       probCounts               count of # of anomalies <= given probability line.
%                                   to get fractions:  probFractions = probCounts / sum(~isnan(anoms(:));

    nprobs = length(probs);
    [yrlen,nbins] = size(histcounts);
    probCounts = zeros(nprobs,1);

    nbins = length(bins);
    cumhist = cumsum(hcounts,2);
    cdf = cumhist ./ cumhist(:,end);
    nprobs = length(probs);
    binix = zeros(1,length(probs));
    if (probs(end)==1)
        nprobs=nprobs-1;
        binix(end)=nbins;
    end
    allflags = histcounts ~= 0;

    for j=1:yrlen
        flags = allflags(j,:);
        ix=max(2,find(flags,1));
        flags(ix-1)=true;
        pvals = interp1(cdf(j,flags),bins(flags),probs,'makima',0);            
        for i=1:nprobs;  binix(i)=find(bins<pvals(i),1,'last'); end 
        pcounts = cumhist(j, binix)';
        fprintf('%3d ',j);
        fprintf('%9.3f ', pcounts);
        fprintf('\n');
        fprintf('%3d ',j);
        fprintf('%9.3f', probCounts);
        fprintf('\n');
        probCounts = probCounts + pcounts;
    end
    totCount = sum(cumhist(:,end));


%         for j=1:yrlen
%             pcounts = interp1(bins,cumhist(j,:),problines(j,:),'makima',0);
%             pcounts(palls)=cumhist(j,end);
%             probCounts = probCounts + pcounts';
%         end
%         totcount = probCounts(end);
    fprintf('find_probcounts:\nend ');
    fprintf('%9.5f ', probCounts/probCounts(end)*100);
    fprintf('\n');
    fprintf('%9.3f', probCounts);
    fprintf('\n');
    
end

