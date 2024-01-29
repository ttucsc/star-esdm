function [probCounts_day, totCount_day] = find_probcounts_day( problines, anoms, hist_edges)
% returns overall counts # of datapoints < each probability line.
%       If pdf was calculated correctly, then there should be .1% of points < .001 probline, 50% of points < .50
%       probline, etc.
%
%   Inputs:
%       problines(yrlen, nprobs) anomaly value for each prob for each day of the year.
%       anoms(ndays)            is column vector of anomaly values,
%                                   where ndays = nyrs * yrlen.
%                             OR
%                               anoms(ndays,nedges) is a histogram of size ndays x nedges (see hist_edges)
%       hist_edges              histogram edges, if anoms is actually a histogram.
%                                   
%   Outputs

%       probCounts_day          probCounts by day of year
%       totCount_day            total count by day of year.
%       

    [yrlen,nprobs] = size(problines);
    probCounts_day  = zeros(yrlen, nprobs);        
    if (~exist('hist_edges','var') || isempty(hist_edges))      % if working with anomalies
        nyrs = numel(anoms(:))/yrlen;
        totCount_day = sum(~isnan(reshape(anoms,yrlen, nyrs)), 2);

        for i=1:nprobs
            amt = repmat(problines(:,i), nyrs, 1);
            keepers = anoms(:) <= amt; 
            probCounts_day(:,i) = nansum(reshape(keepers, yrlen, nyrs),2);
        end    
    else
        error('find_probcounts_day not implemented for histograms yet');
    end

%     fprintf('find_probcounts:\npct ');
%     fprintf('%11.5f ', probCounts/probCounts(end)*100);
%     fprintf('\ncnt ');
%     fprintf('%11.4f ', probCounts);
%     fprintf('\n');
   
end

