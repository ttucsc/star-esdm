function [probCounts, totCount, probCounts_day, totCount_day] = find_probcounts(problines, anoms, season_days, hist_edges)
% returns seasonal counts # of datapoints < each probability line.
%       If pdf was calculated correctly, then there should be .1% of points < .001 probline, 50% of points < .50
%       probline, etc.
%
%   Inputs:
%       problines(yrlen, nprobs) anomaly value for each prob for each day of the year.
%       anoms(ndays,nsets)      is column vector (1 input set) or matrix (multiple input sets, with weights),
%                                   where ndays = nyrs * yrlen.
%                            OR
%                               anoms(ndays,nedges,nsets) is a histogram of size ndays x nedges x nsets (see hist_edges)
%
%       season_days             'month' or 'season', 'total', OR cell array of days to use for each group. 
%                                   default:  'total' : if missing or empty, returns total over all days of year.
%                               OR:  cell array with day groupings desired.
%                                       e.g.: {[335:365,1:59]; 60:151; 152:244; 245:335} for seasonal.
%                                       for a compressed 120 day year, monthly would be {1:10; 11:20; ...};
%       hist_edges              histogram edges, if anoms is actually a histogram.  (NOT IMPLEMENTED YET!)  look at find_probcounts_histcounts??? 
%                                   
%   Outputs
%       probCounts              count of # of anomalies <= given probability line.
%                                   to get fractions:  probFractions = probCounts / sum(~isnan(anoms(:));
%       totCount                total count of data points
%       probCounts_day          probCounts by day of year
%       totCount_day            total count by day of year.
%       

    [yrlen,  nprobs] = size(problines);
    if (~exist('season_days','var') || isempty(season_days)), season_days = 'total'; end
    if (ischar_s(season_days))
        if (strncmpi(season_days,'total',5)) 
            season_days = {1:yrlen};
        elseif (strncmpi(season_days,'month',5)) 
            if (yrlen == 360)
                season_days = {1:30; 31:60; 61:90; 91:120; 121:150; 151:180; 181:210; 211:240; 241:270; 271:300; 301:330; 331:360};
            elseif (yrlen == 365)
                season_days = {1:31; 32:59; 60:90; 91:120; 121:151; 152:181; 182:212; 213:243; 244:273; 274:304; 305:334; 335:365};
            else
                season_days = {1:31; 32:59; 60:90; 91:120; 121:151; 152:181; 182:212; 213:243; 244:273; 274:304; 305:334; 335:365};     % standard 365.25.  kludge for now.  uses 365-day.
            end
        elseif (strncmpi(season_days,'seaso',5)) % seasonal  (DJF, MAM, JJA, SON)
            if (yrlen == 360)
                season_days = {[331:360, 1:60]; 61:150; 151:240; 241:330};
            elseif (yrlen == 365)
                season_days = {[335:365, 1:59]; 60:151; 152:244; 245:365};
            else
                season_days = {[335:365, 1:59]; 60:151; 152:244; 245:365};    % standard 365.25.  kludge for now.  uses 365-day.                
            end
        elseif (strncmpi(season_days,'quart',5))    % quarterly  (JFM, AMJ, JAS, OND)
            if (yrlen == 360)
                season_days = {1:90; 91:180; 181:270; 271:360};
            elseif (yrlen == 365)
                season_days = {1:90; 91:181; 182:273; 274:365};
            else
                season_days = {1:90; 91:181; 182:273; 274:365};    % standard 365.25.  kludge for now.  uses 365-day.                
            end
        else
            error("season_days s/b 'total','monthly', 'season', or 'quarterly'");
        end
    end
    
    if (exist('hist_edges','var'))
        [probCounts_day,totCount_day] = find_probcounts_day(problines, anoms, hist_edges);
    else
        [probCounts_day,totCount_day] = find_probcounts_day(problines, anoms);
    end
    
    nseasons = length(season_days);
    probCounts = zeros(nseasons,nprobs); 
    totCount   = zeros(nseasons,1);
    for i=1:nseasons
        for j=1:nprobs
            probCounts(i,j) = sum(probCounts_day(season_days{i},j));
        end
        totCount(i) = sum(totCount_day(season_days{i}));    
    end
end

