function [d_out, na_map, randix] = weighted_select(s, weights, check_NAs, fieldname, seed_or_randix)   % model time series to downscale
%
%   Returns a column vector of weighted-random selection from each row of s.
%   na_map is a boolean vector the same size as d_out, set to true wherever all points in a row s were NAs
%   Inputs
%       s               data to select from.  Data is assumed to be 1 series per column in s
%                           if s is a struct, data is taken from s.(fieldname)  (usually a DA (disaggregated data)
%                           if s is a cell array of structs, data is taken from s{:}.(fieldname) (cell array of DAs)
%                           if s is a cell array of vectors, data is taken from s{:}
%                           if s is a vector or matrix, data is s
%                               for cell arrays and structs, data must be all the same length.
%       weights         weights to use for selecting.  One weight for each column.  
%                           Weights do not need to sum to 1... they will be normalized to sum to 1.
%                           so one can use a vector of all 1's to use even weighting.
%                           To use 1 set as the default, and the rest to fill in NAs, set the default weight to 1, and
%                           the remaining weights to a very small value e.g. (1e-15);
%       check_NAs       boolean.  If true, if an NA is selected, select from one of the non-NA streams instead. [true]
%                           if false, simply pass back NA if an NA is selected.
%                           Set to false if there are a large number of NAs in one data stream and you want to avoid
%                           biasing the results.
%       fieldname       If s is a struct, take data from s.(fieldname).  Ignored if s is not a struct.
%       seed_or_randix  (optional) if present can be either a seed for the random number generator,
%                           or a vector of column indexes to use when selecting
%                           This provides 2 different ways to use the same selections for multiple runs.
%
%   note:  if length(weights)==1, assumes even weighting.  
%           This is the trivial case for when there's only one data series.

            % if s is a struct or cell array, extract the data into a column matrix
    if (exist('fieldname','var') && isstruct(s)) % s is a DA struct.
        d = s.(fieldname);
    elseif (iscell(s))          % s is a cell array of DA structs.  Extract the field from the cell arrays and make a 2-D array from it.
        if (~exist('fieldname','var')), fieldname=[]; end
        d = DAs2mat(s, fieldname);
    else                    %s is just an array of data values.       
        d=s;
    end
    
            % if only 1 weight, if only 1 data set, return it.  Otherwise, assume even weighting
    if (length(weights)==1)
        if (isrow(d))           % single row
            d_out = d(1,:)';
            return;
        elseif (size(d,2)==1)   % single column
            d_out = d;
            return
        else    % matrix, but only 1 weight, so use even weighting.
            weights = ones(1,size(d,2));
        end
    end
        
    [nr,nc]   = size(d);
    
            % are we given a seed, or a matrix of indexes to select?
    if (~exist('seed_or_randix','var') || isempty(seed_or_randix))
        seed='shuffle';
        randix=[];
    elseif (length(seed_or_randix)==1)
        seed=seed_or_randix;
        randix=[];
    else
        seed='shuffle';
        randix=seed_or_randix;
    end
    if (~exist('check_NAs', 'var') || isempty(check_NAs)), check_NAs = true; end 
    
    weights = weights / sum(weights);
    cumwts = cumsum(weights);

    if (isempty(randix))    % find index of column to keep for each row.
        rstream = RandStream('multFibonacci','Seed',seed);
                % note:     creating a private random stream so we don't interfere with the global one.
                % note 2:   Specifying the RNG type explicitly.  Parfor loops use a different default RNG than serial,
                %           so we want to be sure we use the same RNG under all situations.
        r = rstream.rand(nr,1);
        randix = int8((r<=cumwts(1)));
        nwts = length(cumwts);
        for i=2:nwts
            randix(r>cumwts(i-1) & r<=cumwts(i)) = i;
        end
    end
    
%    fprintf("-----weighted_select:  first 50 selects are:  ");
%    fprintf("%d ",randix(1:min(50,length(randix))));
%    fprintf(" -----\n");
    
    
%   d_out = d(:, randix);                           % select one from each row..  Except this is the wrong syntax in matlab...
    d_out = d(sub2ind([nr,nc],(1:nr)', randix));    % select one from each row. This is the way to do it.
                                                    % sub2ind(...) creates linear indices from the (r,c) inputs.
    if (check_NAs)
        % if any points are NAs, select one of the non-NA points instead.
        mynans = find(isnan(d_out));
        nnans = length(mynans);

            % this is brute force.  Slow if there are lots of NAs, but OK if there are only a few.
        for i=1:nnans
            jx = mynans(i);                         % a row with NAs
            keepers=~isnan(d(jx,:));                % find non-NAs in row jx
            if (any(keepers))
                dd = d(jx,keepers);                 % get the valid data from ith row.
                if (length(dd)==1)                  % if only 1 valid, use it.
                    d_out(jx) = dd;
                else
                    ww=weights(keepers);                % otherwise play the weighted random selection game for the
                    cww=cumsum(ww);                 % valid data.
                    cww=cww/cww(end);
                    r=cww(end)*rand();
                    for j=1:length(cww)
                        if (r <= cww(j))
                            d_out(jx) = dd(j);
                            randix(jx) = -j;        %#ok<AGROW> % update randix for this location.  negative flags it as alternate selection.
                            break;
                        end
                    end
                end
            end
        end
    end    
    if (nargout >1)           % only bother to find the na_map if the caller is keeping it.
        na_map = isnan(d_out);
    end
end
