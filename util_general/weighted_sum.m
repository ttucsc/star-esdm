function d_out = weighted_sum(s, weights, fieldname)
%   returns vector of weighted sums of s or of s.(fieldname), if fieldname is specified.
%   Inputs
%       s           cell array of DA structs with data in fieldname to sum, or 
%                       2-D matrix of data to sum.  Data is assumed to be 1 series per column in s
%                       if s is a struct or cell array of structs, data is taken from s.(fieldname)
%       weights     weights to use for summing.  One weight for each column.  
%                       Weights do not need to sum to 1... they will be normalized to sum to 1.
%                       so one can use a vector of all 1's to use even weighting.
%                       If only 1 sum provided, assumes even weighting.
%       fieldname   If s is a struct, take data from s.(fieldname)

%   note:  d_out(j) will be nan if all values in row j are NAs
%   note:  d_out will be properly scaled in any row where there are NAs.
%               NAs are dropped, and remaining points' weights are rescaled to sum to 1.
%               To preferentially select from 1 column, but use remainder in case of NAs, set preferred column to near 1,
%               and remaining columns to very small values (1e-10, e.g.)
%   NOTE:  weights are normalized to 1, so if you mean to have your values scaled to something other than 1, you'll need
%           to do that yourself after calling this routine.
%               
%   note:  if length(weights)==1, then s (or s.(fieldname) is simply returned.  
%           This is the trivial case when there's only one data series.

            % if s is a struct or cell array, extract the data into a column matrix
    if (exist('fieldname','var') && (isstruct(s) || isa(s,'ARRM_V2_DisaggregateSignal'))) % s is a DA object
        d = s.(fieldname);
%        if (isempty(d)), return; end
    elseif (iscell(s))          % s is a cell array of DA structs.  Extract the field from the cell arrays and make a 2-D array from it.
        if (~exist('fieldname','var')), fieldname=[]; end
%        if (isempty(s.(fieldname))), return; end
        d = DAs2mat(s, fieldname);
    else                    %s is just an array of data values.       
        d=s;
    end
    
    sz    = size(d);
    nsets = sz(end);
    sz    = sz(1:end-1);   % drop last dimension;  we'll sum over that dimension
    nd    = length(sz);
    if (nd>3), error('error:  weighted_sum can''t handle data with %d dimensions',nd+1); end
    
            % if only 1 weight, if only 1 data set, return it.  Otherwise, assume even weighting
    if (length(weights)==1)
        if (isrow(d))           % single row
            d_out = d(1,:)';
            return;
        elseif (size(d,2)==1)   % single column or matrix.
            d_out = d;
            return
        else    % matrix, but only 1 weight, so use even weighting.
            weights = ones(nsets,1);
        end
    end
        
            % OK, now do weighted sum of each row.
    d_out = squeeze(zeros([sz,1]));     % add a singleton dimension for zeros so we only get 1 dimension if data is 1D.
    wsum  = squeeze(zeros([sz,1]));     % note:  zeros(n)  gives us an nxn matrix of zeros...
    
    weights = weights(:);   % in case the weights were 2x2 instead of 1x4 or 4x1...
    for i=1:length(weights)
        if (nd==1)
            dd=d(:,i);
        elseif (nd==2)
            dd=d(:,:,i);
        else
            dd=d(:,:,:,i);
        end
        mynans = isnan(dd);
        ww = weights(i) * (~mynans);        % keep weights for non-NAs only
        dd(mynans)=0;                       % add nothing for points which are nan
        wsum = wsum + ww;                   % sum the weights, to scale the results properly wherever some of the points were NAs.
        d_out = d_out + ww .* dd;
    end
    d_out(wsum==0) = nan;       % these were all NAs.  set their output the NA.
    d_out = d_out ./ wsum;      % normalize by remaining weights.
end

