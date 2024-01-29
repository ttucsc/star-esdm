function d2 = rescale_data_1(data, pwr, minval, direction, offset)
%   Simple shifted power scaling
%
%       pwr can be a prcp_scaling struct w/ fields as given below in the if(...) block
%       And when pwr is a struct, minval should be the direction, "forward" or "reverse"

    % Require user to pass in both direction and offset.  FYI offset of 0.2 is good.
%     if (~exist("direction","var") || isempty(direction)), direction = "forward"; end
%     if (~exist("offset","var") || isempty(offset))
%         offset=.2;      % this seems to be a good number in general.  Maybe want to adjust it based on shape of distribution?
%     end

    if (isstruct(pwr))
        prcp_scaling = pwr;
        if (isstring(minval) || ischar(minval))
            direction    = minval;
        end
        pwr    = prcp_scaling.prcp_scaling;
        minval = prcp_scaling.prcp_min;
        offset = prcp_scaling.prcp_offset;
        
    end
    if (isa(data,"single"))
        off_fudge = offset - 1e-7;
    else
        off_fudge = offset - 1e-12;
    end
    if (strcmp(direction, "forward"))
        d2 = data + offset - minval;        % shift so minval is at offset (.2, usually)
        d2(d2<off_fudge) = nan;             % map anything < minval to nan.  fudge factor to account for arithmetic rounding differences.
        if (isstring(pwr) && strcmp(pwr,"log"))
            d2 = log(d2);
        else
            d2 = d2 .^ pwr;             % scale by power
            d2 = d2 - (offset)^pwr;                % shift result so minval is now at 0.
            d2(d2<0) = 0;                       % and fix any that actually went below 0.
        end
    elseif (strcmp(direction,"reverse"))
        if (isstring(pwr) && strcmp(pwr,"log"))
            d2 = exp(data);
        else        
            d2 = data + (offset^pwr);                % shift back to unscale
            d2 = d2 .^ (1/pwr);              % reverse scaling
        end
        d2 = d2 - offset + minval;       % shift back to original location.  Again, making sure we don't lost data to computer math precions problems.
    end
end
