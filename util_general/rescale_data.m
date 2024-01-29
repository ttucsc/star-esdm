function [data, minval, sminval] = rescale_data(data, data_scaling, minval)
% [obs_data, hist_data, model_data, futobs_data] = ARRM_V2_rescale_data(obs_data, hist_data, model_data, futobs_data, data_scaling, pwr)
%  
%   Rescales data to log or power scale.  Scales data^data_scaling, then subtracts scaled minval.  That way the
%   distribution can be mirrored to get something close to gaussian, without a gap in the middle.
%   Otherwise for an exampl minval of .1, .1^1/4 -> .5623, so mirroring would leave a gap of 1.26 between the two
%   halves.
%   Reverse the scaling by using unscale_data(...).
%
%   NOTE:  values <= 0 are not scaled or rescaled by these functions.

    flags = data>0;
    if (~exist("minval","var") || isempty(minval))
        minval = nanmin(data);
    end

    if (ischar_s(data_scaling))
        if (strcmp(data_scaling,'linear'))
            return;
        elseif (strcmp(data_scaling,'log'))
            sminval = log(1+minval);
            data(flags)    = log(1+data(flags));
        else
            error("bad scaling type:  %s.  must be linear, log or numeric (for power)", data_scaling);
        end
    elseif (data_scaling == 1)
        return
    else
        pwr = data_scaling;
        sminval = minval^pwr;
        data(flags) = ((data(flags)).^(pwr));
    end

end

