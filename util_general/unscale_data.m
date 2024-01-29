function data = unscale_data(data, data_scaling, minval)
%This function does reverse of function rescale_data.
%   Scales data^1/data_scaling.
%   NOTE:  must use same minval as used in rescale_data, because rescale_data subtracts minval after scaling so
%   mirroring the distribution doesn't have a gap in the middle.
%
%   NOTE:  Values <=0 are not scaled.

    ok = data>0;
    if (ischar_s(data_scaling))
        if (strcmp(data_scaling,'linear'))
            return;
        elseif (strcmp(data_scaling,'log'))
            sminval = log(1+minval);
            data(ok) = exp(data(ok)+sminval)-1;
        else
            error("bad scaling type:  %s.  must be linear, log or numeric (for power)", data_scaling);
        end
    else
    pwr = data_scaling;
    sminval = minval ^ pwr;
    data(ok) = (data(ok)+sminval).^ (1/pwr);
    end
end
