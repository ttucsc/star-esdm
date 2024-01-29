function [outdata, offset, scaling] = ncunpack(packed, nbytes, flag_out_of_range, offset, scaling)
% unpacks floating point data to unsigned integer values.
%
%       outdata = data * scaling + offset
%   Inputs:
%       data        data to be packed or unpacked
%       nbytes      size for output.
%                       4 for single, 8 for double
%       flag_out_of_range
%                   boolean.  true:  Out_of_range_value -> NAs
%                       out_of_range_value:  255, 65535, or 2^32-1
%       offset      offset value used in packing
%       scaling     scale value used in packing.
%                       packing:  packed = (orig_data - offset) / scaling;
%
%                       For 8-bit packing w/ out_of_range flagging, use:
%                           offset = min(data);
%                           scaling = (max(data) - offset) / maxval
%                               where maxval is 254, 65534 or 2^32-2
%                               use 255, 65535 or 2^32-1 if not flagging _out_of_range.
%  
%   Outputs:
%       outdata     unpacked data  

    if (nargin < 5), error("ncunpack:  missing arguments"); end
    if (nbytes ~= 4 && nbytes ~= 8), error('ncunpack:  bad nbytes (%d):  must be 4 or 8', nbytes); end
        
    if (nbytes == 4)
        outdata = single(packed) * scaling + offset;
    elseif (nbytes == 8)
        outdata = double(packed) * scaling + offset;
    else
        error('ncpack:  bad nbytes(%d):  must be 4 or 8', nbytes);
    end

    if (flag_out_of_range)
        if (isa(packed,'uint8'))
            nanval = 255;
        elseif (isa(packed,'uint16'))
            nanval = 65535;
        elseif (isa(packed, 'uint32'))
            nanval = 2^32-1;
        else
            error('ncpack:  data must be uint8, uint16 or uint32');
        end
        outdata(packed==nanval) = nan;
    end

end
        