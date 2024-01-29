function [packed, nanval, offset, scaling] = ncpack(data, nbytes, flag_out_of_range, offset, scaling)
% packs floating point data to unsigned integer values.
%
%       outdata = round((data - offset) / scaling)
%
%   Inputs:
%       data        data to be packed or unpacked
%       nbytes      size for output.
%                       nbytes s/b 1,2 or 4 for ubyte, ushort or uint
%       flag_out_of_range
%                   boolean.  true:  values out of range will be stored as nanval
%                                       nanval is 255, 65535 or 2^32-1 for nbyters == 1,2 or 4
%                             false: values out of range truncated at 0 or 255, 65535 or 2^32-1
%                       default:  true
%       offset      offset subtracted from each value before packing
%                       default:  min(data)
%       scaling     scale value to bring all data into range for packed size
%                       default:  (max(data)-offset) / out_maxval
%                           where out_maxval is 254,  65534, or 2^32-2  if flag_out_of_range)
%                                   or          255, 65535,  or 2^32-1  if not flagging out of range.
%  
%   Outputs:
%       outdata     packed or unpacked data
%       offset      offset value used.  
%       scaling     scaling value used.

    if (nargin < 3), error("ncpack:  missing arguments"); end
    if (~exist('flag_out_of_range','var') || isempty(flag_out_of_range)), flag_out_of_range = true; end
    if (~exist('scaling','var')           || isempty(scaling)),           scaling = []; end
    if (~exist('offset','var')            || isempty(offset)),            offset = []; end
    
    if (nbytes == 1)
        if (flag_out_of_range), maxval = 254;  nanval = 255; else, maxval = 255; end
    elseif (nbytes == 2) 
        if (flag_out_of_range), maxval = 65534; nanval = 65535; else, maxval = 65535; end
    elseif (nbytes == 4)
        if (flag_out_of_range), maxval = 2^32-2;  nanval = 2^32-1; else, maxval = s^32-1; end
    else
        error('bad nbytes (%d):  must be 1,2 or 4', nbytes);
    end

    if (isempty(offset)), offset = min(data); end
    if (isempty(scaling)), scaling = (max(data-offset)) / maxval; end

    packed = round((data - offset) / scaling);

    if (flag_out_of_range)
        packed(packed > maxval) = nanval;
        packed(packed < 0)  = nanval;
        packed(isnan(data)) = nanval;
    else
        packed(packed > maxval) = maxval;
        packed(packed < 0) = 0;
    end
    

    if (nbytes == 1)
        packed = uint8(packed);
    elseif (nbytes == 2)
        packed = uint16(packed);
    else
        packed = uint32(packed);
    end

end
        