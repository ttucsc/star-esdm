function ix = binsearch(vals, val)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    ix1 = 1;
    ix2 = length(vals);
    mid = floor((ix1+ix2)/2);
    done=false;
    while (~done)
        if (val==vals(mid))
            ix=mid;
            return;
        elseif (val > vals(mid))
            ix1=mid;
        else
            ix2 = mid;
        end
        mid=ceil((ix1+ix2)/2);
        delta=ix2-ix1;
        if (delta<=1)
            if (val==vals{ix1})
                ix=ix1;
                return;
            elseif (val==vals{ix2})
                ix=ix2;
                return;
            else
                ix=[];
                return;
            end
        end
    end
end

