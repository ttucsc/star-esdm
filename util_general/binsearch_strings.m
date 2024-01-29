function ix = binsearch_strings(names, name)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    ix1 = 1;
    ix2 = length(names);
    mid = floor((ix1+ix2)/2);
    done=false;
    sname=string(name);
    while (~done)
        smid=string(names{mid});
        if (sname==smid)
            ix=mid;
            return;
        elseif (sname > smid)
            ix1=mid;
        else
            ix2 = mid;
        end
        mid=floor((ix1+ix2)/2);
        delta=ix2-ix1;
        if (delta<=1)
            if (sname==string(names{ix1}))
                ix=ix1;
                return;
            elseif (sname==string(names{ix2}))
                ix=ix2;
                return;
            else
                ix=[];
                return;
            end
        end
    end
end

