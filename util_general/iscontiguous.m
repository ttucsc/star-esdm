function yesno = iscontiguous( vals )
    %returns true if elements of vals are consecutive (increasing or decreasing)
    
    if (length(vals)==1)
        yesno = true;
        return;
    end
    difs = vals(2:end) - vals(1:(end-1));
    if (vals(end)>vals(1))        % increasing
        if (any(difs ~= 1))
            yesno = false;
        else
            yesno = true;
        end
%         for i=2:length(vals)
%             if (vals(i)~=vals(i-1)+1)
%                 yesno = false;
%                 return;
%             end
%         yesno = true;
%         end
%        return;
    else
        if (any(difs ~= -1))
            yesno = false;
        else
            yesno = true;
        end
%         for i=2:length(vals)
%             if (vals(i)~=vals(i-1)-1)
%                 yesno = false;
%                 return;
%             end
%         end
%         yesno = true;
%         return;
%         end
    end
    
end

