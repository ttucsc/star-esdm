function  excess_threshold  = calc_excess_precip_threshold(keep_frac, prob_surf, pr_vals, prcp_min)
% returns estimated threshold by linearly interpolating the fraction to keep to the probability distribution.
% pr_vals are the precip values (1 x n) for each probability given in prob_surf  (yrlen x n)   may be shorter than 365 if
% year is compressed to fewer days.
    yrlen = size(prob_surf,1);
    thresh = zeros(yrlen,1);
    for i=1:yrlen

        if (keep_frac(i)>=1.0)
            thresh(i)= 0;   %  we'll bump these up to prcp_min at the end;
        else
    %       thresh(i) = interp1(1-prob_surf(i,:), pr_vals, 1-keep_frac(i), "linear",0);  % could use interp1 here, except interp 1 doesn't like duplicate values
            ix = find(prob_surf(i,:) < keep_frac(i),1);
            if (isempty(ix) || ix==1)
                thresh(i) = pr_vals(end);
            else

                        % interpolate between pr_vals(ix-1) and pr_vals(ix)
                thresh(i) = pr_vals(ix) + (keep_frac(i)-prob_surf(i,ix))/(prob_surf(i,ix-1)-prob_surf(i,ix))*(pr_vals(ix-1)-pr_vals(ix));
        %        thresh(i) = interp1(prob_surf(i,20:-1:2), pr_vals(20:-1:2), keep_frac(i));
                if (isnan(thresh(i)))
                    oops("thresh(%d) is nan", i);
                end
            end
        end
    end
%   excess_threshold = thresh;
    excess_threshold = climatology(thresh,12,4);
    excess_threshold(keep_frac == 1 | excess_threshold < prcp_min) = prcp_min;
end

