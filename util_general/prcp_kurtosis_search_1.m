function [pwr_best, found, kurt_best, sig_best, mu_best, offset, bins, hc2] = prcp_kurtosis_search_1(prcp, prcp_min, min_pwr, lbl, offset, fignum, mnp)
%   Does a quick-and-dirty search for an approximate best scaling to result in an excess-kurtosis of 0 -- i.e., as close
%   to gaussian as we can get using power scaling.
%   Steps through range 1 to maxscaling in steps of 0.5, finds two points straddling zero kurtosis, and steps through that
%   range steps of .01.  Binary search would be faster...
%
% maybe best to to limit this to 1.5 - , or 7?  If outside this range, distribution is likely rather quirky, and maybe
% best to just go with a central number?

%     if (~exist("koffset","var") || isempty(koffset)), koffset = 0; end    % koffset removed for now.
    
        % defaults:
    if (~exist("lbl",   "var")), lbl    = "";  end    
    if (~exist("offset","var")), offset = 0.2; end
    if (~exist("fignum","var")), fignum = [];  end
    if (~exist("mnp",   "var")), mnp    = [];  end
    
    prcp(prcp<prcp_min) = nan;
    pmax = nanmax(prcp);
    binmax = ceil(pmax/100)*100;
    nedges = min(500, binmax);
    prcp = prcp(:);
    kurt = nan(3,1);
    sigs = nan(3,1);
    mus  = nan(3,1);
    bins = cell(3,1);
    hc2  = cell(3,1);
    pwr   = [min_pwr; (1+min_pwr)/2; 1.0];
    step = 0;
    
%     binmax = 2000;
%     nedges = 1000;
    
    end_delta = .001;
    
%     pwr_best = .2;
%     [kkk, ~,    ~,     sss, mmm, bbb, hhh] = calc_precip_kurt_1(prcp, pwr_best, prcp_min, binmax, offset, nedges);
%     kurt_best = kkk;
%     sig_best  = sss;
%     mu_best = mmm;
%     bins = bbb;
%     hc2  = hhh;
%     found = false;
%     return;
    
    
    for i=1:3
%       [kurt,   hcnts, edges, sig,     mu,     bins,    hc2]    = calc_precip_kurt_1(pprp, pwr,    prcp_min, binmax, offset, nedges, fignum)
        [kurt(i), ~,    ~,     sigs(i), mus(i), bins{i}, hc2{i}] = calc_precip_kurt_1(prcp, pwr(i), prcp_min, binmax, offset, nedges);
    end
    del_pwr = pwr(3) - pwr(1);
    if (~isempty(fignum))
        fprintf("step %3d:  pwr:  %6.4f %6.4f %6.4f  del:  %6.4f kurt:  %6.4f %6.4f %6.4f\n", step, pwr(1:3), del_pwr, kurt(1:3));
    end
    
    if (sign(kurt(1)) == sign(kurt(3)))   % problem if best scaling outside range.
        fprintf(2, "\nkurtosis out of range:  %s:  step %3d:  pwr:  %6.4f %6.4f %6.4f  del:  %6.4f kurt:  %6.4f %6.4f %6.4f\n", lbl, step, pwr(1:3), del_pwr, kurt(1:3));
        ix = find(abs(kurt) == min(abs(kurt)),1);
        pwr_best = pwr(ix);
        kurt_best = kurt(ix);
        sig_best  = sigs(ix);
        mu_best = mus(ix);
        bins = bins{ix};
        hc2  = hc2{ix};
        found = false;
        return;
    end
    
    
                % keep the two which straddle 0.
    while (del_pwr > end_delta)                      
        step = step + 1;
        if (sign(kurt(1)) ~= sign(kurt(2)))
            pwr(3)  = pwr(2);
            kurt(3) = kurt(2);
            sigs(3) = sigs(2);
            mus(3)  = mus(2);
            bins{3} = bins{2};
            hc2{3}  = hc2{2};
        else
            pwr(1)  = pwr(2);
            kurt(1) = kurt(2);
            sigs(1) = sigs(2);
            mus(1)  = mus(2);
            bins{1} = bins{2};
            hc2{1}  = hc2{2};
        end
        
        pwr(2) = mean([pwr(1),pwr(3)]);
        del_pwr = pwr(3) - pwr(1);
        
        if (del_pwr < end_delta && ~isempty(fignum))
%           [kurt,   hcnts, edges, sig,     mu,     bins, hc2] = calc_precip_kurt_1(pprp, pwr,    prcp_min, binmax, offset, nedges, fignum)
            [kurt(2), ~,    ~,     sigs(2), mus(2), bins{2}, hc2{2}] = calc_precip_kurt_1(prcp, pwr(2), prcp_min, binmax, offset, nedges, lbl, fignum, mnp);
        else
            [kurt(2), ~,    ~,     sigs(2), mus(2), bins{2}, hc2{2}] = calc_precip_kurt_1(prcp, pwr(2), prcp_min, binmax, offset, nedges);
        end        
        if (~isempty(fignum))
            fprintf("step %3d:  pwr:  %6.4f %6.4f %6.4f  del:  %6.4f kurt:  %6.4f %6.4f %6.4f\n", step, pwr(1:3), del_pwr, kurt(1:3));
        end
    end

    ix = find(abs(kurt) == min(abs(kurt)),1);
    pwr_best = pwr(ix);
    kurt_best  = kurt(ix);
    sig_best   = sigs(ix);
    mu_best    = mus(ix);
    bins = bins{ix};
    hc2  = hc2{ix};
    found = true;
    
%     nedges = length(edges);
%     bins = (edges(1:end-1)+edges(2:end))/2;
%     
%     if (exist("fignum","var") && ~isempty(fignum) && fignum > 0)
% 
%         figure(fignum);
%         clf;
% 
%         gsig = sig_best * length(edges)/bins(end)/2;
%         mygauss = gauss(nedges-1, gsig, mid-1);
%         bar(bins, hc2);
%         hold on;
%         plot(bins, mygauss,"r-","linewidth",2);
%         jx = find(hc2>0,1);
%         kx = find(mygauss > .0001, 1);
%         xmin = bins(min(jx,kx));
%         xlim([xmin, abs(xmin)]);
%         hold off;
%     end
%     
end    