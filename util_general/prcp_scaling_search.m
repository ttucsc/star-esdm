function [pwr_best, found, skew_best, kurt_best, sig_best, mu_best, offset, bins, hcnts] = prcp_scaling_search(prcp, prcp_min, min_pwr, minimizer, frac, offset, lbl, fignum, mnp, fidlog, lbl2)
% finds power that minimizes term "minimizer"  (either "skew" or "kurt") for distribution of pcrp.
%   minimizer   is the name of the parameter (skew or kurtosis) that we're trying to minimize/find zero-crossing for.
%               This search doesn't expand its range if zero crossing outside initial range.  
%               When that happens, it simply returns the best of the range (left, middle or center), with found set to
%               false.
    
        % defaults:
    if (~exist("lbl",   "var")), lbl    = "";  end    
    if (~exist("lbl2",  "var")), lbl2   = "";  end    
    if (~exist("offset","var")), offset = 0.2; end
    if (~exist("fignum","var")), fignum = [];  end
    if (~exist("mnp",   "var")), mnp    = [];  end
    if (~exist("fidlog","var")), fidlog = [];  end
    
    prcp(prcp<prcp_min) = nan;
    pmax = max(prcp);
    binmax = ceil(pmax/100)*100;
    nedges = max(500, binmax);

%           I think we want to use the middle 95% of precips to calculate the scaling.  That would help avoid squirrelly stuff?
    prcp = prcp(:);
%     skew = nan(3,1);
%     kurt = nan(3,1);
%     sigs = nan(3,1);
%     mus  = nan(3,1);
%     bins = cell(3,1);
%     hcnts= cell(3,1);
%     
                % if user selected log scaling, we don't try to minimize.  Just calculate the log
    if (strcmp(string(min_pwr),"log") || strcmp(minimizer,"log"))
        if (offset <= 0), error("error:  prcp_scaling_search:  log scaling, but offset not positive:  %d", offset); end
        [skew_best, kurt_best, hcnts,    ~, sig_best, mu_best, bins] = calc_precip_kurt_skew(prcp, "log", prcp_min, binmax, offset, nedges, frac, lbl, fignum, mnp, lbl2);
        pwr_best = "log";
        found = false;
        return;
    end
    sk = table('Size',[3,7],'VariableTypes', {'double','double','double','double','double','cell', 'cell'}, 'VariableNames',{'skew','kurt','mu','sig','pwr', 'bins','hcnts'});
    
    max_pwr = 1;
    sk.pwr  = [min_pwr; (max_pwr+min_pwr)/2; max_pwr];
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
%       [  skew(i),    kurt(i),    hcnts{i},    ~,   sigs(i),   mus(i),    bins{i}] = calc_precip_kurt_skew(prcp, pwr(i), prcp_min, binmax, offset, nedges);
       [sk.skew(i), sk.kurt(i), sk.hcnts{i},    ~, sk.sig(i), sk.mu(i), sk.bins{i}] = calc_precip_kurt_skew(prcp, sk.pwr(i), prcp_min, binmax, offset, nedges, frac);
    end
    del_pwr = sk.pwr(3) - sk.pwr(1);
    if (~isempty(fignum))
        fprintf("step %4d:  pwr:  %7.4f %7.4f %7.4f  del:  %6.4f skew: %7.4f %7.4f %7.4f   kurt:  %7.4f %7.4f %7.4f\n", step, sk.pwr(1:3), del_pwr, sk.skew(1:3), sk.kurt(1:3));
    end
    
%   if (sign(kurt(1)) == sign(kurt(3)))   % problem if best scaling outside range.
    if (any(isnan(sk.(minimizer))) || sign(sk.(minimizer)(1)) == sign(sk.(minimizer)(3)))   % problem if best scaling outside range.
        fprintf(         2, "\nprcp_scaling:  %s,  step %3d:  pwr:  %6.4f %6.4f %6.4f  del:  %6.4f   skew: %6.4f %6.4f %6.4f   kurt:  %6.4f %6.4f %6.4f  %s out of range\n", lbl, step, sk.pwr(1:3), del_pwr, sk.skew(1:3), sk.kurt(1:3), minimizer);
        if (~isempty(fidlog))
            fprintf(fidlog, "\nprcp_scaling:  %s,  step %3d:  pwr:  %6.4f %6.4f %6.4f  del:  %6.4f   skew: %6.4f %6.4f %6.4f   kurt:  %6.4f %6.4f %6.4f  %s out of range\n", lbl, step, sk.pwr(1:3), del_pwr, sk.skew(1:3), sk.kurt(1:3), minimizer);
        end
%       ix = find(abs(kurt) == min(abs(kurt)),1);
        ix = find(abs(sk.(minimizer)) == min(abs(sk.(minimizer))),1);
        pwr_best = sk.pwr(ix);
        skew_best  = sk.skew(ix);
        kurt_best = sk.kurt(ix);
        sig_best  = sk.sig(ix);
        mu_best = sk.mu(ix);
        bins = sk.bins{ix};
        hcnts  = sk.hcnts{ix};
        found = false;
        calc_precip_kurt_skew(prcp, pwr_best, prcp_min, binmax, offset, nedges, frac, sprintf("SKEW OUT %s", lbl), fignum, mnp, lbl2);
        return;
    end
    
 % thoughts:  minimize skew+kurtosis (if possible);  but may need some intelligence to get started.
 % this assumes we can straddle 0...probably needs a different criteria for the sum in order to decide how to move.
 %  perhaps simulated annealing?
 
                % keep the two which straddle 0.
    while (~any(isnan(sk.(minimizer))) && del_pwr > end_delta)                      
        step = step + 1;
%       if (sign(kurt(1)) ~= sign(kurt(2)))
        if (sign(sk.(minimizer)(1)) ~= sign(sk.(minimizer)(2)))
%             pwr(3)  = pwr(2);
%             skew(3) = skew(2);
%             kurt(3) = kurt(2);
%             sigs(3) = sigs(2);
%             mus(3)  = mus(2);
%             bins{3} = bins{2};
%             hcnts{3}  = hcnts{2};
            sk(3,:) = sk(2,:);
        else
            sk(1,:) = sk(2,:);
%             pwr(1)  = pwr(2);
%             skew(1) = skew(2);
%             kurt(1) = kurt(2);
%             sigs(1) = sigs(2);
%             mus(1)  = mus(2);
%             bins{1} = bins{2};
%             hcnts{1}  = hcnts{2};
        end
        
        sk.pwr(2) = mean([sk.pwr(1),sk.pwr(3)]);
        del_pwr = sk.pwr(3) - sk.pwr(1);
        
        if (del_pwr < end_delta && ~isempty(fignum))
%           [kurt,    hcnts,   edges,  sig,     mu,     bins, hc2] = calc_precip_kurt_1(pprp, pwr,    prcp_min, binmax, offset, nedges, fignum)
            [sk.skew(2), sk.kurt(2), sk.hcnts{2},    ~,  sk.sig(2), sk.mu(2), sk.bins{2}] = calc_precip_kurt_skew(prcp, sk.pwr(2), prcp_min, binmax, offset, nedges, frac, lbl, fignum, mnp, lbl2);
        else
            [sk.skew(2), sk.kurt(2), sk.hcnts{2},    ~,  sk.sig(2), sk.mu(2), sk.bins{2}] = calc_precip_kurt_skew(prcp, sk.pwr(2), prcp_min, binmax, offset, nedges, frac);
        end 


        if (~isempty(fignum))
            fprintf("step %4d:  pwr:  %7.4f %7.4f %7.4f  del:  %6.4f skew: %7.4f %7.4f %7.4f   kurt:  %7.4f %7.4f %7.4f\n", step, sk.pwr(1:3), del_pwr, sk.skew(1:3), sk.kurt(1:3));
        end
    end

%    ix = find(abs(kurt) == min(abs(kurt)),1);
    ix = find(abs(sk.(minimizer)) == min(abs(sk.(minimizer))),1);
    pwr_best   = sk.pwr(ix);
    skew_best  = sk.skew(ix);
    kurt_best  = sk.kurt(ix);
    sig_best   = sk.sig(ix);
    mu_best    = sk.mu(ix);
    bins       = sk.bins{ix};
    hcnts      = sk.hcnts{ix};
    found      = del_pwr < end_delta;

        if (~isempty(fidlog) && fidlog>0)
            fprintf(fidlog, "prcp_scaling:  %s: final %4d:  pwr:  %7.4f   del:  %6.4f skew: %7.4f    kurt:  %7.4f \n", lbl, step, pwr_best, del_pwr, skew_best, kurt_best);
        end   
        if (~isempty(fignum))
%           sk(2,:) = sk(ix,:);
%           fprintf("prcp_scaling:  %s: final %4d:  pwr:  %7.4f %7.4f %7.4f  del:  %6.4f skew: %7.4f %7.4f %7.4f   kurt:  %7.4f %7.4f %7.4f\n", lbl, step, sk.pwr(1:3), del_pwr, sk.skew(1:3), sk.kurt(1:3));
            fprintf("prcp_scaling:  %s: final %4d:  pwr:  %7.4f  del:  %6.4f skew: %7.4f   kurt:  %7.4f\n", lbl, step, pwr_best, del_pwr, skew_best, kurt_best);
        end
    
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
