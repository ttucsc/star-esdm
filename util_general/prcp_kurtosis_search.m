function [scale_best, found, kurt_best, sig_best] = prcp_kurtosis_search(prcp, prcp_min, max_scaling)
%   Does a quick-and-dirty search for an approximate best scaling to result in an excess-kurtosis of 0 -- i.e., as close
%   to gaussian as we can get using power scaling.
%   Steps through range 1 to maxscaling in steps of 0.5, finds two points straddling zero kurtosis, and steps through that
%   range steps of .01.  Binary search would be faster...

%     if (~exist("koffset","var") || isempty(koffset)), koffset = 0; end    % koffset removed for now.
    
    pmax = nanmax(prcp);
    binrange = [prcp_min, pmax];
    nedges = 5*pmax;
    prcp = prcp(:);
%   wet_days = prcp>prcp_min;
    big_steps = 1:.5:max_scaling;
    nsteps = length(big_steps);
    kurt = nan(nsteps,1);
    sigs  = nan(nsteps,1);
    
%   [kurt(1), hcnts, edges, sig,      mu] = calc_precip_kurt(prcp, 1.0, binrange, nedges,  false, []);
    [kurt(1), hcnts, edges, sigs(1), ~  ] = calc_precip_kurt(prcp, 1.0, binrange, nedges, false, []);
%   fprintf("scaling:  %6.3f kurt %6.3f mu %10.5f sig %10.5f\n", big_steps(1), kurt(1), mu, sig);        
   
    if (kurt(1) < 0)
        scale_best = 1;
        kurt_best  = kurt(1);
        sig_best   = sigs(1);
        return;
    end
    found = false;
    for i=2:nsteps
        [kurt(i), ~, ~, sigs(i), ~] = calc_precip_kurt(hcnts, big_steps(i), binrange, nedges, true, edges);
%       fprintf("scaling:  %6.3f kurt %6.3f mu %10.5f sig %10.5f\n", big_steps(i), kurt(i), mu, sig);        
        if (kurt(i)        <=0 || abs(kurt(i))         < 1e-3), found = true; break; end
%       if (kurt(i)+koffset<=0 || abs(kurt(i))+koffset < 1e-3), found = true; break; end
    end
    if (~found || abs(kurt(i)) < 1e-4)
        scale_best = big_steps(i);
        kurt_best = kurt(i);
        return;
    end

    small_steps = big_steps(i-1):.01:big_steps(i);
    nsteps = length(small_steps)-1;
    skurt = nan(nsteps+2, 1);
    ssigs = nan(nsteps+2, 1);
    skurt(1) = kurt(i-1);
    skurt(end) = kurt(i);
    ssigs(1)   = sigs(i-1);
    ssigs(end) = sigs(i);
%   fprintf("scaling:  %6.3f kurt %6.3f\n", small_steps(1), skurt(1));        
    for j=2:nsteps
        [skurt(j), ~, ~, ssigs(j), ~] = calc_precip_kurt(hcnts, small_steps(j), binrange, nedges, true, edges);
%       fprintf("scaling:  %6.3f kurt %6.3f mu %10.5f sig %10.5f\n", small_steps(j), skurt(j), mu, sig);        
        if (skurt(j)        <=0 || abs(skurt(j))         < 1e-3), break; end
%       if (skurt(j)+koffset<=0 || abs(skurt(j))+koffset < 1e-3), break; end
    end

            % keep the better of the two points
    if (abs(skurt(j-1)) < abs(skurt(j)))
        jbest = j-1;
    else
        jbest = j;
    end
    scale_best = small_steps(jbest);
    kurt_best = skurt(jbest);
    sig_best  = ssigs(jbest);
end    