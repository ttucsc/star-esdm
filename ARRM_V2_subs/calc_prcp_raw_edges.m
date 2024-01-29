function [edges, bins] = calc_prcp_raw_edges(prcp, prcp_min)
% Set up binning for trimming excess precip.  We'll use a nonlinlear (log-scale) binning here with more bins at the low end.  
    
    if (isa(prcp, "single")), fudge = 1e-7;  else, fudge = 1e-16; end  % fudge used to avoid copmuter arithmetic errors.
    
    pwr = log(2000)/log(10);    % pwr is approx. 3.301.
    prcp_max = max(prcp(:));
    medges = 1001;
    pmax = 2000^(1/pwr);    % max daily precip ever recorded is ~1800 mm/day...

    if (prcp_min > 0)
        edge_1 = prcp_min^(1/pwr);           % make 1st bin catch everything up to prcp_min.
        pedges = [0,linspace(edge_1,pmax, medges)]; % gives us up to 1000 bins, covering .1 mm (or pcrp_min) out to 2000 mm (global max is ~1800 mm/day, per Wikipedia...)
    else
        pedges = linspace( 0,pmax, medges);
    end
    edges = pedges.^pwr - fudge;            % subtract a fudge factor in case computer arithmetic puts minimum precip a fraction below the intended prcp_min.
    nedges = find(edges > prcp_max, 1);
    if (isempty(nedges))
        nedges = medges;
    else
        edges = edges(1:nedges);
    end
    bins  = edges(1:(end-1));       % maybe need to offset by half?
    edges(nedges) = inf;            % set the last bin to infinite, so we don't miss some counts.
                                    % *shouldn't* be necessary...we already made sure our edges cover the full range.
end

%     prcp_max_mdl = max(prcp_mdl(:));
%     if (mdl_prcp_min > 0)
%         mdl_edge_1 = mdl_prcp_min^(1/pwr);
%         pedges_mdl = [0,linspace(mdl_edge_1,pmax, nedges)]; 
%     else
%         pedges_mdl = linspace( 0,pmax, nedges);
%     end
%     edges_mdl = pedges_mdl.^pwr - fudge;                    % but make 1st bin start at prcp_min.
%     nedges_mdl = min(nedges, find(edges_mdl > prcp_max_mdl, 1)+1);
%     edges_mdl = edges_mdl(1:nedges_mdl);
%     edges_mdl(nedges_mdl) = inf;
%     bins_mdl  = edges_mdl(1:(end-1));     % maybe need to offset by half?
