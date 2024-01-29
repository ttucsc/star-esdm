function [skew, kurt, hcnts, edges, sig, mu, bins]= calc_precip_kurt_skew(p, pwr, minval, binmax, offset, nedges, frac, lbl, fignum, mnp, lbl2)
% returns kurtosis and (mirrored) histogram of input precip p, scaled by p^(1/scaling).
%   Inputs:
%       p               raw precip data
%       pwr             power to scale the precip by.  Will calculate p^(1`/scaling), so scaling of 2 uses sqrt(p).
%       minval, binmax  range of precip values in data  (usually [0.1, 1200] good for CONUS, [.1, 2000] for global
%                           max precip recorded in 1 day in CONUS is 1070 mm in one day  (Alvin,  TX, 7/25-26, 1979.  Tropical Storm Claudette) 
%                           max precip recorded in 1 day globally is 1825 mm.
%       offset
%       nedges          # of edges to use in scaled histogram
%                               Also good:  binrange=[.1,max(p)], and nedges = 5*max(p).
%       frac            fractional portion of histogrammed precip to use to calculate skew & kurtosis over
%                           if single value, histogram central portion (.5:  use central .25 to .75)
%                           Or can be 2-element vector,  [min_prob, max_prob]
%   optional:
%       lbl, fignum, mnp, lbl2  info to create plot from histogrammed data.
%
%   Outputs:
%       skew            skew of histogram
%       kurt            "excess kurtosis" of distribution
%       hcnts           scaled histogram of input data
%       edges           edges used in scaled histogram
%       sig, mu         std deviation & mean of distribution
%       bins            Bins for distribution hcounts
%
%   NOTE:  

    nedges = nedges + 1-mod(nedges,2);  % make nedges odd.
    brange = rescale_data_1([minval+offset,binmax], pwr, minval, "forward", offset);      % rescale the binmax

    edges = linspace(brange(1), brange(2), nedges);                              % create edges linear in scaled space.
    dx = edges(2)-edges(1);
    bins = edges(1:end-1) + dx/2;%    mid = ceil(length(edges)/2);

    ps = rescale_data_1(p, pwr, minval, "forward", offset);
    if (~strcmpi(string(pwr),"log"))
        ps(ps<0) = nan;
    end
    hcnts = histcounts(ps, edges);
    hc2 = hcnts/sum(hcnts);
    
    if (length(frac) == 2 || frac < 1)
        if (length(frac)==1), frac = [frac/2, 1-frac/2]; end 
        cdf = cumsum(hc2);
        ix1 = max(1, find(cdf >= frac(1),1)-1);
        ix2 = min(length(bins), find(cdf <= frac(2),1,"last")+1);
    else
        ix1=1;
        ix2=length(bins);
        frac = [0,1];
    end
        

    [mu, sig, skew,    kurt] = pdf_stats(hc2(ix1:ix2), bins(ix1:ix2), [], true);

    if (exist("fignum","var") && ~isempty(fignum) && fignum > 0)

        h=figure(fignum);
        if (exist("mnp","var") && ~isempty(mnp))
            subplot(mnp(1),mnp(2),mnp(3));
            cla ;
        end

        my_mu = (mu-edges(1))/dx;
        gsig = sig/dx;
        mygauss = gauss(nedges-1, gsig, my_mu+.5) * (frac(2)-frac(1));
        bar(bins, hc2);
        hold on;
        plot(bins, mygauss,"r-","linewidth",2);
        ix = find(hc2==0,1,"last");
        jx = find(mygauss > .00001, 1);
        xmin = bins(min(ix,jx));
        ix = find(hc2>0,1, "last");
        jx = find(mygauss > .00001, 1, "last");
        xmax = bins(max(ix,jx));
        xlim([xmin, xmax]);
        hold off;
        title(sprintf("%s, pwr=%s, sig=%.3f skew=%.3f kurt=%.3f", lbl, string(pwr), sig, skew, kurt));
        xlabel("scaled prcip (1/pwr)","interpreter","none");
        ylabel("probability");
        drawnow();
        if (strlength(lbl2)>0 && ~stsrcmp(lbl2,"test"))
            figname=sprintf("/Volumes/lacie_1/data/prcp_figs/calc_prcp_kurt_skew_%s.png", fix_filename(lbl2));
            saveas(h, figname);        
        end
    end                           
end

