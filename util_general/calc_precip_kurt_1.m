function [kurt, hcnts, edges, sig, mu, bins, hc2]= calc_precip_kurt_1(p, pwr, minval, binmax, offset, nedges, lbl, fignum, mnp)
% returns kurtosis and (mirrored) histogram of input precip p, scaled by p^(1/scaling).
%   Inputs:
%       p               either raw precip data, or already (linearly) histogrammed data, histogrammed using linear_edges
%       scaling         power to scale the precip by.  Will calculate p^(1`/scaling), so scaling of 2 uses sqrt(p).
%       binrange        range of precip values in data  (usually [0.1, 1200] good for CONUS, [.1, 2000] for global
%                           max precip recorded in 1 day in CONUS is 1070 mm in one day  (Alvin,  TX, 7/25-26, 1979.  Tropical Storm Claudette) 
%                           max precip recorded in 1 day globally is 1825 mm.
%       nedges          # of edges to use in scaled histogram
%                               Also good:  binrange=[.1,max(p)], and nedges = 5*max(p).
%       ishist          boolean.  set to true if p is already histogrammed, using linear_edges
%       linear_edges    edges already used for binning p into a histogram.
%       binstep         stepsize between bins for binning
%
%   Outputs:
%       kurt            "excess kurtosis" of distribution
%       hcnts           scaled histogram of input data
%       edges           edges used in scaled histogram
%       sig, mu         std deviation & mean of mirrored distribution
%       bins            Bins for distribution hc2  (mirrored, halfway between edges)
%       hc2             mirrored hcnts, normalized to 1, so is ready for running KDE to estimate the underlying pdf.
%
%   NOTE:  

    nedges = nedges + 1-mod(nedges,2);  % make nedges odd.
    bmax = rescale_data_1(binmax, pwr, minval, "forward", offset);      % rescale the binmax
    edges = linspace(-bmax, bmax, nedges);                              % create edges linear in scaled space.
    dx = edges(2)-edges(1);
    bins = edges(1:end-1) + dx/2;
    mid = ceil(length(edges)/2);

    ps = rescale_data_1(p, pwr, minval, "forward", offset);
    ps(ps<0) = nan;
    hcnts = histcounts(ps, edges(mid:end));

    hc2 = [fliplr(hcnts), hcnts];   % make it a double-tailed pdf.
    hc2 = hc2 / sum(hc2);

    [mu, sig, ~,    kurt] = pdf_stats(hc2, bins, [], true);

    if (exist("fignum","var") && ~isempty(fignum) && fignum > 0)

        figure(fignum);
        if (exist("mnp","var") && ~isempty(mnp))
            subplot(mnp(1),mnp(2),mnp(3));
            cla ;
        end

        gsig = sig * length(edges)/bins(end)/2;
        mygauss = gauss(nedges-1, gsig, mid-.5);
        bar(bins, hc2);
        hold on;
        plot(bins, mygauss,"r-","linewidth",2);
        ix = find(hc2>0,1);
        jx = find(mygauss > .0001, 1);
        xmin = bins(min(ix,jx));
        xlim([xmin, abs(xmin)]);
        hold off;
        title(sprintf("%s, pwr=%.3f, sig=%.3f", lbl, pwr, sig));
        xlabel("scaled prcip (1/pwr)","interpreter","none");
        ylabel("probability");
    end                           
end

