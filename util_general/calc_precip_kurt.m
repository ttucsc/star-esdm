function [kurt, hcnts, edges, sig, mu, bins, hc2]= calc_precip_kurt(p, scaling, binrange, nedges, ishist, linear_edges)
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

        binrange = sort(binrange);
        myscaling = 1/scaling;
        minval = binrange(1);        
%         binrange = rescale_data(binrange, myscaling, minval);
%         edges = linspace(binrange(1), binrange(2), nedges);
        binrange = rescale_data(binrange, myscaling, minval);
        edges = linspace(binrange(1), binrange(2), nedges);
        dx = edges(2)-edges(1);
        bb = edges(1:end-1)+dx;
        bins = [fliplr(-bb),bb];
        
        if (ishist)         % if data is already histogrammed into a linear histogram, sum it into scaled histogram.
            linear_hist = p;
            linear_bins = (linear_edges(1:end-1) + linear_edges(2:end))/2;
            scaled_bins = rescale_data(linear_bins, myscaling, minval);
                % make output histogram, and sum each linear bin into output histogram.
            hcnts = zeros(1,length(edges)-1);
            for i=1:length(linear_bins)
                if (~isnan(linear_hist(i)))
                    ix = find(bb < scaled_bins(i), 1, 'last');
                    if (isempty(ix))
                        ix = 1;
                        w1 = 1;
                        w2 = 0;
                    elseif (ix == length(bb))
                        w1 = 1;
                        w2 = 0;
                    else
                        w1 = (scaled_bins(i)-bb(ix))/(bb(ix+1)-bb(ix));
                        w2 = 1-w1;
                    end
                    hcnts(ix)   = hcnts(ix)   + w1 * linear_hist(i);
                    if (w2 ~= 0)
                        hcnts(ix+1) = hcnts(ix+1) + w2 * linear_hist(i);
                    end
                end
            end
        else
                        % data isn't histogrammed yet.  scale it and histogram it
                        
            p(p==0) = nan;
%             ps = rescale_data(p, myscaling, minval);
%             ps(ps<binrange(1)) = nan;
            ps = rescale_data_1(p, myscaling, minval,"forward") - 1;
            ps(ps<0) = nan;
            hcnts = histcounts(ps, edges);
        end
        hc2 = [fliplr(hcnts), hcnts];
        hc2 = hc2 / sum(hc2);
%       kde_pdf = ic_kde(hc2, npts, true);
                
%       [mu, sig, skew, kurt] = pdf_stats(hc2, bins, [], true);
        [mu, sig, ~,    kurt] = pdf_stats(hc2, bins, [], true);
        
%         nbins = length(bins);
%         gscale = nbins/(bins(end)-bins(1));
%         gsig = sig * gscale;
%         gequiv = gauss(nbins, gsig);
        
%         subplot(3,4,i);
%         plot(bins, hc2,'c-');
%         hold on;
%         plot(bins, kde_pdf,"b-","linewidth",3);
%         plot(bins, g_equiv, "r-","linewidth",1.5);
%         hold off;
%         title(ttl(i),"fontsize",15);
%         if (i==1)
%             legend("raw", "kde","gaussian");
%         end
%         grid on;
%         drawnow();

end

