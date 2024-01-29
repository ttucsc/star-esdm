function [kurt, hcnts, edges, sig, mu]= calc_kurt(p, scaling, binrange, nedges, npts, ishist, linear_edges)

        binrange = sort(binrange);      % this needs figuring out, Ian!
        myscaling = 1/scaling;
        minval = binrange(1);        
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
            ps = rescale_data(p, myscaling, minval);
            try
                hcnts = histcounts(ps, edges);
            catch
%               oops();
            end
        end
        hc2 = [fliplr(hcnts), hcnts];
        mpts = 2*npts;
        hc2 = hc2 / mpts;
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
        