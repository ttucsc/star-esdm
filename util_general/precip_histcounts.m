function [hcnts, edges, bins, sig, kurt, npts, ntrace]= precip_histcounts(prcp, pwr, binrange, nedges)
% returns histogram counts for precip data, rescaled to 1/pwr.
% binrange(1) sets the min value, which is subtracted off before mirroring the histogram to a shape close to normal.
%  Because histogram is mirrored, the mean and skewness will be zero.
%   For most locations, a power of between 2 and 3 should produce a fairly small skewness.  
%       avg for 18 test cities in CONUS & southern Canada is ~2.24
%       avg for closest Livneh gridpoints to the same 18 site sfrom the Livneh data is 2.64.
%   A reasonable value to use is 2.5.
%
%   However, some locations are vastly different, and need ranges up to 10.
%   Model precip data should have trace amounts removed before finding an appropriate scaling value.
%   Use prcp_kurtosis_search(...) after either removing trace precip events at a threshold of, 0.1 mm to 0.254 mm (.01in) 
%   or use trim_excess_prcp(...) to trim precip to appropriate level for corresponding observation data.
%
%   Hcnts is mirrored;  run KDE on the full hcnts, then keep only the second half to create a pdf. 
%   
%   if prcp is a matrix, it works on each column separately, and returns a matrix of histcounts, of size 2*(nedges-1) x ncols
%
%-------------

    invpwr = 1/pwr;
    minval = binrange(1);
    binrange = rescale_data(binrange, invpwr, minval);


    edges = linspace(binrange(1), binrange(2), nedges);
    dx = edges(2)-edges(1);
    bb = edges(1:end-1)+dx;
    bins = [fliplr(-bb),bb];
    nbins = length(bins);
        
        
    if (isrow(prcp)), prcp = prcp'; end  % is 1D and a row vector, make it a column vector.
    [~,nc] = size(prcp);
    
    hcnts  = zeros(nbins,nc);
    sig    = zeros(1,nc);
    kurt   = zeros(1,nc);

    ntrace = nansum(prcp > 0 & prcp < minval);
    prcp(prcp < minval) = nan;
    npts   = sum(~isnan(prcp),1);
    
    ps = rescale_data(prcp, invpwr, minval);
    
    for c = 1:nc
        hc = histcounts(ps(:,c), edges);
        
        hcnts(:,c) = [fliplr(hc), hc];
        
        [~, sig(c), ~, kurt(c)] = pdf_stats(hcnts(:,c), bins, [], true);
        
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

end
 